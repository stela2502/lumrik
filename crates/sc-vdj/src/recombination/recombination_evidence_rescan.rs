use super::{
    local_alignment, refresh_recombination_from_observed, ConstantRegionEvidence, Recombination,
    RecombinationId,
};
use crate::index::{reverse_complement, Chain, SegmentId, SegmentKind, VdjIndex};
use crate::runner::BamIdentityResolver;
use anyhow::{Context, Result};
use fast_tag_mapper::{FastLocusMapper, FeatureEntry, MapStatus};
use int_to_str::IntToStr;
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use scdata::CellHash;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::Path;

const CDR3_UPSTREAM_BASES: usize = 32;
const J_BAIT_BASES: usize = 24;
const RESCAN_EVIDENCE_BATCH_SIZE: usize = 200_000;
const RESCAN_PROGRESS_EVERY_BAM_RECORDS: usize = 100_000;
const MIN_JUNCTION_OVERLAP: usize = 8;
const MIN_REFINE_DEPTH: u32 = 3;

#[derive(Debug, Clone, Copy)]
struct ReceptorTarget {
    call_group: usize,
    call_index: usize,
}

#[derive(Debug, Default, Clone, Copy)]
struct CandidateSupport {
    fragments: u32,
    spanning_reads: u32,
}

#[derive(Debug, Clone)]
pub struct RecombinationRescanCandidate {
    pub segment: SegmentId,
    pub fragments: u32,
    pub same_read: u32,
}

#[derive(Debug, Clone)]
pub struct RecombinationRescanCall {
    pub cell_id: u64,
    pub stable_id: RecombinationId,
    pub bait: Vec<u8>,
    pub candidates: Vec<RecombinationRescanCandidate>,
    pub rescued_segment: Option<SegmentId>,
}

#[derive(Debug, Clone, Default)]
pub struct RecombinationEvidenceRescanReport {
    pub bam_records_scanned: usize,
    pub wanted_cell_records: usize,
    pub receptor_hit_records: usize,
    pub constant_hit_records: usize,
    pub linked_fragments: usize,
    pub junction_support_reads: usize,
    pub junction_spanning_reads: usize,
    pub junction_conflicting_reads: usize,
    pub junction_refined_calls: usize,
    pub junction_refined_bases: usize,
    pub rescued: usize,
    pub calls: Vec<RecombinationRescanCall>,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct RecombinationEvidenceRescanProgress {
    pub bam_records_scanned: usize,
    pub wanted_cell_records: usize,
    pub batches_completed: usize,
    pub receptor_hit_records: usize,
    pub constant_hit_records: usize,
    pub linked_fragments: usize,
    pub junction_support_reads: usize,
    pub junction_spanning_reads: usize,
    pub junction_conflicting_reads: usize,
}

#[derive(Debug, Default)]
struct FragmentHits {
    receptors: HashSet<(usize, usize, Chain)>,
    constants: HashSet<(SegmentId, Chain)>,
    spanning_reads: HashMap<(usize, usize, SegmentId), u32>,
}

#[derive(Debug)]
struct RescanRecord {
    fragment_id: u64,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
}

#[derive(Debug, Clone, Default)]
struct JunctionPileup {
    reference_start: Option<usize>,
    counts: Vec<[u32; 4]>,
    quality_sums: Vec<[u32; 4]>,
    support_reads: u32,
    spanning_reads: u32,
    conflicting_reads: u32,
}

#[derive(Debug, Default)]
struct CellRescanEvidence {
    support: HashMap<(usize, usize, SegmentId), CandidateSupport>,
    rediscovery_reads: HashMap<(usize, usize), u32>,
    junction_pileups: HashMap<(usize, usize), JunctionPileup>,
    receptor_hit_records: usize,
    constant_hit_records: usize,
    linked_fragments: usize,
    junction_support_reads: usize,
    junction_spanning_reads: usize,
    junction_conflicting_reads: usize,
}

#[derive(Debug, Default)]
struct RescanEvidenceVdj {
    cells: CellHash<CellRescanEvidence>,
    receptor_hit_records: usize,
    constant_hit_records: usize,
    linked_fragments: usize,
    junction_support_reads: usize,
    junction_spanning_reads: usize,
    junction_conflicting_reads: usize,
}

#[derive(Debug, Default)]
struct RescanTotals {
    support: HashMap<(usize, usize, SegmentId), CandidateSupport>,
    rediscovery_reads: HashMap<(usize, usize), u32>,
    junction_pileups: HashMap<(usize, usize), JunctionPileup>,
    receptor_hit_records: usize,
    constant_hit_records: usize,
    linked_fragments: usize,
    junction_support_reads: usize,
    junction_spanning_reads: usize,
    junction_conflicting_reads: usize,
}

impl RescanEvidenceVdj {
    fn consume_batch(
        &mut self,
        batch: Vec<(u64, RescanRecord)>,
        receptor_mapper: &FastLocusMapper,
        receptor_target_by_cell: &HashMap<(usize, u64), ReceptorTarget>,
        constant_mapper: &FastLocusMapper,
        constant_targets: &[SegmentId],
        index: &VdjIndex,
        calls: &[(u64, Vec<Recombination>)],
    ) {
        let mut by_cell = HashMap::<u64, Vec<RescanRecord>>::new();
        for (cell_id, record) in batch {
            by_cell.entry(cell_id).or_default().push(record);
        }

        let deltas: Vec<_> = by_cell
            .into_par_iter()
            .map(|(cell_id, records)| {
                let mut fragments = HashMap::<u64, FragmentHits>::new();
                let mut delta = CellRescanEvidence::default();

                for record in records {
                    let receptor_hit = map_receptor_oriented(
                        receptor_mapper,
                        cell_id,
                        &record.sequence,
                    );
                    let receptor = receptor_hit
                        .as_ref()
                        .and_then(|(feature_index, _, _)| {
                            receptor_target_by_cell
                                .get(&(*feature_index, cell_id))
                                .copied()
                        })
                        .map(|target| {
                            let chain = calls[target.call_group].1[target.call_index].chain;
                            (target.call_group, target.call_index, chain)
                        });
                    let constant = map_constant(constant_mapper, &record.sequence)
                        .and_then(|feature_index| constant_targets.get(feature_index).copied())
                        .and_then(|segment| {
                            index.segment(segment).map(|entry| (segment, entry.chain))
                        });

                    if let Some((call_group, call_index, _)) = receptor {
                        delta.receptor_hit_records = delta.receptor_hit_records.saturating_add(1);
                        let reads = delta
                            .rediscovery_reads
                            .entry((call_group, call_index))
                            .or_default();
                        *reads = reads.saturating_add(1);

                        if let Some((_, oriented_sequence, reverse)) = receptor_hit.as_ref() {
                            let oriented_qualities = if *reverse {
                                record.qualities.iter().rev().copied().collect::<Vec<_>>()
                            } else {
                                record.qualities.clone()
                            };
                            let call = &calls[call_group].1[call_index];
                            if let Some(read_evidence) = junction_read_evidence(
                                oriented_sequence,
                                &oriented_qualities,
                                call,
                                index,
                            ) {
                                delta.junction_support_reads = delta
                                    .junction_support_reads
                                    .saturating_add(1);
                                if read_evidence.spans_junction {
                                    delta.junction_spanning_reads = delta
                                        .junction_spanning_reads
                                        .saturating_add(1);
                                }
                                if read_evidence.conflicts {
                                    delta.junction_conflicting_reads = delta
                                        .junction_conflicting_reads
                                        .saturating_add(1);
                                }
                                merge_junction_read(
                                    delta
                                        .junction_pileups
                                        .entry((call_group, call_index))
                                        .or_default(),
                                    read_evidence,
                                );
                            }
                        }
                    }
                    if constant.is_some() {
                        delta.constant_hit_records = delta.constant_hit_records.saturating_add(1);
                    }
                    if receptor.is_none() && constant.is_none() {
                        continue;
                    }

                    let fragment = fragments.entry(record.fragment_id).or_default();
                    if let Some(receptor) = receptor {
                        fragment.receptors.insert(receptor);
                    }
                    if let Some(constant) = constant {
                        fragment.constants.insert(constant);
                    }
                    if let (
                        Some((call_group, call_index, receptor_chain)),
                        Some((segment, constant_chain)),
                    ) = (receptor, constant)
                    {
                        if receptor_chain == constant_chain {
                            let reads = fragment
                                .spanning_reads
                                .entry((call_group, call_index, segment))
                                .or_default();
                            *reads = reads.saturating_add(1);
                        }
                    }
                }

                for fragment in fragments.values() {
                    if !fragment.receptors.is_empty() && !fragment.constants.is_empty() {
                        delta.linked_fragments = delta.linked_fragments.saturating_add(1);
                    }
                    for &(call_group, call_index, receptor_chain) in &fragment.receptors {
                        for &(constant, constant_chain) in &fragment.constants {
                            if receptor_chain != constant_chain {
                                continue;
                            }
                            let candidate = delta
                                .support
                                .entry((call_group, call_index, constant))
                                .or_default();
                            candidate.fragments = candidate.fragments.saturating_add(1);
                            candidate.spanning_reads = candidate.spanning_reads.saturating_add(
                                fragment
                                    .spanning_reads
                                    .get(&(call_group, call_index, constant))
                                    .copied()
                                    .unwrap_or(0),
                            );
                        }
                    }
                }

                (cell_id, delta)
            })
            .collect();

        for (cell_id, mut delta) in deltas {
            self.receptor_hit_records = self
                .receptor_hit_records
                .saturating_add(delta.receptor_hit_records);
            self.constant_hit_records = self
                .constant_hit_records
                .saturating_add(delta.constant_hit_records);
            self.linked_fragments = self
                .linked_fragments
                .saturating_add(delta.linked_fragments);
            self.junction_support_reads = self
                .junction_support_reads
                .saturating_add(delta.junction_support_reads);
            self.junction_spanning_reads = self
                .junction_spanning_reads
                .saturating_add(delta.junction_spanning_reads);
            self.junction_conflicting_reads = self
                .junction_conflicting_reads
                .saturating_add(delta.junction_conflicting_reads);

            let cell = self.cells.entry_cell(cell_id).or_default();
            cell.receptor_hit_records = cell
                .receptor_hit_records
                .saturating_add(delta.receptor_hit_records);
            cell.constant_hit_records = cell
                .constant_hit_records
                .saturating_add(delta.constant_hit_records);
            cell.linked_fragments = cell.linked_fragments.saturating_add(delta.linked_fragments);
            cell.junction_support_reads = cell
                .junction_support_reads
                .saturating_add(delta.junction_support_reads);
            cell.junction_spanning_reads = cell
                .junction_spanning_reads
                .saturating_add(delta.junction_spanning_reads);
            cell.junction_conflicting_reads = cell
                .junction_conflicting_reads
                .saturating_add(delta.junction_conflicting_reads);

            for (key, pileup) in delta.junction_pileups.drain() {
                merge_junction_pileup(cell.junction_pileups.entry(key).or_default(), pileup);
            }
            for (key, count) in delta.rediscovery_reads.drain() {
                let total = cell.rediscovery_reads.entry(key).or_default();
                *total = total.saturating_add(count);
            }
            for (key, counts) in delta.support.drain() {
                let total = cell.support.entry(key).or_default();
                total.fragments = total.fragments.saturating_add(counts.fragments);
                total.spanning_reads = total.spanning_reads.saturating_add(counts.spanning_reads);
            }
        }
    }

    fn progress_counts(&self) -> (usize, usize, usize, usize, usize, usize) {
        (
            self.receptor_hit_records,
            self.constant_hit_records,
            self.linked_fragments,
            self.junction_support_reads,
            self.junction_spanning_reads,
            self.junction_conflicting_reads,
        )
    }

    fn into_totals(self) -> RescanTotals {
        let mut totals = RescanTotals::default();
        for bucket in self.cells.into_iter() {
            for (_, cell) in bucket {
                totals.receptor_hit_records = totals
                    .receptor_hit_records
                    .saturating_add(cell.receptor_hit_records);
                totals.constant_hit_records = totals
                    .constant_hit_records
                    .saturating_add(cell.constant_hit_records);
                totals.linked_fragments = totals
                    .linked_fragments
                    .saturating_add(cell.linked_fragments);
                totals.junction_support_reads = totals
                    .junction_support_reads
                    .saturating_add(cell.junction_support_reads);
                totals.junction_spanning_reads = totals
                    .junction_spanning_reads
                    .saturating_add(cell.junction_spanning_reads);
                totals.junction_conflicting_reads = totals
                    .junction_conflicting_reads
                    .saturating_add(cell.junction_conflicting_reads);

                for (key, pileup) in cell.junction_pileups {
                    merge_junction_pileup(totals.junction_pileups.entry(key).or_default(), pileup);
                }
                for (key, count) in cell.rediscovery_reads {
                    let total = totals.rediscovery_reads.entry(key).or_default();
                    *total = total.saturating_add(count);
                }
                for (key, counts) in cell.support {
                    let total = totals.support.entry(key).or_default();
                    total.fragments = total.fragments.saturating_add(counts.fragments);
                    total.spanning_reads =
                        total.spanning_reads.saturating_add(counts.spanning_reads);
                }
            }
        }
        totals
    }
}

pub(crate) fn rescue_missing_constants_from_bam<P: AsRef<Path>, R: BamIdentityResolver>(
    path: P,
    resolver: &R,
    index: &VdjIndex,
    calls: &mut [(u64, Vec<Recombination>)],
    threads: usize,
) -> Result<usize> {
    Ok(
        rescue_missing_constants_from_bam_with_report(path, resolver, index, calls, threads)?
            .rescued,
    )
}

pub(crate) fn rescue_missing_constants_from_bam_with_report<
    P: AsRef<Path>,
    R: BamIdentityResolver,
>(
    path: P,
    resolver: &R,
    index: &VdjIndex,
    calls: &mut [(u64, Vec<Recombination>)],
    threads: usize,
) -> Result<RecombinationEvidenceRescanReport> {
    rescue_missing_constants_from_bam_with_report_and_progress(
        path,
        resolver,
        index,
        calls,
        threads,
        |_| {},
    )
}

pub(crate) fn rescue_missing_constants_from_bam_with_report_and_progress<P, R, F>(
    path: P,
    resolver: &R,
    index: &VdjIndex,
    calls: &mut [(u64, Vec<Recombination>)],
    threads: usize,
    mut progress: F,
) -> Result<RecombinationEvidenceRescanReport>
where
    P: AsRef<Path>,
    R: BamIdentityResolver,
    F: FnMut(RecombinationEvidenceRescanProgress),
{
    let mut report = RecombinationEvidenceRescanReport::default();
    let mut batches_completed = 0usize;
    let mut receptor_mapper = FastLocusMapper::new().with_min_hits(3);
    let mut receptor_targets = Vec::<ReceptorTarget>::new();
    let mut receptor_target_by_cell = HashMap::<(usize, u64), ReceptorTarget>::new();
    let mut receptor_feature_by_bait = HashMap::<Vec<u8>, usize>::new();
    let mut wanted_cells = HashSet::<u64>::new();
    let mut wanted_chains = BTreeSet::new();

    for (call_group, (cell_id, recombinations)) in calls.iter().enumerate() {
        for (call_index, recombination) in recombinations.iter().enumerate() {
            let Some(bait) = receptor_bait(index, recombination) else {
                continue;
            };
            let feature_index = if let Some(&feature_index) = receptor_feature_by_bait.get(&bait) {
                receptor_mapper.add_locus_alias(feature_index, *cell_id);
                feature_index
            } else {
                let feature_index = receptor_mapper.add_feature(
                    *cell_id,
                    &bait,
                    FeatureEntry::new(
                        receptor_feature_by_bait.len() as u64 + 1,
                        recombination.stable_id.to_string(),
                        "vdj_recombination_bait",
                    ),
                );
                receptor_feature_by_bait.insert(bait.clone(), feature_index);
                feature_index
            };
            let target = ReceptorTarget {
                call_group,
                call_index,
            };
            receptor_target_by_cell.insert((feature_index, *cell_id), target);
            receptor_targets.push(target);
            report.calls.push(RecombinationRescanCall {
                cell_id: *cell_id,
                stable_id: recombination.stable_id.clone(),
                bait,
                candidates: Vec::new(),
                rescued_segment: None,
            });
            wanted_cells.insert(*cell_id);
            wanted_chains.insert(recombination.chain);
        }
    }

    if receptor_targets.is_empty() {
        return Ok(report);
    }

    let mut constant_mapper = FastLocusMapper::new().with_min_hits(4);
    let mut constant_targets = Vec::<SegmentId>::new();
    let mut seen_constants = BTreeSet::<SegmentId>::new();
    for chain in wanted_chains {
        for segment in index.segments_for(chain, SegmentKind::C) {
            if !seen_constants.insert(segment.id) {
                continue;
            }
            let feature_index = constant_mapper.add_feature(
                0,
                &segment.sequence,
                FeatureEntry::new(
                    constant_targets.len() as u64 + 1,
                    segment.name.clone(),
                    "vdj_constant_region",
                ),
            );
            debug_assert_eq!(feature_index, constant_targets.len());
            constant_targets.push(segment.id);
        }
    }
    if constant_targets.is_empty() {
        return Ok(report);
    }

    let mut reader = bam::Reader::from_path(path.as_ref()).with_context(|| {
        format!(
            "opening {} for VDJ evidence rescan",
            path.as_ref().display()
        )
    })?;
    if threads > 1 {
        reader
            .set_threads(threads)
            .context("configuring multithreaded BAM decoding for VDJ evidence rescan")?;
    }
    let mut evidence = RescanEvidenceVdj::default();
    let mut batch = Vec::<(u64, RescanRecord)>::with_capacity(RESCAN_EVIDENCE_BATCH_SIZE);
    let mut last_query: Option<(u64, Vec<u8>)> = None;
    let mut current_fragment_id = 0u64;
    let mut next_fragment_id = 0u64;

    // This is intentionally a true second pass: reopen the BAM and visit every
    // record to EOF. Expensive recombination/constant matching is restricted to
    // cells that already have reconstructed calls, but the BAM denominator and
    // live progress always reflect the complete file scan.
    for record in reader.records() {
        let record = record?;
        report.bam_records_scanned += 1;
        if report.bam_records_scanned % RESCAN_PROGRESS_EVERY_BAM_RECORDS == 0 {
            let (receptor_hit_records, constant_hit_records, linked_fragments, junction_support_reads, junction_spanning_reads, junction_conflicting_reads) =
                evidence.progress_counts();
            progress(RecombinationEvidenceRescanProgress {
                bam_records_scanned: report.bam_records_scanned,
                wanted_cell_records: report.wanted_cell_records,
                batches_completed,
                receptor_hit_records,
                constant_hit_records,
                linked_fragments,
                junction_support_reads,
                junction_spanning_reads,
                junction_conflicting_reads,
            });
        }
        let Some(cell) = resolver.cell(&record) else {
            continue;
        };
        let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
        if !wanted_cells.contains(&cell_id) {
            continue;
        }
        report.wanted_cell_records += 1;

        let is_new_query = last_query.as_ref().is_none_or(|(last_cell, last_qname)| {
            *last_cell != cell_id || last_qname.as_slice() != record.qname()
        });
        if is_new_query {
            if batch.len() >= RESCAN_EVIDENCE_BATCH_SIZE {
                let full =
                    std::mem::replace(&mut batch, Vec::with_capacity(RESCAN_EVIDENCE_BATCH_SIZE));
                evidence.consume_batch(
                    full,
                    &receptor_mapper,
                    &receptor_target_by_cell,
                    &constant_mapper,
                    &constant_targets,
                    index,
                    calls,
                );
                batches_completed = batches_completed.saturating_add(1);
                let (receptor_hit_records, constant_hit_records, linked_fragments, junction_support_reads, junction_spanning_reads, junction_conflicting_reads) =
                    evidence.progress_counts();
                progress(RecombinationEvidenceRescanProgress {
                    bam_records_scanned: report.bam_records_scanned,
                    wanted_cell_records: report.wanted_cell_records,
                    batches_completed,
                    receptor_hit_records,
                    constant_hit_records,
                    linked_fragments,
                    junction_support_reads,
                    junction_spanning_reads,
                    junction_conflicting_reads,
                });
            }
            current_fragment_id = next_fragment_id;
            next_fragment_id = next_fragment_id.wrapping_add(1);
            last_query = Some((cell_id, record.qname().to_vec()));
        }

        let sequence = record.seq().as_bytes();
        if sequence.len() < 8 {
            continue;
        }
        batch.push((
            cell_id,
            RescanRecord {
                fragment_id: current_fragment_id,
                sequence,
                qualities: record.qual().to_vec(),
            },
        ));
    }

    if !batch.is_empty() {
        evidence.consume_batch(
            batch,
            &receptor_mapper,
            &receptor_target_by_cell,
            &constant_mapper,
            &constant_targets,
            index,
            calls,
        );
        batches_completed = batches_completed.saturating_add(1);
        let (receptor_hit_records, constant_hit_records, linked_fragments, junction_support_reads, junction_spanning_reads, junction_conflicting_reads) =
            evidence.progress_counts();
        progress(RecombinationEvidenceRescanProgress {
            bam_records_scanned: report.bam_records_scanned,
            wanted_cell_records: report.wanted_cell_records,
            batches_completed,
            receptor_hit_records,
            constant_hit_records,
            linked_fragments,
            junction_support_reads,
            junction_spanning_reads,
            junction_conflicting_reads,
        });
    }

    let totals = evidence.into_totals();
    report.receptor_hit_records = totals.receptor_hit_records;
    report.constant_hit_records = totals.constant_hit_records;
    report.linked_fragments = totals.linked_fragments;
    report.junction_support_reads = totals.junction_support_reads;
    report.junction_spanning_reads = totals.junction_spanning_reads;
    report.junction_conflicting_reads = totals.junction_conflicting_reads;

    let junction_pileups = totals.junction_pileups;
    for (&(call_group, call_index), pileup) in &junction_pileups {
        let call = &mut calls[call_group].1[call_index];
        call.receptor_linkage.junction_support_reads = pileup.support_reads;
        call.receptor_linkage.junction_spanning_reads = pileup.spanning_reads;
        call.receptor_linkage.junction_conflicting_reads = pileup.conflicting_reads;

        let mut candidate = call.clone();
        let refined = refine_junction_from_pileup(&mut candidate, pileup);
        if refined > 0 && refresh_recombination_from_observed(&mut candidate, index) {
            candidate.receptor_linkage.junction_support_reads = pileup.support_reads;
            candidate.receptor_linkage.junction_spanning_reads = pileup.spanning_reads;
            candidate.receptor_linkage.junction_conflicting_reads = pileup.conflicting_reads;
            candidate.receptor_linkage.junction_refined_bases = refined as u16;
            *call = candidate;
            report.junction_refined_calls = report.junction_refined_calls.saturating_add(1);
            report.junction_refined_bases = report.junction_refined_bases.saturating_add(refined);
        }
    }

    for ((call_group, call_index), count) in totals.rediscovery_reads {
        calls[call_group].1[call_index]
            .receptor_linkage
            .rediscovery_reads = calls[call_group].1[call_index]
            .receptor_linkage
            .rediscovery_reads
            .saturating_add(count);
    }
    let support = totals.support;

    for (target_index, target) in receptor_targets.iter().enumerate() {
        let mut candidates: Vec<_> = support
            .iter()
            .filter_map(|(&(group, index_in_group, segment), &candidate)| {
                (group == target.call_group && index_in_group == target.call_index).then_some(
                    RecombinationRescanCandidate {
                        segment,
                        fragments: candidate.fragments,
                        same_read: candidate.spanning_reads,
                    },
                )
            })
            .collect();
        candidates.sort_by_key(|candidate| {
            (
                std::cmp::Reverse(candidate.same_read),
                std::cmp::Reverse(candidate.fragments),
                candidate.segment,
            )
        });
        report.calls[target_index].candidates = candidates;
    }

    let mut rescued = 0usize;
    for call_group in 0..calls.len() {
        for call_index in 0..calls[call_group].1.len() {
            let existing_constant = calls[call_group].1[call_index]
                .constant
                .as_ref()
                .map(|constant| constant.segment);
            let mut ranked: Vec<_> = support
                .iter()
                .filter_map(|(&(group, index_in_group, segment), &candidate)| {
                    (group == call_group && index_in_group == call_index)
                        .then_some((segment, candidate))
                })
                .collect();
            ranked.sort_by_key(|(_, candidate)| {
                (
                    std::cmp::Reverse(candidate.spanning_reads),
                    std::cmp::Reverse(candidate.fragments),
                )
            });

            let supported = if let Some(existing) = existing_constant {
                ranked
                    .iter()
                    .copied()
                    .find(|(segment, _)| *segment == existing)
            } else {
                let Some(best) = ranked.first().copied() else {
                    continue;
                };
                if ranked.get(1).is_some_and(|(_, second)| {
                    (second.spanning_reads, second.fragments)
                        == (best.1.spanning_reads, best.1.fragments)
                }) {
                    continue;
                }
                Some(best)
            };

            let Some((supported_id, supported_counts)) = supported else {
                continue;
            };
            calls[call_group].1[call_index]
                .receptor_linkage
                .constant_link_fragments = supported_counts.fragments;
            calls[call_group].1[call_index]
                .receptor_linkage
                .constant_spanning_reads = supported_counts.spanning_reads;
            calls[call_group].1[call_index]
                .receptor_linkage
                .constant_segment = Some(supported_id);

            if existing_constant.is_none() {
                let Some(best_segment) = index.segment(supported_id) else {
                    continue;
                };
                calls[call_group].1[call_index].constant = Some(ConstantRegionEvidence {
                    segment: supported_id,
                    supporting_features: supported_counts.fragments,
                    sequence: best_segment.sequence.clone(),
                });
                if let Some(target_index) = receptor_targets.iter().position(|target| {
                    target.call_group == call_group && target.call_index == call_index
                }) {
                    report.calls[target_index].rescued_segment = Some(supported_id);
                }
                rescued += 1;
            }
        }
    }

    for (target_index, target) in receptor_targets.iter().enumerate() {
        report.calls[target_index].stable_id =
            calls[target.call_group].1[target.call_index].stable_id.clone();
    }
    report.rescued = rescued;
    Ok(report)
}

fn receptor_bait(index: &VdjIndex, recombination: &Recombination) -> Option<Vec<u8>> {
    let observed = &recombination.observed_rearrangement;
    let j = index.segment(recombination.j)?;
    if observed.len() < 8 || j.sequence.len() < 8 {
        return None;
    }

    let j_alignment = local_alignment(observed, &j.sequence);
    let j_start = if j_alignment.score > 0 {
        j_alignment.query_start.min(observed.len())
    } else {
        observed.len().saturating_sub(J_BAIT_BASES)
    };
    let start = j_start.saturating_sub(CDR3_UPSTREAM_BASES);
    let end = (j_start + J_BAIT_BASES).min(observed.len());
    (end.saturating_sub(start) >= 8).then(|| observed[start..end].to_vec())
}

fn map_receptor_oriented(
    mapper: &FastLocusMapper,
    cell_id: u64,
    sequence: &[u8],
) -> Option<(usize, Vec<u8>, bool)> {
    if let MapStatus::Hit { feature_index, .. } = mapper.map_status(cell_id, sequence) {
        return Some((feature_index, sequence.to_vec(), false));
    }
    let reverse = reverse_complement(sequence);
    match mapper.map_status(cell_id, &reverse) {
        MapStatus::Hit { feature_index, .. } => Some((feature_index, reverse, true)),
        MapStatus::NoHit | MapStatus::Tie { .. } => None,
    }
}

#[derive(Debug)]
struct JunctionReadEvidence {
    reference_start: usize,
    counts: Vec<[u32; 4]>,
    quality_sums: Vec<[u32; 4]>,
    spans_junction: bool,
    conflicts: bool,
}

fn junction_read_evidence(
    sequence: &[u8],
    qualities: &[u8],
    call: &Recombination,
    index: &VdjIndex,
) -> Option<JunctionReadEvidence> {
    const ANCHOR_K: usize = 9;
    const ANCHOR_FLANK: usize = 18;

    let (junction_start, junction_end) = refinement_region(call, index)?;
    if junction_end.saturating_sub(junction_start) < MIN_JUNCTION_OVERLAP {
        return None;
    }
    let anchor_start = junction_start.saturating_sub(ANCHOR_FLANK);
    let anchor_end = (junction_end + ANCHOR_FLANK).min(call.observed_rearrangement.len());
    let offset = infer_ungapped_offset(
        sequence,
        &call.observed_rearrangement,
        anchor_start,
        anchor_end,
        ANCHOR_K,
    )?;

    let region_len = junction_end - junction_start;
    let mut counts = vec![[0u32; 4]; region_len];
    let mut quality_sums = vec![[0u32; 4]; region_len];
    let mut covered = 0usize;
    let mut conflicts = false;
    for reference_pos in junction_start..junction_end {
        let query_pos = reference_pos as isize + offset;
        if query_pos < 0 {
            continue;
        }
        let query_pos = query_pos as usize;
        let Some(&base) = sequence.get(query_pos) else { continue; };
        let Some(base_index) = base_index(base) else { continue; };
        let junction_pos = reference_pos - junction_start;
        counts[junction_pos][base_index] = counts[junction_pos][base_index].saturating_add(1);
        let quality = qualities.get(query_pos).copied().unwrap_or(0).min(60) as u32 + 1;
        quality_sums[junction_pos][base_index] = quality_sums[junction_pos][base_index]
            .saturating_add(quality);
        covered += 1;
        if call.observed_rearrangement[reference_pos].to_ascii_uppercase()
            != base.to_ascii_uppercase()
        {
            conflicts = true;
        }
    }
    if covered < MIN_JUNCTION_OVERLAP {
        return None;
    }

    Some(JunctionReadEvidence {
        reference_start: junction_start,
        counts,
        quality_sums,
        spans_junction: covered == region_len,
        conflicts,
    })
}

fn refinement_region(call: &Recombination, index: &VdjIndex) -> Option<(usize, usize)> {
    if !call.airr_junction.is_empty() {
        if let Some(start) = find_subslice(&call.observed_rearrangement, &call.airr_junction) {
            return Some((start, start + call.airr_junction.len()));
        }
    }

    // A Stage-2 call can lack an AIRR CDR3 when its initial junction is too
    // incomplete to establish the CDS-fixed conserved anchors. Still let the
    // rescan sharpen the observed V/J boundary; successful refinement is then
    // re-annotated and may acquire a valid CDR3/productivity call.
    let v = index.segment(call.v)?;
    let j = index.segment(call.j)?;
    let va = local_alignment(&call.observed_rearrangement, &v.sequence);
    let ja = local_alignment(&call.observed_rearrangement, &j.sequence);
    if va.score <= 0 || ja.score <= 0 || va.query_end > ja.query_start {
        return None;
    }
    let start = va.query_end.saturating_sub(12);
    let end = (ja.query_start + 12).min(call.observed_rearrangement.len());
    (end > start).then_some((start, end))
}

/// Infer a read-to-receptor ungapped coordinate offset from exact short anchors.
/// This is intentionally much cheaper than dynamic-programming alignment for
/// every Stage-3 receptor hit. A unique offset supported by at least two exact
/// anchors is required before a read is allowed to vote on junction bases.
fn infer_ungapped_offset(
    query: &[u8],
    reference: &[u8],
    start: usize,
    end: usize,
    k: usize,
) -> Option<isize> {
    if k == 0 || query.len() < k || end.saturating_sub(start) < k || end > reference.len() {
        return None;
    }
    let mut votes = HashMap::<isize, u16>::new();
    let last = end - k;
    let mut reference_pos = start;
    while reference_pos <= last {
        let anchor = &reference[reference_pos..reference_pos + k];
        for (query_pos, window) in query.windows(k).enumerate() {
            if window == anchor {
                let offset = query_pos as isize - reference_pos as isize;
                let vote = votes.entry(offset).or_default();
                *vote = vote.saturating_add(1);
            }
        }
        reference_pos = reference_pos.saturating_add(4);
    }
    let mut ranked: Vec<_> = votes.into_iter().collect();
    ranked.sort_by_key(|(offset, count)| (std::cmp::Reverse(*count), *offset));
    let (best_offset, best_votes) = ranked.first().copied()?;
    if best_votes < 2 {
        return None;
    }
    if ranked.get(1).is_some_and(|(_, second_votes)| *second_votes == best_votes) {
        return None;
    }
    Some(best_offset)
}

fn merge_junction_read(dst: &mut JunctionPileup, src: JunctionReadEvidence) {
    if dst.reference_start.is_some_and(|start| start != src.reference_start) {
        return;
    }
    dst.reference_start = Some(src.reference_start);
    ensure_pileup_len(dst, src.counts.len());
    dst.support_reads = dst.support_reads.saturating_add(1);
    if src.spans_junction {
        dst.spanning_reads = dst.spanning_reads.saturating_add(1);
    }
    if src.conflicts {
        dst.conflicting_reads = dst.conflicting_reads.saturating_add(1);
    }
    for pos in 0..src.counts.len() {
        for base in 0..4 {
            dst.counts[pos][base] = dst.counts[pos][base].saturating_add(src.counts[pos][base]);
            dst.quality_sums[pos][base] = dst.quality_sums[pos][base]
                .saturating_add(src.quality_sums[pos][base]);
        }
    }
}

fn merge_junction_pileup(dst: &mut JunctionPileup, src: JunctionPileup) {
    if let Some(start) = src.reference_start {
        if dst.reference_start.is_some_and(|existing| existing != start) {
            return;
        }
        dst.reference_start = Some(start);
    }
    ensure_pileup_len(dst, src.counts.len());
    dst.support_reads = dst.support_reads.saturating_add(src.support_reads);
    dst.spanning_reads = dst.spanning_reads.saturating_add(src.spanning_reads);
    dst.conflicting_reads = dst.conflicting_reads.saturating_add(src.conflicting_reads);
    for pos in 0..src.counts.len() {
        for base in 0..4 {
            dst.counts[pos][base] = dst.counts[pos][base].saturating_add(src.counts[pos][base]);
            dst.quality_sums[pos][base] = dst.quality_sums[pos][base]
                .saturating_add(src.quality_sums[pos][base]);
        }
    }
}

fn ensure_pileup_len(pileup: &mut JunctionPileup, len: usize) {
    if pileup.counts.len() < len {
        pileup.counts.resize(len, [0; 4]);
        pileup.quality_sums.resize(len, [0; 4]);
    }
}

fn refine_junction_from_pileup(call: &mut Recombination, pileup: &JunctionPileup) -> usize {
    let Some(junction_start) = pileup.reference_start else { return 0; };
    let mut refined = 0usize;
    for (pos, counts) in pileup.counts.iter().enumerate() {
        let depth: u32 = counts.iter().sum();
        if depth < MIN_REFINE_DEPTH {
            continue;
        }
        let Some((best_index, &best_count)) = counts
            .iter()
            .enumerate()
            .max_by_key(|(base, count)| (**count, pileup.quality_sums[pos][*base]))
        else {
            continue;
        };
        if best_count < MIN_REFINE_DEPTH || best_count.saturating_mul(3) < depth.saturating_mul(2) {
            continue;
        }
        let best = index_base(best_index);
        let absolute = junction_start + pos;
        if call.observed_rearrangement.get(absolute).copied() == Some(best) {
            continue;
        }
        call.observed_rearrangement[absolute] = best;
        if absolute < call.observed_receptor_sequence.len() {
            call.observed_receptor_sequence[absolute] = best;
        }
        refined += 1;
    }
    refined
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack.windows(needle.len()).position(|window| window == needle)
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn index_base(index: usize) -> u8 {
    [b'A', b'C', b'G', b'T'][index]
}

fn map_constant(mapper: &FastLocusMapper, sequence: &[u8]) -> Option<usize> {
    if let MapStatus::Hit { feature_index, .. } = mapper.map_status(0, sequence) {
        return Some(feature_index);
    }
    let reverse = reverse_complement(sequence);
    match mapper.map_status(0, &reverse) {
        MapStatus::Hit { feature_index, .. } => Some(feature_index),
        MapStatus::NoHit | MapStatus::Tie { .. } => None,
    }
}
