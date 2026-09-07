use super::{
    local_alignment, ConstantRegionEvidence, Recombination, RecombinationId,
};
use crate::index::{reverse_complement, SegmentId, SegmentKind, VdjIndex};
use crate::runner::BamIdentityResolver;
use anyhow::{Context, Result};
use fast_tag_mapper::{FastLocusMapper, FeatureEntry, MapStatus};
use int_to_str::IntToStr;
use rust_htslib::bam::{self, Read};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::Path;

const CDR3_UPSTREAM_BASES: usize = 32;
const J_BAIT_BASES: usize = 24;

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
    pub rescued: usize,
    pub calls: Vec<RecombinationRescanCall>,
}

#[derive(Debug, Default)]
struct FragmentHits {
    receptors: HashSet<(usize, usize)>,
    constants: HashSet<SegmentId>,
    spanning_reads: HashMap<(usize, usize, SegmentId), u32>,
}

pub(crate) fn rescue_missing_constants_from_bam<P: AsRef<Path>, R: BamIdentityResolver>(
    path: P,
    resolver: &R,
    index: &VdjIndex,
    calls: &mut [(u64, Vec<Recombination>)],
    threads: usize,
) -> Result<usize> {
    Ok(rescue_missing_constants_from_bam_with_report(
        path, resolver, index, calls, threads,
    )?
    .rescued)
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
    let mut report = RecombinationEvidenceRescanReport::default();
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

    let mut reader = bam::Reader::from_path(path.as_ref())
        .with_context(|| format!("opening {} for VDJ evidence rescan", path.as_ref().display()))?;
    if threads > 1 {
        reader
            .set_threads(threads)
            .context("configuring multithreaded BAM decoding for VDJ evidence rescan")?;
    }
    let mut fragments = HashMap::<(u64, Vec<u8>), FragmentHits>::new();
    for record in reader.records() {
        let record = record?;
        report.bam_records_scanned += 1;
        let Some(cell) = resolver.cell(&record) else {
            continue;
        };
        let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
        if !wanted_cells.contains(&cell_id) {
            continue;
        }
        report.wanted_cell_records += 1;
        let sequence = record.seq().as_bytes();
        if sequence.len() < 8 {
            continue;
        }

        let receptor = map_receptor(&receptor_mapper, cell_id, &sequence)
            .and_then(|feature_index| receptor_target_by_cell.get(&(feature_index, cell_id)).copied())
            .map(|target| (target.call_group, target.call_index));
        let constant = map_constant(&constant_mapper, &sequence)
            .and_then(|feature_index| constant_targets.get(feature_index).copied());
        if let Some((call_group, call_index)) = receptor {
            report.receptor_hit_records += 1;
            calls[call_group].1[call_index].receptor_linkage.rediscovery_reads = calls[call_group].1
                [call_index]
                .receptor_linkage
                .rediscovery_reads
                .saturating_add(1);
        }
        if constant.is_some() {
            report.constant_hit_records += 1;
        }

        if receptor.is_none() && constant.is_none() {
            continue;
        }

        let fragment = fragments
            .entry((cell_id, record.qname().to_vec()))
            .or_default();
        if let Some(target) = receptor {
            fragment.receptors.insert(target);
        }
        if let Some(segment) = constant {
            fragment.constants.insert(segment);
        }
        if let (Some((call_group, call_index)), Some(segment)) = (receptor, constant) {
            *fragment
                .spanning_reads
                .entry((call_group, call_index, segment))
                .or_insert(0) += 1;
        }
    }

    report.linked_fragments = fragments
        .values()
        .filter(|fragment| !fragment.receptors.is_empty() && !fragment.constants.is_empty())
        .count();

    let mut support = HashMap::<(usize, usize, SegmentId), CandidateSupport>::new();
    for fragment in fragments.values() {
        for &(call_group, call_index) in &fragment.receptors {
            let chain = calls[call_group].1[call_index].chain;
            for &constant in &fragment.constants {
                if index.segment(constant).is_none_or(|segment| segment.chain != chain) {
                    continue;
                }
                let candidate = support
                    .entry((call_group, call_index, constant))
                    .or_default();
                candidate.fragments += 1;
                candidate.spanning_reads += fragment
                    .spanning_reads
                    .get(&(call_group, call_index, constant))
                    .copied()
                    .unwrap_or(0);
            }
        }
    }

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
                ranked.iter().copied().find(|(segment, _)| *segment == existing)
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
            calls[call_group].1[call_index].receptor_linkage.constant_link_fragments =
                supported_counts.fragments;
            calls[call_group].1[call_index].receptor_linkage.constant_spanning_reads =
                supported_counts.spanning_reads;
            calls[call_group].1[call_index].receptor_linkage.constant_segment = Some(supported_id);

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

fn map_receptor(mapper: &FastLocusMapper, cell_id: u64, sequence: &[u8]) -> Option<usize> {
    if let MapStatus::Hit { feature_index, .. } = mapper.map_status(cell_id, sequence) {
        return Some(feature_index);
    }
    let reverse = reverse_complement(sequence);
    match mapper.map_status(cell_id, &reverse) {
        MapStatus::Hit { feature_index, .. } => Some(feature_index),
        MapStatus::NoHit | MapStatus::Tie { .. } => None,
    }
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
