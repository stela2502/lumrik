use crate::cellrep::{
    AlignmentGeometry, BamFeatureEvidence, BamFeatureSequenceParts, CellEvidence, CellEvidenceVdj,
    EvidenceId, MapperEvidence, SequencePart,
};
use crate::index::{SegmentKind, VdjIndex};
use crate::recombination::{
    process_chain_work, rescue_missing_constants_from_bam,
    rescue_missing_constants_from_bam_with_report,
    rescue_missing_constants_from_bam_with_report_and_progress, ChainWork, Recombination,
    RecombinationEvidenceRescanProgress, RecombinationEvidenceRescanReport,
};
use anyhow::{Context, Result};
use int_to_str::IntToStr;
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use std::collections::{HashMap, HashSet};
use std::path::Path;

const EVIDENCE_BATCH_SIZE: usize = 200_000;
const UNMAPPED_IGH_SEED_LEN: usize = 13;
const UNMAPPED_IGH_MIN_J_SEEDS: usize = 2;
const UNMAPPED_IGH_MIN_J_SEED_SPAN: usize = UNMAPPED_IGH_SEED_LEN;
const UNMAPPED_IGH_MIN_HARVEST_SEEDS: usize = 2;
const UNMAPPED_CANDIDATE_BATCH_SIZE: usize = 100_000;

#[derive(Debug, Clone)]
struct UnmappedReadCandidate {
    cell_id: u64,
    id: EvidenceId,
    bases: Vec<u8>,
    qualities: Vec<u8>,
    is_last_in_template: bool,
    is_secondary: bool,
    is_supplementary: bool,
}

#[derive(Default)]
struct SeedSupport {
    distinct: HashSet<Vec<u8>>,
    min_read_pos: usize,
    max_read_pos: usize,
    initialized: bool,
}

impl SeedSupport {
    fn observe(&mut self, seed: &[u8], read_pos: usize) {
        self.distinct.insert(seed.to_vec());
        if self.initialized {
            self.min_read_pos = self.min_read_pos.min(read_pos);
            self.max_read_pos = self.max_read_pos.max(read_pos);
        } else {
            self.min_read_pos = read_pos;
            self.max_read_pos = read_pos;
            self.initialized = true;
        }
    }

    fn distinct_hits(&self) -> usize {
        self.distinct.len()
    }

    fn start_span(&self) -> usize {
        self.max_read_pos.saturating_sub(self.min_read_pos)
    }
}

/// Exact IGH germline seeds used only for the deferred unmapped-read rescue.
///
/// BAM ingestion deliberately retains every barcoded unmapped record without
/// deciding whether it is receptor evidence. After the mapped evidence pass is
/// complete, candidates are searched in parallel. Admission is strict and
/// IGHJ-driven: at least two distinct J-specific seeds must support the same J
/// segment across a real span of the read. Once admitted, V/D/J/C support is
/// harvested more broadly and handed to the ordinary CellEvidenceVdj compactor.
struct UnmappedIghSeeds {
    by_seed: HashMap<Vec<u8>, Vec<crate::index::SegmentId>>,
    j_specific_by_seed: HashMap<Vec<u8>, Vec<crate::index::SegmentId>>,
}

impl UnmappedIghSeeds {
    fn new(index: &VdjIndex) -> Self {
        let mut all_vdj_by_seed = HashMap::<Vec<u8>, Vec<crate::index::SegmentId>>::new();
        for segment in &index.segments {
            if segment.sequence.len() < UNMAPPED_IGH_SEED_LEN {
                continue;
            }
            for seed in segment.sequence.windows(UNMAPPED_IGH_SEED_LEN) {
                let ids = all_vdj_by_seed.entry(seed.to_vec()).or_default();
                if !ids.contains(&segment.id) {
                    ids.push(segment.id);
                }
            }
        }

        let mut by_seed = HashMap::<Vec<u8>, Vec<crate::index::SegmentId>>::new();
        for kind in [SegmentKind::V, SegmentKind::D, SegmentKind::J, SegmentKind::C] {
            for segment in index.segments_for(crate::index::Chain::Igh, kind) {
                if segment.sequence.len() < UNMAPPED_IGH_SEED_LEN {
                    continue;
                }
                for seed in segment.sequence.windows(UNMAPPED_IGH_SEED_LEN) {
                    let ids = by_seed.entry(seed.to_vec()).or_default();
                    if !ids.contains(&segment.id) {
                        ids.push(segment.id);
                    }
                }
            }
        }

        let j_specific_by_seed = all_vdj_by_seed
            .into_iter()
            .filter_map(|(seed, ids)| {
                let all_igh_j = ids.iter().all(|id| {
                    index.segment(*id).is_some_and(|segment| {
                        segment.chain == crate::index::Chain::Igh
                            && segment.kind == SegmentKind::J
                    })
                });
                all_igh_j.then_some((seed, ids))
            })
            .collect();

        Self {
            by_seed,
            j_specific_by_seed,
        }
    }

    fn rescue_candidate(
        &self,
        candidate: UnmappedReadCandidate,
        index: &VdjIndex,
    ) -> Option<(u64, BamFeatureEvidence)> {
        let reverse = crate::index::reverse_complement(&candidate.bases);

        let forward_j = self.j_support(&candidate.bases);
        let reverse_j = self.j_support(&reverse);
        let (read_is_reverse_to_transcript, mut support, j_support) =
            if self.has_strong_j(&forward_j) {
                (false, self.all_support(&candidate.bases), forward_j)
            } else if self.has_strong_j(&reverse_j) {
                (true, self.all_support(&reverse), reverse_j)
            } else {
                return None;
            };
        // Preserve every J segment that passed the strict admission gate even
        // if family-shared seeds make its broad-support count less impressive.
        for (id, j) in j_support {
            support.entry(id).or_insert(j);
        }

        let mut segment_ids = Vec::new();
        for kind in [SegmentKind::V, SegmentKind::D, SegmentKind::J, SegmentKind::C] {
            let best = support
                .iter()
                .filter(|(id, evidence)| {
                    evidence.distinct_hits() >= UNMAPPED_IGH_MIN_HARVEST_SEEDS
                        && index.segment(**id).is_some_and(|s| s.kind == kind)
                })
                .map(|(_, evidence)| evidence.distinct_hits())
                .max()
                .unwrap_or(0);
            if best == 0 {
                continue;
            }
            segment_ids.extend(support.iter().filter_map(|(id, evidence)| {
                (evidence.distinct_hits() == best
                    && index.segment(*id).is_some_and(|s| s.kind == kind))
                .then_some(*id)
            }));
        }
        segment_ids.sort_unstable();
        segment_ids.dedup();
        if !segment_ids
            .iter()
            .any(|id| index.segment(*id).is_some_and(|s| s.kind == SegmentKind::J))
        {
            return None;
        }

        let mappings = segment_ids
            .into_iter()
            .map(|segment_id| {
                let segment_is_reverse = index.segment(segment_id).is_some_and(|segment| {
                    matches!(segment.strand, crate::index::Strand::Minus)
                });
                MapperEvidence {
                    segment_id,
                    alignment: AlignmentGeometry {
                        tid: -1,
                        start: 0,
                        end: 0,
                        // summarize_chain_work converts mapper orientation back
                        // to transcript orientation by XORing the segment strand.
                        is_reverse: read_is_reverse_to_transcript ^ segment_is_reverse,
                        is_secondary: candidate.is_secondary,
                        is_supplementary: candidate.is_supplementary,
                        mapq: 0,
                        ref_blocks: Vec::new(),
                    },
                }
            })
            .collect();
        let cell_id = candidate.cell_id;
        let part = SequencePart {
            bases: candidate.bases,
            qualities: candidate.qualities,
        };
        let mut sequence_parts = BamFeatureSequenceParts::default();
        if candidate.is_last_in_template {
            sequence_parts.r2 = Some(part);
        } else {
            sequence_parts.r1 = Some(part);
        }

        Some((
            cell_id,
            BamFeatureEvidence {
                id: candidate.id,
                sequence: sequence_parts,
                mappings,
                intronic_constant_segments: Vec::new(),
            },
        ))
    }

    fn has_strong_j(&self, support: &HashMap<crate::index::SegmentId, SeedSupport>) -> bool {
        support.values().any(|evidence| {
            evidence.distinct_hits() >= UNMAPPED_IGH_MIN_J_SEEDS
                && evidence.start_span() >= UNMAPPED_IGH_MIN_J_SEED_SPAN
        })
    }

    fn j_support(&self, sequence: &[u8]) -> HashMap<crate::index::SegmentId, SeedSupport> {
        self.support_from_table(sequence, &self.j_specific_by_seed)
    }

    fn all_support(&self, sequence: &[u8]) -> HashMap<crate::index::SegmentId, SeedSupport> {
        self.support_from_table(sequence, &self.by_seed)
    }

    fn support_from_table(
        &self,
        sequence: &[u8],
        table: &HashMap<Vec<u8>, Vec<crate::index::SegmentId>>,
    ) -> HashMap<crate::index::SegmentId, SeedSupport> {
        let mut support = HashMap::<crate::index::SegmentId, SeedSupport>::new();
        if sequence.len() < UNMAPPED_IGH_SEED_LEN {
            return support;
        }
        for (read_pos, seed) in sequence.windows(UNMAPPED_IGH_SEED_LEN).enumerate() {
            if let Some(ids) = table.get(seed) {
                for &id in ids {
                    support.entry(id).or_default().observe(seed, read_pos);
                }
            }
        }
        support
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChainKneeSelection {
    pub chain: crate::index::Chain,
    pub threshold_records: usize,
    pub evidence_cells: usize,
    pub selected_cells: usize,
}

pub trait BamIdentityResolver {
    fn cell(&self, record: &bam::Record) -> Option<String>;
}
#[derive(Debug, Clone, Copy, Default)]
pub struct NelruneIdentityResolver;
impl BamIdentityResolver for NelruneIdentityResolver {
    fn cell(&self, record: &bam::Record) -> Option<String> {
        if let Ok(Aux::String(s)) = record.aux(b"CB") {
            return Some(s.to_string());
        }
        let fields: Vec<_> = record.qname().split(|b| *b == b'|').collect();
        decode_hex_ascii(fields.get(1).copied()?)
    }
}

#[derive(Debug, Clone)]
pub struct VdjRunnerConfig {
    pub min_sequence_overlap: usize,
}
impl Default for VdjRunnerConfig {
    fn default() -> Self {
        Self {
            min_sequence_overlap: 12,
        }
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct BamIngestProgress {
    pub bam_records: usize,
    pub allowed_cell_records: usize,
    pub receptor_overlap_records: usize,
    pub unmapped_candidates: usize,
    pub unmapped_igh_admitted: usize,
    pub unmapped_igh_rescued_cells: usize,
    pub unmapped_igh_v_mappings: usize,
    pub unmapped_igh_d_mappings: usize,
    pub unmapped_igh_j_mappings: usize,
    pub unmapped_igh_c_mappings: usize,
}

#[derive(Debug, Default)]
struct UnmappedIghRescueStats {
    candidates: usize,
    admitted: usize,
    rescued_cells: HashSet<u64>,
    v_mappings: usize,
    d_mappings: usize,
    j_mappings: usize,
    c_mappings: usize,
}

impl UnmappedIghRescueStats {
    fn add_assign(&mut self, other: Self) {
        self.candidates = self.candidates.saturating_add(other.candidates);
        self.admitted = self.admitted.saturating_add(other.admitted);
        self.rescued_cells.extend(other.rescued_cells);
        self.v_mappings = self.v_mappings.saturating_add(other.v_mappings);
        self.d_mappings = self.d_mappings.saturating_add(other.d_mappings);
        self.j_mappings = self.j_mappings.saturating_add(other.j_mappings);
        self.c_mappings = self.c_mappings.saturating_add(other.c_mappings);
    }

    fn progress(&self, bam_records: usize, allowed_cell_records: usize, receptor_overlap_records: usize) -> BamIngestProgress {
        BamIngestProgress {
            bam_records,
            allowed_cell_records,
            receptor_overlap_records,
            unmapped_candidates: self.candidates,
            unmapped_igh_admitted: self.admitted,
            unmapped_igh_rescued_cells: self.rescued_cells.len(),
            unmapped_igh_v_mappings: self.v_mappings,
            unmapped_igh_d_mappings: self.d_mappings,
            unmapped_igh_j_mappings: self.j_mappings,
            unmapped_igh_c_mappings: self.c_mappings,
        }
    }
}

pub struct VdjRunner {
    pub index: VdjIndex,
    pub evidence: CellEvidenceVdj,
    pub cell_names: HashMap<u64, String>,
    pub config: VdjRunnerConfig,
    flush_id: u32,
    threads: usize,
}
impl VdjRunner {
    pub fn new(index: VdjIndex, config: VdjRunnerConfig) -> Self {
        Self {
            index,
            evidence: CellEvidenceVdj::new(),
            cell_names: HashMap::new(),
            config,
            flush_id: 0,
            threads: 1,
        }
    }
    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads.max(1);
    }

    fn consume_unmapped_candidate_batch(
        &mut self,
        candidates: &mut Vec<UnmappedReadCandidate>,
        seeds: &UnmappedIghSeeds,
        pool: Option<&rayon::ThreadPool>,
    ) -> UnmappedIghRescueStats {
        if candidates.is_empty() {
            return UnmappedIghRescueStats::default();
        }

        let batch = std::mem::replace(
            candidates,
            Vec::with_capacity(UNMAPPED_CANDIDATE_BATCH_SIZE),
        );
        let candidate_count = batch.len();
        let rescue_one = |candidate| seeds.rescue_candidate(candidate, &self.index);
        let rescued: Vec<(u64, BamFeatureEvidence)> = if let Some(pool) = pool {
            pool.install(|| batch.into_par_iter().filter_map(rescue_one).collect())
        } else {
            batch.into_iter().filter_map(rescue_one).collect()
        };

        let mut stats = UnmappedIghRescueStats {
            candidates: candidate_count,
            admitted: rescued.len(),
            ..Default::default()
        };
        let mut rescued_cells = HashSet::new();
        for (cell_id, feature) in &rescued {
            rescued_cells.insert(*cell_id);
            for mapping in &feature.mappings {
                match self.index.segment(mapping.segment_id).map(|s| s.kind) {
                    Some(SegmentKind::V) => stats.v_mappings = stats.v_mappings.saturating_add(1),
                    Some(SegmentKind::D) => stats.d_mappings = stats.d_mappings.saturating_add(1),
                    Some(SegmentKind::J) => stats.j_mappings = stats.j_mappings.saturating_add(1),
                    Some(SegmentKind::C) => stats.c_mappings = stats.c_mappings.saturating_add(1),
                    None => {}
                }
            }
        }
        stats.rescued_cells = rescued_cells;

        if !rescued.is_empty() {
            self.evidence.consume_batch(
                rescued,
                &self.index,
                self.config.min_sequence_overlap,
            );
        }
        stats
    }

    pub fn read_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &mut self,
        path: P,
        resolver: &R,
    ) -> Result<usize> {
        self.read_bam_with_progress(path, resolver, |_, _, _| {})
    }

    pub fn read_bam_with_progress<P, R, F>(
        &mut self,
        path: P,
        resolver: &R,
        progress: F,
    ) -> Result<usize>
    where
        P: AsRef<Path>,
        R: BamIdentityResolver,
        F: FnMut(BamIngestProgress, &CellEvidenceVdj, &VdjIndex),
    {
        self.read_bam_with_progress_for_cells(path, resolver, None, progress)
    }

    /// Read receptor evidence while optionally restricting ingestion to a
    /// preliminary set of valid cell IDs. The allowed set is only a gate: it
    /// does not pre-populate CellEvidenceVdj. Cells enter CellHash only after
    /// an allowed BAM record contributes actual receptor evidence.
    pub fn read_bam_with_progress_for_cells<P, R, F>(
        &mut self,
        path: P,
        resolver: &R,
        allowed_cells: Option<&HashSet<u64>>,
        mut progress: F,
    ) -> Result<usize>
    where
        P: AsRef<Path>,
        R: BamIdentityResolver,
        F: FnMut(BamIngestProgress, &CellEvidenceVdj, &VdjIndex),
    {
        let mut reader = bam::Reader::from_path(path.as_ref())
            .with_context(|| format!("opening {}", path.as_ref().display()))?;
        if self.threads > 1 {
            reader
                .set_threads(self.threads)
                .context("configuring multithreaded BAM decoding")?;
        }
        let header = reader.header().to_owned();
        let mut n = 0usize;
        let mut bam_records = 0usize;
        let mut allowed_cell_records = 0usize;
        let mut unmapped_rescue = UnmappedIghRescueStats::default();
        let mut entry = 0u32;
        let mut batch = Vec::<(u64, BamFeatureEvidence)>::with_capacity(EVIDENCE_BATCH_SIZE);
        let mut unmapped_candidates =
            Vec::<UnmappedReadCandidate>::with_capacity(UNMAPPED_CANDIDATE_BATCH_SIZE);
        let unmapped_igh_seeds = UnmappedIghSeeds::new(&self.index);
        let unmapped_pool = (self.threads > 1)
            .then(|| {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(self.threads)
                    .build()
                    .expect("building sc-vdj unmapped rescue Rayon pool")
            });

        // Mapper BAMs emit the records belonging to one physical query together.
        // Keep only the immediately preceding query key so paired/supplementary
        // records share an EvidenceId without retaining a BAM-sized QNAME map.
        // A full batch is flushed only when a new query begins, so one physical
        // fragment is never split merely because it crossed the evidence-batch boundary.
        let mut last_query: Option<(u64, Vec<u8>)> = None;
        let mut current_evidence_id: Option<EvidenceId> = None;

        for rec in reader.records() {
            let rec = rec?;
            bam_records = bam_records.saturating_add(1);
            let Some(cell) = resolver.cell(&rec) else {
                continue;
            };
            let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
            let query_key = (cell_id, rec.qname().to_vec());

            if last_query.as_ref() != Some(&query_key) {
                if batch.len() >= EVIDENCE_BATCH_SIZE {
                    let full =
                        std::mem::replace(&mut batch, Vec::with_capacity(EVIDENCE_BATCH_SIZE));
                    self.evidence.consume_batch(
                        full,
                        &self.index,
                        self.config.min_sequence_overlap,
                    );
                    progress(
                        unmapped_rescue.progress(bam_records, allowed_cell_records, n),
                        &self.evidence,
                        &self.index,
                    );
                }
                if unmapped_candidates.len() >= UNMAPPED_CANDIDATE_BATCH_SIZE {
                    let rescue = self.consume_unmapped_candidate_batch(
                        &mut unmapped_candidates,
                        &unmapped_igh_seeds,
                        unmapped_pool.as_ref(),
                    );
                    n = n.saturating_add(rescue.admitted);
                    unmapped_rescue.add_assign(rescue);
                    progress(
                        unmapped_rescue.progress(bam_records, allowed_cell_records, n),
                        &self.evidence,
                        &self.index,
                    );
                }
                last_query = Some(query_key);
                current_evidence_id = None;
            }

            // Do not ask STAR to decide whether an unmapped molecule deserves
            // to participate in receptor discovery. Retain every barcoded
            // unmapped record now; classify it only after the mapped evidence
            // pass, in parallel, against the receptor index. This also bypasses
            // the preliminary exonic-cell gate for genuinely rescued cells.
            if rec.is_unmapped() {
                let id = *current_evidence_id.get_or_insert_with(|| {
                    let id = EvidenceId {
                        flush: self.flush_id,
                        entry,
                    };
                    entry = entry.wrapping_add(1);
                    id
                });
                unmapped_candidates.push(UnmappedReadCandidate {
                    cell_id,
                    id,
                    bases: rec.seq().as_bytes(),
                    qualities: rec.qual().to_vec(),
                    is_last_in_template: rec.is_last_in_template(),
                    is_secondary: rec.is_secondary(),
                    is_supplementary: rec.is_supplementary(),
                });
                self.cell_names.entry(cell_id).or_insert(cell);
                continue;
            }

            if allowed_cells.is_some_and(|allowed| !allowed.contains(&cell_id)) {
                continue;
            }
            allowed_cell_records = allowed_cell_records.saturating_add(1);

            let tid = rec.tid();
            if tid < 0 {
                continue;
            }
            let chr = String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned();
            let blocks = aligned_blocks(&rec);
            let segment_ids = self.index.overlapping(&chr, &blocks);
            let intronic_constant_segments =
                self.index.intronic_constant_overlapping(&chr, &blocks, 8);
            let has_rearrangement_segment = segment_ids.iter().any(|id| {
                self.index
                    .segment(*id)
                    .is_some_and(|segment| segment.kind != SegmentKind::C)
            });
            // C-exon-only reads are ambiguous: they may belong to a normal
            // rearranged transcript whose fragment simply does not reach J.
            // They neither classify biology nor seed expensive reconstruction.
            if !has_rearrangement_segment && intronic_constant_segments.is_empty() {
                continue;
            }

            let id = *current_evidence_id.get_or_insert_with(|| {
                let id = EvidenceId {
                    flush: self.flush_id,
                    entry,
                };
                entry = entry.wrapping_add(1);
                id
            });

            let geometry = AlignmentGeometry {
                tid,
                start: blocks.first().map_or(0, |x| x.0),
                end: blocks.last().map_or(0, |x| x.1),
                is_reverse: rec.is_reverse(),
                is_secondary: rec.is_secondary(),
                is_supplementary: rec.is_supplementary(),
                mapq: rec.mapq(),
                ref_blocks: blocks,
            };
            let mappings = segment_ids
                .into_iter()
                .map(|segment_id| MapperEvidence {
                    segment_id,
                    alignment: geometry.clone(),
                })
                .collect();
            let mut sequence = BamFeatureSequenceParts::default();
            if has_rearrangement_segment {
                let part = SequencePart {
                    bases: rec.seq().as_bytes(),
                    qualities: rec.qual().to_vec(),
                };
                if rec.is_last_in_template() {
                    sequence.r2 = Some(part)
                } else {
                    sequence.r1 = Some(part)
                };
            }

            batch.push((
                cell_id,
                BamFeatureEvidence {
                    id,
                    sequence,
                    mappings,
                    intronic_constant_segments,
                },
            ));
            self.cell_names.entry(cell_id).or_insert(cell);
            n += 1;
        }

        if !batch.is_empty() {
            self.evidence
                .consume_batch(batch, &self.index, self.config.min_sequence_overlap);
            progress(
                unmapped_rescue.progress(bam_records, allowed_cell_records, n),
                &self.evidence,
                &self.index,
            );
        }

        if !unmapped_candidates.is_empty() {
            let rescue = self.consume_unmapped_candidate_batch(
                &mut unmapped_candidates,
                &unmapped_igh_seeds,
                unmapped_pool.as_ref(),
            );
            n = n.saturating_add(rescue.admitted);
            unmapped_rescue.add_assign(rescue);
            progress(
                unmapped_rescue.progress(bam_records, allowed_cell_records, n),
                &self.evidence,
                &self.index,
            );
        }

        self.flush_id = self.flush_id.wrapping_add(1);
        Ok(n)
    }
    /// Re-scan the retained BAM with cell-specific reconstructed receptor baits.
    ///
    /// This independently rediscovers receptor reads and records direct
    /// receptor-to-constant linkage support. Missing constant calls may be
    /// filled when one constant segment has uniquely stronger direct support.
    pub fn rediscover_receptor_linkage_from_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<usize> {
        rescue_missing_constants_from_bam(path, resolver, &self.index, calls, self.threads)
    }

    pub fn rediscover_receptor_linkage_from_bam_with_report<
        P: AsRef<Path>,
        R: BamIdentityResolver,
    >(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<RecombinationEvidenceRescanReport> {
        rescue_missing_constants_from_bam_with_report(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
        )
    }

    pub fn rediscover_receptor_linkage_from_bam_with_report_and_progress<P, R, F>(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
        progress: F,
    ) -> Result<RecombinationEvidenceRescanReport>
    where
        P: AsRef<Path>,
        R: BamIdentityResolver,
        F: FnMut(RecombinationEvidenceRescanProgress),
    {
        rescue_missing_constants_from_bam_with_report_and_progress(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
            progress,
        )
    }

    pub fn rescue_missing_constants_from_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<usize> {
        rescue_missing_constants_from_bam(path, resolver, &self.index, calls, self.threads)
    }

    pub fn rescue_missing_constants_from_bam_with_report<P: AsRef<Path>, R: BamIdentityResolver>(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<RecombinationEvidenceRescanReport> {
        rescue_missing_constants_from_bam_with_report(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
        )
    }
    pub fn receptor_knee_selection(&self) -> Vec<ChainKneeSelection> {
        crate::index::Chain::ALL
            .into_iter()
            .map(|chain| {
                let mut counts: Vec<usize> = self
                    .evidence
                    .cells()
                    .filter_map(|(_, cell)| {
                        let v_mappings = cell.segment_mappings(chain, SegmentKind::V);
                        let j_mappings = cell.segment_mappings(chain, SegmentKind::J);
                        qualifies_reconstruction_candidate(v_mappings, j_mappings)
                            .then(|| cell.reconstruction_records(chain))
                    })
                    .collect();
                counts.sort_unstable_by(|a, b| b.cmp(a));

                let threshold_records = receptor_knee_threshold(&counts);
                let selected_cells = counts
                    .iter()
                    .filter(|count| **count >= threshold_records)
                    .count();

                ChainKneeSelection {
                    chain,
                    threshold_records,
                    evidence_cells: counts.len(),
                    selected_cells,
                }
            })
            .collect()
    }

    pub fn identify_with_receptor_knees(
        &self,
    ) -> (Vec<(u64, Vec<Recombination>)>, Vec<ChainKneeSelection>) {
        let selection = self.receptor_knee_selection();
        let out = self.identify_with_receptor_knee_selection(&selection);
        (out, selection)
    }

    pub fn identify_with_receptor_knee_selection(
        &self,
        selection: &[ChainKneeSelection],
    ) -> Vec<(u64, Vec<Recombination>)> {
        let thresholds: HashMap<_, _> = selection
            .iter()
            .map(|x| (x.chain, x.threshold_records))
            .collect();

        let cells: Vec<_> = self.evidence.cells().collect();
        let process_cell = |(cell_id, cell): (u64, &CellEvidence)| {
            let recombinations = cell
                .chains(&self.index)
                .into_iter()
                .filter(|chain| {
                    let threshold = thresholds.get(chain).copied().unwrap_or(1);
                    let v_mappings = cell.segment_mappings(*chain, SegmentKind::V);
                    let j_mappings = cell.segment_mappings(*chain, SegmentKind::J);
                    qualifies_reconstruction_candidate(v_mappings, j_mappings)
                        && cell.reconstruction_records(*chain) >= threshold
                })
                .flat_map(|chain| {
                    process_chain_work(ChainWork {
                        cell,
                        index: &self.index,
                        chain,
                        min_overlap: self.config.min_sequence_overlap,
                    })
                })
                .collect::<Vec<_>>();

            (cell_id, recombinations)
        };

        let mut out = if self.threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .expect("building sc-vdj Rayon pool")
                .install(|| cells.into_par_iter().map(process_cell).collect::<Vec<_>>())
        } else {
            cells.into_iter().map(process_cell).collect::<Vec<_>>()
        };

        out.sort_by_key(|x| x.0);
        out
    }

    pub fn identify(&self) -> Vec<(u64, Vec<Recombination>)> {
        // Library callers keep the unrestricted behavior. The production
        // nelrune-vdj binary explicitly opts into per-locus knee selection.
        let cells: Vec<_> = self.evidence.cells().collect();
        let process_cell = |(cell_id, cell): (u64, &CellEvidence)| {
            let recombinations = cell
                .chains(&self.index)
                .into_iter()
                .flat_map(|chain| {
                    process_chain_work(ChainWork {
                        cell,
                        index: &self.index,
                        chain,
                        min_overlap: self.config.min_sequence_overlap,
                    })
                })
                .collect::<Vec<_>>();
            (cell_id, recombinations)
        };

        let mut out = if self.threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .expect("building sc-vdj Rayon pool")
                .install(|| cells.into_par_iter().map(process_cell).collect::<Vec<_>>())
        } else {
            cells.into_iter().map(process_cell).collect::<Vec<_>>()
        };
        out.sort_by_key(|x| x.0);
        out
    }
}

fn qualifies_reconstruction_candidate(v_mappings: usize, j_mappings: usize) -> bool {
    // Direct V+J evidence remains sufficient.  V-only cells are also allowed
    // when V support is repeated, because fragmented receptor reconstruction can
    // recover a J that was not directly mapped in the first pass.  Requiring at
    // least two V mappings keeps singleton/ambient V hits out of the expensive
    // reconstruction population.
    v_mappings > 0 && (j_mappings > 0 || v_mappings >= 2)
}

fn receptor_knee_threshold(counts_desc: &[usize]) -> usize {
    if counts_desc.is_empty() {
        return 0;
    }
    if counts_desc.len() < 4 || counts_desc[0] == *counts_desc.last().unwrap() {
        return *counts_desc.last().unwrap();
    }

    let n = counts_desc.len();
    let x_max = (n as f64).ln();
    let y_max = (counts_desc[0] as f64).ln();
    let y_min = (*counts_desc.last().unwrap() as f64).ln();
    let y_span = y_max - y_min;
    if x_max <= 0.0 || y_span <= f64::EPSILON {
        return *counts_desc.last().unwrap();
    }

    let mut best_idx = 0usize;
    let mut best_distance = 0.0f64;
    for (i, count) in counts_desc
        .iter()
        .enumerate()
        .skip(1)
        .take(n.saturating_sub(2))
    {
        let x = ((i + 1) as f64).ln() / x_max;
        let y = ((*count as f64).ln() - y_min) / y_span;
        let distance = (1.0 - x) - y;
        if distance > best_distance {
            best_distance = distance;
            best_idx = i;
        }
    }

    // No visible bend: keep every observed cell for this receptor class.
    if best_idx == 0 || best_distance < 0.05 {
        return *counts_desc.last().unwrap();
    }

    // The maximum chord distance lies at the beginning of the low-depth tail.
    // Keep the last depth immediately before that tail, including all ties.
    counts_desc[best_idx - 1]
}

#[cfg(test)]
mod knee_tests {
    use super::{
        qualifies_reconstruction_candidate, receptor_knee_threshold, UnmappedIghSeeds,
        UnmappedReadCandidate, UNMAPPED_IGH_SEED_LEN,
    };
    use crate::cellrep::EvidenceId;
    use crate::index::{Chain, SegmentKind, Strand, VdjIndex, VdjSegment};

    fn test_segment(name: &str, kind: SegmentKind, sequence: &[u8]) -> VdjSegment {
        VdjSegment {
            id: 0,
            name: name.to_string(),
            transcript_id: name.to_string(),
            gene_id: name.to_string(),
            chain: Chain::Igh,
            kind,
            chromosome: "chr12".to_string(),
            start: 0,
            end: sequence.len() as u32,
            strand: Strand::Plus,
            exon_blocks: vec![(0, sequence.len() as u32)],
            coding_start: (kind == SegmentKind::V).then_some(0),
            sequence: sequence.to_vec(),
        }
    }

    #[test]
    fn unmapped_igh_rescue_is_j_gated_and_harvests_v_support() {
        let v = b"ACGTTGCAACCTGATCGTACCGATGCTAGCATGGA";
        let j = b"TTGACCGTATCGGATCCGATGACCTGGA";
        let index = VdjIndex::from_segments(vec![
            test_segment("IGHV-test", SegmentKind::V, v),
            test_segment("IGHJ-test", SegmentKind::J, j),
        ])
        .unwrap();
        let seeds = UnmappedIghSeeds::new(&index);

        let mut read = b"NNNN".to_vec();
        read.extend_from_slice(v);
        read.extend_from_slice(b"GGTACCTAAC");
        read.extend_from_slice(j);
        let candidate = |bases: Vec<u8>| UnmappedReadCandidate {
            cell_id: 7,
            id: EvidenceId { flush: 0, entry: 1 },
            qualities: vec![30; bases.len()],
            bases,
            is_last_in_template: false,
            is_secondary: false,
            is_supplementary: false,
        };

        let (_, rescued) = seeds.rescue_candidate(candidate(read.clone()), &index).unwrap();
        assert!(rescued.mappings.iter().any(|m| {
            index.segment(m.segment_id).unwrap().kind == SegmentKind::V
        }));
        assert!(rescued.mappings.iter().any(|m| {
            index.segment(m.segment_id).unwrap().kind == SegmentKind::J
        }));

        let reverse = crate::index::reverse_complement(&read);
        assert!(seeds.rescue_candidate(candidate(reverse), &index).is_some());

        let mut v_without_j = v.to_vec();
        v_without_j.extend_from_slice(b"AAAAAAAAAAAAAAAAAAAA");
        assert!(seeds
            .rescue_candidate(candidate(v_without_j), &index)
            .is_none());

        let one_j_seed = j[..UNMAPPED_IGH_SEED_LEN].to_vec();
        assert!(seeds
            .rescue_candidate(candidate(one_j_seed), &index)
            .is_none());
    }

    #[test]
    fn reconstruction_candidate_accepts_vj_or_repeated_v_only() {
        assert!(qualifies_reconstruction_candidate(1, 1));
        assert!(qualifies_reconstruction_candidate(2, 0));
        assert!(qualifies_reconstruction_candidate(56, 0));
        assert!(!qualifies_reconstruction_candidate(1, 0));
        assert!(!qualifies_reconstruction_candidate(0, 10));
        assert!(!qualifies_reconstruction_candidate(0, 0));
    }

    #[test]
    fn knee_keeps_sparse_high_depth_head_and_drops_single_read_tail() {
        let counts = vec![100, 50, 10, 2, 1, 1, 1, 1, 1, 1];
        assert_eq!(receptor_knee_threshold(&counts), 2);
    }

    #[test]
    fn knee_keeps_all_when_distribution_has_no_bend() {
        let counts = vec![100, 90, 80, 70, 60, 50, 40, 30, 20, 10];
        assert_eq!(receptor_knee_threshold(&counts), 10);
    }

    #[test]
    fn knee_does_not_delete_tiny_receptor_classes() {
        let counts = vec![5, 2, 1];
        assert_eq!(receptor_knee_threshold(&counts), 1);
    }
}

fn aligned_blocks(record: &bam::Record) -> Vec<(u32, u32)> {
    let mut ref_pos = record.pos().max(0) as u32;
    let mut out = Vec::new();
    let mut block_start = None;
    for c in record.cigar().iter() {
        match *c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                if block_start.is_none() {
                    block_start = Some(ref_pos)
                }
                ref_pos = ref_pos.saturating_add(n)
            }
            Cigar::Del(n) => {
                if block_start.is_none() {
                    block_start = Some(ref_pos)
                }
                ref_pos = ref_pos.saturating_add(n)
            }
            Cigar::RefSkip(n) => {
                if let Some(a) = block_start.take() {
                    if a < ref_pos {
                        out.push((a, ref_pos))
                    }
                }
                ref_pos = ref_pos.saturating_add(n)
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if let Some(a) = block_start {
        if a < ref_pos {
            out.push((a, ref_pos))
        }
    }
    out
}
fn decode_hex_ascii(input: &[u8]) -> Option<String> {
    if input.is_empty() || input.len() % 2 != 0 {
        return None;
    }
    let mut o = Vec::with_capacity(input.len() / 2);
    for x in input.chunks_exact(2) {
        o.push((hex(x[0])? << 4) | hex(x[1])?)
    }
    String::from_utf8(o).ok()
}
fn hex(x: u8) -> Option<u8> {
    match x {
        b'0'..=b'9' => Some(x - b'0'),
        b'a'..=b'f' => Some(x - b'a' + 10),
        b'A'..=b'F' => Some(x - b'A' + 10),
        _ => None,
    }
}
