//! BAM-derived per-cell receptor evidence.
//!
//! `BamFeatureEvidence` is deliberately short-lived. The runner collects a
//! bounded batch, folds that batch into compact per-cell/per-chain receptor
//! summaries and linkage counters, then drops the raw sequence, qualities and
//! mapper geometry before reading the next batch.

use crate::index::{Chain, SegmentId, SegmentKind, VdjIndex};
use rayon::prelude::*;
use run_stats::RunStats;
use rust_htslib::bam::record::Cigar;
use scdata::CellHash;
use std::collections::{HashMap, HashSet};
use std::time::Instant;

mod evidence;
mod summary;

pub use evidence::{RawEvidenceDisplay, SummarizedEvidenceDisplay};
use summary::{flush_refinement_buffers, merge_completed_models, ChainSummaryWork};
pub use summary::{GermlineAnchor, ReceptorSequenceEvidence};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EvidenceId {
    pub flush: u32,
    pub entry: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequencePart {
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
}
impl SequencePart {
    pub fn new(bases: Vec<u8>, qualities: Vec<u8>) -> Self {
        Self { bases, qualities }
    }
    pub fn uniform(bases: &[u8], q: u8) -> Self {
        Self {
            bases: bases.to_vec(),
            qualities: vec![q; bases.len()],
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BamFeatureSequenceParts {
    pub r1: Option<SequencePart>,
    pub r2: Option<SequencePart>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignmentGeometry {
    pub tid: i32,
    pub start: u32,
    pub end: u32,
    pub is_reverse: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub mapq: u8,
    pub ref_blocks: Vec<(u32, u32)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MapperEvidence {
    pub segment_id: SegmentId,
    pub alignment: AlignmentGeometry,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawBamFeatureEvidence {
    pub id: EvidenceId,
    pub umi: Option<Vec<u8>>,
    pub tid: i32,
    pub pos: u32,
    pub cigar: Vec<Cigar>,
    pub bases: Vec<u8>,
    pub qualities: Vec<u8>,
    pub is_last_in_template: bool,
    pub is_reverse: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub mapq: u8,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RawCellEvidence {
    pub records: Vec<RawBamFeatureEvidence>,
}

const UMI_REGION_JOIN_GAP: u32 = 8;

#[derive(Debug, Clone, Default)]
struct UmiRegionState {
    regions: HashMap<Vec<u8>, Vec<(u32, u32)>>,
}

impl UmiRegionState {
    #[inline]
    fn read_span(read: &RawBamFeatureEvidence) -> (u32, u32) {
        let start = read.pos;
        let end = start.saturating_add(read.bases.len() as u32);
        (start, end)
    }

    /// Add one observed UMI span. Regions separated by at most 8 bp are treated
    /// as connected informative territory; larger holes remain explicit.
    fn add_span(regions: &mut Vec<(u32, u32)>, start: u32, end: u32) {
        let mut merged_start = start;
        let mut merged_end = end;
        regions.retain(|&(known_start, known_end)| {
            let connected = known_start <= merged_end.saturating_add(UMI_REGION_JOIN_GAP)
                && merged_start <= known_end.saturating_add(UMI_REGION_JOIN_GAP);
            if connected {
                merged_start = merged_start.min(known_start);
                merged_end = merged_end.max(known_end);
                false
            } else {
                true
            }
        });
        regions.push((merged_start, merged_end));
    }

    fn absorb_read(&mut self, read: &RawBamFeatureEvidence) {
        let Some(umi) = read.umi.as_ref() else {
            return;
        };
        let (start, end) = Self::read_span(read);
        Self::add_span(self.regions.entry(umi.clone()).or_default(), start, end);
    }
}

impl RawCellEvidence {
    /// Copy-on-first-touch UMI reduction against persistent genomic territory.
    /// Only a touched UMI's tiny region list is cloned. Reads wholly contained
    /// in already observed territory are discarded before receptor/model work.
    /// Reads that extend, bridge <=8 bp, or create a new region are retained.
    fn retain_umi_region_extensions(&mut self, persistent: Option<&UmiRegionState>) {
        let mut touched = HashMap::<Vec<u8>, Vec<(u32, u32)>>::new();

        self.records.retain(|read| {
            let Some(umi) = read.umi.as_ref() else {
                return true;
            };
            let (start, end) = UmiRegionState::read_span(read);
            let regions = touched.entry(umi.clone()).or_insert_with(|| {
                persistent
                    .and_then(|state| state.regions.get(umi))
                    .cloned()
                    .unwrap_or_default()
            });

            if regions
                .iter()
                .any(|&(known_start, known_end)| known_start <= start && end <= known_end)
            {
                return false;
            }

            UmiRegionState::add_span(regions, start, end);
            true
        });
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BamFeatureEvidence {
    pub id: EvidenceId,
    pub sequence: BamFeatureSequenceParts,
    /// Corrected molecule barcode (BAM UB tag), retained so receptor support can
    /// be reported as independent UMIs rather than amplified/read-level evidence.
    pub umi: Option<Vec<u8>>,
    pub mappings: Vec<MapperEvidence>,
    /// Constant-gene spans covered outside annotated exons. These are retained
    /// only as compact intronic-state evidence and never seed reconstruction.
    pub intronic_constant_segments: Vec<SegmentId>,
}

/// Compact aggregate of one kind of direct physical-fragment receptor -> C link.
/// Many fragments with the same segment signature collapse into one counter.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FragmentLinkSignature {
    pub chain: Chain,
    pub receptor_segments: Vec<SegmentId>,
    pub constant_segments: Vec<SegmentId>,
}

#[derive(Debug, Clone, Default)]
pub struct CellEvidenceStats {
    accepted_records: usize,
    physical_fragments: usize,
    single_segment_records: usize,
    multi_segment_records: usize,
    linked_fragments: usize,
    chain_records: HashMap<Chain, usize>,
    reconstruction_records: HashMap<Chain, usize>,
    intronic_constant_records: HashMap<Chain, usize>,
    intronic_constant_segments: HashMap<SegmentId, usize>,
    segment_mappings: HashMap<(Chain, SegmentKind), usize>,
}

#[derive(Debug, Clone, Default)]
pub struct CellEvidence {
    summaries: HashMap<Chain, Vec<ReceptorSequenceEvidence>>,
    fragment_link_support: HashMap<FragmentLinkSignature, u32>,
    stats: CellEvidenceStats,
}

#[derive(Debug, Default)]
pub struct CellEvidenceVdj {
    cells: CellHash<CellEvidence>,
    /// Incrementally maintained live model count. Batch merges update this from
    /// local deltas; the hot path never scans the persistent CellHash to count models.
    live_model_count: usize,
    /// Number of raw evidence batches consumed; used only for compact telemetry.
    batch_number: usize,
    /// Persistent per-cell/per-UMI connected genomic territory. This is a cheap
    /// pre-model memory: repeated PCR observations that add no new territory can
    /// be rejected before receptor interpretation.
    umi_regions: CellHash<UmiRegionState>,
}

impl CellEvidenceVdj {
    pub fn new() -> Self {
        Self {
            cells: CellHash::new(),
            live_model_count: 0,
            batch_number: 0,
            umi_regions: CellHash::new(),
        }
    }

    /// Interpret one minimal CellHash BAM batch entirely against fresh batch-local
    /// model state, then reconcile the resulting compact model tree with persistent
    /// evidence in one separate bucket-parallel merge. Raw reads never interrogate
    /// the progressively growing historical model tree.
    pub fn consume_raw_batch(
        &mut self,
        batch: CellHash<RawCellEvidence>,
        contig_names: &[String],
        index: &VdjIndex,
        min_overlap: usize,
    ) -> (usize, RunStats) {
        self.batch_number = self.batch_number.saturating_add(1);
        let batch_number = self.batch_number;
        let input_records: usize = batch
            .buckets()
            .iter()
            .flat_map(|bucket| bucket.values())
            .map(|cell| cell.records.len())
            .sum();

        // Phase 1: reject same-UMI observations that add no genomic territory.
        // Persistent state is read-only here. Each worker clones only the tiny
        // region list for a UMI on first touch, then mutates that local copy.
        let umi_reduce_started = Instant::now();
        let reduced_batch =
            batch.par_map_owned_with(&self.umi_regions, |raw_bucket, region_bucket| {
                let mut reduced_bucket = HashMap::with_capacity(raw_bucket.len());
                for (cell_id, mut raw_cell) in raw_bucket {
                    raw_cell.retain_umi_region_extensions(region_bucket.get(&cell_id));
                    if !raw_cell.records.is_empty() {
                        reduced_bucket.insert(cell_id, raw_cell);
                    }
                }
                reduced_bucket
            });

        let umi_reduce_elapsed = umi_reduce_started.elapsed();
        let reduced_records: usize = reduced_batch
            .buckets()
            .iter()
            .flat_map(|bucket| bucket.values())
            .map(|cell| cell.records.len())
            .sum();
        let umi_contained_skips = input_records.saturating_sub(reduced_records);

        // Persist only information that survived the local reduction. Replaying
        // those retained spans against the global state reproduces the worker's
        // copy-on-first-touch result without sharing mutable state in Rayon.
        let umi_merge_started = Instant::now();
        let _ = self.umi_regions.par_merge_with_mapping(
            &reduced_batch,
            |region_bucket, raw_bucket, _local_stats| {
                for (&cell_id, raw_cell) in raw_bucket {
                    let state = region_bucket.entry(cell_id).or_default();
                    for read in &raw_cell.records {
                        state.absorb_read(read);
                    }
                }
            },
        );

        let umi_merge_elapsed = umi_merge_started.elapsed();

        // Phase 2: build a completely independent model tree from only reads that
        // add molecule territory. No historical depth/saturation analysis is used.
        let model_build_started = Instant::now();
        let empty_context = CellHash::<()>::new();
        let (batch_models, mut batch_stats) = reduced_batch.par_map_owned_with_stats(
            &empty_context,
            |raw_bucket, _unused_bucket, local_stats| {
                let mut model_bucket = HashMap::with_capacity(raw_bucket.len());
                for (cell_id, raw_cell) in raw_bucket {
                    let features: Vec<BamFeatureEvidence> = raw_cell
                        .records
                        .into_iter()
                        .filter_map(|raw| interpret_raw_feature(raw, contig_names, index))
                        .collect();
                    if features.is_empty() {
                        continue;
                    }

                    let feature_refs: Vec<_> = features.iter().collect();
                    let mut cell = CellEvidence::default();
                    cell.consume_compact_metadata(&feature_refs, index);

                    for chain in Chain::ALL {
                        let chain_features: Vec<_> = feature_refs
                            .iter()
                            .copied()
                            .filter(|feature| {
                                feature_has_rearrangement_segment(feature, index, chain)
                            })
                            .collect();
                        if chain_features.is_empty() {
                            continue;
                        }
                        let summaries = ChainSummaryWork {
                            features: &chain_features,
                            index,
                            chain,
                            min_overlap,
                        }
                        .summarize(local_stats);
                        if !summaries.is_empty() {
                            cell.summaries.insert(chain, summaries);
                        }
                    }

                    model_bucket.insert(cell_id, cell);
                }
                model_bucket
            },
        );

        let model_build_elapsed = model_build_started.elapsed();
        let batch_model_count: usize = batch_models
            .buckets()
            .iter()
            .flat_map(|bucket| bucket.values())
            .map(CellEvidence::summary_count)
            .sum();

        // Phase 3: merge compact batch-local models into the persistent tree once.
        // UMI provenance lives on ReceptorSequenceEvidence and is merged by the
        // existing summary merge machinery, so local collapsing does not turn PCR
        // copies into independent molecule support.
        let model_merge_started = Instant::now();
        let merge_mapping = self.cells.par_merge_owned_with_mapping(
            batch_models,
            |self_bucket, model_bucket, mapping| {
                for (cell_id, delta) in model_bucket {
                    mapping.ok_reads = mapping
                        .ok_reads
                        .saturating_add(delta.stats.accepted_records);
                    let cell = self_bucket.entry(cell_id).or_default();
                    let before = cell.summary_count();
                    cell.merge_compact_delta(delta, index, min_overlap, mapping);
                    let after = cell.summary_count();
                    if after >= before {
                        mapping.report_n("vdj.batch.models_added", after - before);
                    } else {
                        mapping.report_n("vdj.batch.models_removed", before - after);
                    }
                }
            },
        );

        let model_merge_elapsed = model_merge_started.elapsed();
        let models_merged = merge_mapping.get_issue_count("vdj.model_merge.merged");
        let models_inserted = merge_mapping.get_issue_count("vdj.model_merge.inserted");
        let merge_accounted = models_merged.saturating_add(models_inserted);
        let merge_us_per_merged = if models_merged == 0 {
            0.0
        } else {
            model_merge_elapsed.as_secs_f64() * 1_000_000.0 / models_merged as f64
        };
        debug_assert_eq!(
            batch_model_count, merge_accounted,
            "every completed batch model must be either inserted or merged"
        );

        let models_added = merge_mapping.get_issue_count("vdj.batch.models_added");
        let models_removed = merge_mapping.get_issue_count("vdj.batch.models_removed");
        self.live_model_count = self
            .live_model_count
            .saturating_add(models_added)
            .saturating_sub(models_removed);
        batch_stats.merge(&merge_mapping);

        batch_stats.reads_log.insert(
            "vdj.umi_region.contained_skips".to_string(),
            umi_contained_skips,
        );
        batch_stats.add_timer_duration("vdj.batch.umi_reduce", umi_reduce_elapsed);
        batch_stats.add_timer_duration("vdj.batch.umi_persist", umi_merge_elapsed);
        batch_stats.add_timer_duration("vdj.batch.model_build", model_build_elapsed);
        batch_stats.add_timer_duration("vdj.batch.model_merge", model_merge_elapsed);

        if batch_number == 1 {
            eprintln!(
                "[nelrune vdj] batch performance\n\
 batch |      raw | retained |  skipped | incoming | inserted |  merged | build_s | merge_s | us/merge\n\
-------+----------+----------+----------+----------+----------+---------+---------+---------+---------"
            );
        }
        eprintln!(
            "{:>6} | {:>8} | {:>8} | {:>8} | {:>8} | {:>8} | {:>7} | {:>7.3} | {:>7.3} | {:>8.1}",
            batch_number,
            input_records,
            reduced_records,
            umi_contained_skips,
            batch_model_count,
            models_inserted,
            models_merged,
            model_build_elapsed.as_secs_f64(),
            model_merge_elapsed.as_secs_f64(),
            merge_us_per_merged,
        );

        (merge_mapping.ok_reads, batch_stats)
    }

    /// Existing interpreted-evidence entry point retained for unmapped rescue.
    pub fn consume_batch(
        &mut self,
        batch: Vec<(u64, BamFeatureEvidence)>,
        index: &VdjIndex,
        min_overlap: usize,
        run_stats: &mut RunStats,
    ) {
        let mut by_cell = CellHash::<Vec<BamFeatureEvidence>>::new();
        for (cell_id, feature) in batch {
            by_cell.entry_cell(cell_id).or_default().push(feature);
        }
        let batch_stats = self.cells.par_merge_owned_with_mapping(
            by_cell,
            |self_bucket, other_bucket, mapping| {
                for (cell_id, features) in other_bucket {
                    let refs: Vec<_> = features.iter().collect();
                    let mut delta = CellEvidence::default();
                    delta.consume_compact_metadata(&refs, index);
                    for chain in Chain::ALL {
                        let chain_features: Vec<_> = refs
                            .iter()
                            .copied()
                            .filter(|feature| {
                                feature_has_rearrangement_segment(feature, index, chain)
                            })
                            .collect();
                        if chain_features.is_empty() {
                            continue;
                        }
                        let summaries = ChainSummaryWork {
                            features: &chain_features,
                            index,
                            chain,
                            min_overlap,
                        }
                        .summarize(mapping);
                        if !summaries.is_empty() {
                            delta.summaries.insert(chain, summaries);
                        }
                    }
                    mapping.ok_reads = mapping
                        .ok_reads
                        .saturating_add(delta.stats.accepted_records);
                    let cell = self_bucket.entry(cell_id).or_default();
                    let before = cell.summary_count();
                    cell.merge_compact_delta(delta, index, min_overlap, mapping);
                    let after = cell.summary_count();
                    if after >= before {
                        mapping.report_n("vdj.batch.models_added", after - before);
                    } else {
                        mapping.report_n("vdj.batch.models_removed", before - after);
                    }
                }
            },
        );
        self.live_model_count = self
            .live_model_count
            .saturating_add(batch_stats.get_issue_count("vdj.batch.models_added"))
            .saturating_sub(batch_stats.get_issue_count("vdj.batch.models_removed"));
        run_stats.merge(&batch_stats);
    }

    /// Finalize every partial per-model refinement buffer after the complete
    /// BAM has been consumed. This is intentionally separate from consume_batch:
    /// BAM/Rayon transport boundaries must not force 1..19-read refinements.
    pub fn flush_pending_refinement(
        &mut self,
        index: &VdjIndex,
        min_overlap: usize,
        run_stats: &mut RunStats,
    ) {
        let stats = self
            .cells
            .buckets_mut()
            .par_iter_mut()
            .map(|bucket| {
                let mut local_stats = RunStats::new(None, 0.0, 0);
                for cell in bucket.values_mut() {
                    for summaries in cell.summaries.values_mut() {
                        flush_refinement_buffers(
                            summaries,
                            index,
                            min_overlap,
                            true,
                            &mut local_stats,
                        );
                    }
                }
                local_stats
            })
            .reduce(
                || RunStats::new(None, 0.0, 0),
                |mut stats, local_stats| {
                    stats.merge(&local_stats);
                    stats
                },
            );
        run_stats.merge(&stats);
    }

    pub fn get(&self, cell_id: &u64) -> Option<&CellEvidence> {
        self.cells.get_cell(cell_id)
    }

    pub fn cell_count(&self) -> usize {
        self.cells.cell_count()
    }

    /// O(1) live model count maintained by batch-local delta merges.
    pub fn live_model_count(&self) -> usize {
        self.live_model_count
    }
    pub fn cells(&self) -> impl Iterator<Item = (u64, &CellEvidence)> {
        self.cells
            .buckets()
            .iter()
            .flat_map(|bucket| bucket.iter().map(|(id, e)| (*id, e)))
    }
    pub fn into_cells(self) -> impl Iterator<Item = (u64, CellEvidence)> {
        self.cells.into_iter().flat_map(|bucket| bucket.into_iter())
    }
}

impl CellEvidence {
    pub fn chains(&self, _index: &VdjIndex) -> Vec<Chain> {
        let mut out: Vec<_> = self.summaries.keys().copied().collect();
        out.sort();
        out
    }

    pub fn summaries_for_chain(&self, chain: Chain) -> &[ReceptorSequenceEvidence] {
        self.summaries.get(&chain).map(Vec::as_slice).unwrap_or(&[])
    }

    pub fn summary_count(&self) -> usize {
        self.summaries.values().map(Vec::len).sum()
    }

    pub fn fragment_link_support(&self) -> &HashMap<FragmentLinkSignature, u32> {
        &self.fragment_link_support
    }

    pub fn accepted_records(&self) -> usize {
        self.stats.accepted_records
    }
    pub fn physical_fragments(&self) -> usize {
        self.stats.physical_fragments
    }
    pub fn single_segment_records(&self) -> usize {
        self.stats.single_segment_records
    }
    pub fn multi_segment_records(&self) -> usize {
        self.stats.multi_segment_records
    }
    pub fn linked_fragments(&self) -> usize {
        self.stats.linked_fragments
    }
    pub fn chain_records(&self, chain: Chain) -> usize {
        self.stats.chain_records.get(&chain).copied().unwrap_or(0)
    }

    /// Records carrying V or J evidence for this receptor class. Constant-only
    /// transcription is deliberately excluded so noisy/sterile C-region signal
    /// cannot decide which cell/locus combinations enter reconstruction.
    pub fn reconstruction_records(&self, chain: Chain) -> usize {
        self.stats
            .reconstruction_records
            .get(&chain)
            .copied()
            .unwrap_or(0)
    }

    pub fn intronic_constant_records(&self, chain: Chain) -> usize {
        self.stats
            .intronic_constant_records
            .get(&chain)
            .copied()
            .unwrap_or(0)
    }

    pub fn intronic_constant_segment_records(&self, segment: SegmentId) -> usize {
        self.stats
            .intronic_constant_segments
            .get(&segment)
            .copied()
            .unwrap_or(0)
    }

    /// Number of physical fragments with both J and constant-region evidence.
    /// This is stronger evidence for an expressed rearranged receptor than an
    /// unlinked pile of C-region mappings.
    pub fn j_constant_linked_fragments(&self, index: &VdjIndex) -> usize {
        self.fragment_link_support
            .iter()
            .filter_map(|(signature, &count)| {
                let has_j = signature.receptor_segments.iter().any(|id| {
                    index
                        .segment(*id)
                        .is_some_and(|segment| segment.kind == SegmentKind::J)
                });
                (has_j && !signature.constant_segments.is_empty()).then_some(count as usize)
            })
            .sum()
    }
    pub fn segment_mappings(&self, chain: Chain, kind: SegmentKind) -> usize {
        self.stats
            .segment_mappings
            .get(&(chain, kind))
            .copied()
            .unwrap_or(0)
    }

    fn merge_compact_delta(
        &mut self,
        mut delta: CellEvidence,
        index: &VdjIndex,
        min_overlap: usize,
        stats: &mut RunStats,
    ) {
        for (chain, batch_summaries) in delta.summaries.drain() {
            let summaries = self.summaries.entry(chain).or_default();
            merge_completed_models(summaries, batch_summaries, index, min_overlap, stats, true);
        }

        for (signature, count) in delta.fragment_link_support.drain() {
            let support = self.fragment_link_support.entry(signature).or_default();
            *support = support.saturating_add(count);
        }

        self.stats.accepted_records = self
            .stats
            .accepted_records
            .saturating_add(delta.stats.accepted_records);
        self.stats.physical_fragments = self
            .stats
            .physical_fragments
            .saturating_add(delta.stats.physical_fragments);
        self.stats.single_segment_records = self
            .stats
            .single_segment_records
            .saturating_add(delta.stats.single_segment_records);
        self.stats.multi_segment_records = self
            .stats
            .multi_segment_records
            .saturating_add(delta.stats.multi_segment_records);
        self.stats.linked_fragments = self
            .stats
            .linked_fragments
            .saturating_add(delta.stats.linked_fragments);

        for (chain, count) in delta.stats.chain_records.drain() {
            let n = self.stats.chain_records.entry(chain).or_default();
            *n = n.saturating_add(count);
        }
        for (chain, count) in delta.stats.reconstruction_records.drain() {
            let n = self.stats.reconstruction_records.entry(chain).or_default();
            *n = n.saturating_add(count);
        }
        for (chain, count) in delta.stats.intronic_constant_records.drain() {
            let n = self
                .stats
                .intronic_constant_records
                .entry(chain)
                .or_default();
            *n = n.saturating_add(count);
        }
        for (segment, count) in delta.stats.intronic_constant_segments.drain() {
            let n = self
                .stats
                .intronic_constant_segments
                .entry(segment)
                .or_default();
            *n = n.saturating_add(count);
        }
        for (key, count) in delta.stats.segment_mappings.drain() {
            let n = self.stats.segment_mappings.entry(key).or_default();
            *n = n.saturating_add(count);
        }
    }

    fn consume_compact_metadata(&mut self, features: &[&BamFeatureEvidence], index: &VdjIndex) {
        self.stats.accepted_records = self.stats.accepted_records.saturating_add(features.len());

        let mut fragments = HashMap::<EvidenceId, Vec<&BamFeatureEvidence>>::new();
        for feature in features {
            if feature.mappings.len() <= 1 {
                self.stats.single_segment_records =
                    self.stats.single_segment_records.saturating_add(1);
            } else {
                self.stats.multi_segment_records =
                    self.stats.multi_segment_records.saturating_add(1);
            }

            let mut record_chains = HashSet::<Chain>::new();
            let mut reconstruction_chains = HashSet::<Chain>::new();
            for mapping in &feature.mappings {
                let Some(segment) = index.segment(mapping.segment_id) else {
                    continue;
                };
                record_chains.insert(segment.chain);
                if matches!(segment.kind, SegmentKind::V | SegmentKind::J) {
                    reconstruction_chains.insert(segment.chain);
                }
                *self
                    .stats
                    .segment_mappings
                    .entry((segment.chain, segment.kind))
                    .or_default() += 1;
            }
            for chain in record_chains {
                *self.stats.chain_records.entry(chain).or_default() += 1;
            }
            for chain in reconstruction_chains {
                *self.stats.reconstruction_records.entry(chain).or_default() += 1;
            }

            let mut intronic_chains = HashSet::<Chain>::new();
            for &segment_id in &feature.intronic_constant_segments {
                let Some(segment) = index.segment(segment_id) else {
                    continue;
                };
                intronic_chains.insert(segment.chain);
                *self
                    .stats
                    .intronic_constant_segments
                    .entry(segment_id)
                    .or_default() += 1;
            }
            for chain in intronic_chains {
                *self
                    .stats
                    .intronic_constant_records
                    .entry(chain)
                    .or_default() += 1;
            }

            fragments.entry(feature.id).or_default().push(feature);
        }

        self.stats.physical_fragments = self
            .stats
            .physical_fragments
            .saturating_add(fragments.len());

        for fragment in fragments.values() {
            let mut all_kinds = HashSet::<SegmentKind>::new();
            let mut validated_links = HashSet::<FragmentLinkSignature>::new();

            for feature in fragment {
                for mapping in &feature.mappings {
                    if let Some(segment) = index.segment(mapping.segment_id) {
                        all_kinds.insert(segment.kind);
                    }
                }

                for chain in Chain::ALL {
                    for (j, c) in splice_supported_jc_links(feature, index, chain, 8) {
                        validated_links.insert(FragmentLinkSignature {
                            chain,
                            receptor_segments: vec![j],
                            constant_segments: vec![c],
                        });
                    }
                }
            }

            if all_kinds.len() > 1 {
                self.stats.linked_fragments = self.stats.linked_fragments.saturating_add(1);
            }

            for signature in validated_links {
                let n = self.fragment_link_support.entry(signature).or_default();
                *n = n.saturating_add(1);
            }
        }
    }
}

fn interpret_raw_feature(
    raw: RawBamFeatureEvidence,
    contig_names: &[String],
    index: &VdjIndex,
) -> Option<BamFeatureEvidence> {
    let chr = contig_names.get(raw.tid as usize)?;
    let blocks = aligned_blocks_from_cigar(raw.pos, &raw.cigar);
    let segment_ids = index.overlapping(chr, &blocks);
    let intronic_constant_segments = index.intronic_constant_overlapping(chr, &blocks, 8);
    let rearrangement_chains: HashSet<_> = segment_ids
        .iter()
        .filter_map(|id| index.segment(*id))
        .filter(|segment| segment.kind != SegmentKind::C)
        .map(|segment| segment.chain)
        .collect();
    if rearrangement_chains.is_empty() && intronic_constant_segments.is_empty() {
        return None;
    }

    let geometry = AlignmentGeometry {
        tid: raw.tid,
        start: blocks.first().map_or(raw.pos, |x| x.0),
        end: blocks.last().map_or(raw.pos, |x| x.1),
        is_reverse: raw.is_reverse,
        is_secondary: raw.is_secondary,
        is_supplementary: raw.is_supplementary,
        mapq: raw.mapq,
        ref_blocks: blocks,
    };
    let mappings = segment_ids
        .into_iter()
        .filter(|segment_id| {
            index.segment(*segment_id).is_some_and(|segment| {
                segment.kind == SegmentKind::C || rearrangement_chains.contains(&segment.chain)
            })
        })
        .map(|segment_id| MapperEvidence {
            segment_id,
            alignment: geometry.clone(),
        })
        .collect();

    let mut sequence = BamFeatureSequenceParts::default();
    if !rearrangement_chains.is_empty() {
        let part = SequencePart {
            bases: raw.bases,
            qualities: raw.qualities,
        };
        if raw.is_last_in_template {
            sequence.r2 = Some(part);
        } else {
            sequence.r1 = Some(part);
        }
    }
    Some(BamFeatureEvidence {
        id: raw.id,
        sequence,
        umi: raw.umi,
        mappings,
        intronic_constant_segments,
    })
}

fn aligned_blocks_from_cigar(start: u32, cigar: &[Cigar]) -> Vec<(u32, u32)> {
    let mut ref_pos = start;
    let mut out = Vec::new();
    let mut block_start = None;
    for c in cigar {
        match *c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                if block_start.is_none() {
                    block_start = Some(ref_pos);
                }
                ref_pos = ref_pos.saturating_add(n);
            }
            Cigar::Del(n) => {
                if block_start.is_none() {
                    block_start = Some(ref_pos);
                }
                ref_pos = ref_pos.saturating_add(n);
            }
            Cigar::RefSkip(n) => {
                if let Some(a) = block_start.take() {
                    if a < ref_pos {
                        out.push((a, ref_pos));
                    }
                }
                ref_pos = ref_pos.saturating_add(n);
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if let Some(a) = block_start {
        if a < ref_pos {
            out.push((a, ref_pos));
        }
    }
    out
}

fn feature_has_rearrangement_segment(
    feature: &BamFeatureEvidence,
    index: &VdjIndex,
    chain: Chain,
) -> bool {
    feature.mappings.iter().any(|mapping| {
        index
            .segment(mapping.segment_id)
            .is_some_and(|segment| segment.chain == chain && segment.kind != SegmentKind::C)
    })
}

fn splice_supported_jc_links(
    feature: &BamFeatureEvidence,
    index: &VdjIndex,
    chain: Chain,
    min_anchor: u32,
) -> Vec<(SegmentId, SegmentId)> {
    let mut js = Vec::new();
    let mut cs = Vec::new();
    for mapping in &feature.mappings {
        let Some(segment) = index.segment(mapping.segment_id) else {
            continue;
        };
        if segment.chain != chain {
            continue;
        }
        match segment.kind {
            SegmentKind::J => js.push(mapping.segment_id),
            SegmentKind::C => cs.push(mapping.segment_id),
            _ => {}
        }
    }

    let Some(geometry) = feature.mappings.first().map(|m| &m.alignment) else {
        return Vec::new();
    };
    if geometry.ref_blocks.len() < 2 {
        return Vec::new();
    }

    let mut out = Vec::new();
    for j in js {
        let Some(jseg) = index.segment(j) else {
            continue;
        };
        for &c in &cs {
            let Some(cseg) = index.segment(c) else {
                continue;
            };
            if jseg.chromosome != cseg.chromosome || jseg.strand != cseg.strand {
                continue;
            }
            if has_splice_bridge(geometry, jseg, cseg, min_anchor) {
                out.push((j, c));
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn has_splice_bridge(
    geometry: &AlignmentGeometry,
    j: &crate::index::VdjSegment,
    c: &crate::index::VdjSegment,
    min_anchor: u32,
) -> bool {
    let block_overlap = |block: (u32, u32), exons: &[(u32, u32)]| -> u32 {
        exons
            .iter()
            .map(|&(a, b)| block.1.min(b).saturating_sub(block.0.max(a)))
            .sum()
    };

    for (ji, &jb) in geometry.ref_blocks.iter().enumerate() {
        if block_overlap(jb, &j.exon_blocks) < min_anchor {
            continue;
        }
        for (ci, &cb) in geometry.ref_blocks.iter().enumerate() {
            if ji == ci || block_overlap(cb, &c.exon_blocks) < min_anchor {
                continue;
            }
            let transcript_order_ok = match j.strand {
                crate::index::Strand::Plus => ji < ci,
                crate::index::Strand::Minus => ji > ci,
            };
            if transcript_order_ok {
                return true;
            }
        }
    }
    false
}
