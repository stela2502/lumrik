//! BAM-derived per-cell receptor evidence.
//!
//! `BamFeatureEvidence` is deliberately short-lived. The runner collects a
//! bounded batch, folds that batch into compact per-cell/per-chain receptor
//! summaries and linkage counters, then drops the raw sequence, qualities and
//! mapper geometry before reading the next batch.

use crate::index::{Chain, SegmentId, SegmentKind, VdjIndex};
use rayon::prelude::*;
use scdata::CellHash;
use std::collections::{HashMap, HashSet};

mod evidence;
mod summary;

pub use evidence::{RawEvidenceDisplay, SummarizedEvidenceDisplay};
pub use summary::{
    summarize_chain_work, ChainSummaryWork, GermlineAnchor, ReceptorSequenceEvidence,
};
use summary::merge_summary_batch;

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
pub struct BamFeatureEvidence {
    pub id: EvidenceId,
    pub sequence: BamFeatureSequenceParts,
    pub mappings: Vec<MapperEvidence>,
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
}

impl CellEvidenceVdj {
    pub fn new() -> Self {
        Self {
            cells: CellHash::new(),
        }
    }

    /// Fold a bounded raw-evidence batch into persistent compact state.
    /// Nothing from `batch` is retained after this function returns.
    pub fn consume_batch(
        &mut self,
        batch: Vec<(u64, BamFeatureEvidence)>,
        index: &VdjIndex,
        min_overlap: usize,
    ) {
        let mut by_cell = HashMap::<u64, Vec<&BamFeatureEvidence>>::new();
        for (cell_id, feature) in &batch {
            by_cell.entry(*cell_id).or_default().push(feature);
        }

        // Cell-local batch compaction is the expensive part and cells are
        // independent. Build one compact delta per cell in parallel, without
        // touching the persistent CellHash from worker threads.
        let deltas: Vec<_> = by_cell
            .into_par_iter()
            .map(|(cell_id, features)| {
                let mut delta = CellEvidence::default();
                delta.consume_compact_metadata(&features, index);

                for chain in Chain::ALL {
                    let chain_features: Vec<_> = features
                        .iter()
                        .copied()
                        .filter(|feature| feature_has_chain(feature, index, chain))
                        .collect();
                    if chain_features.is_empty() {
                        continue;
                    }

                    let batch_summaries = summarize_chain_work(ChainSummaryWork {
                        features: &chain_features,
                        index,
                        chain,
                        min_overlap,
                    });
                    if !batch_summaries.is_empty() {
                        delta.summaries.insert(chain, batch_summaries);
                    }
                }

                (cell_id, delta)
            })
            .collect();

        // Only compact deltas cross back into the persistent scdata-like
        // CellHash. This merge is intentionally serial; the expensive raw
        // read summarization above is already complete and lock-free.
        for (cell_id, delta) in deltas {
            self.cells
                .entry_cell(cell_id)
                .or_default()
                .merge_compact_delta(delta, index, min_overlap);
        }
        // `batch` dies here: sequence, qualities, mappings and ref_blocks are
        // returned to the allocator before the next BAM batch is accumulated.
    }

    pub fn get(&self, cell_id: &u64) -> Option<&CellEvidence> {
        self.cells.get_cell(cell_id)
    }
    pub fn cell_count(&self) -> usize {
        self.cells.cell_count()
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
    ) {
        for (chain, batch_summaries) in delta.summaries.drain() {
            let summaries = self.summaries.entry(chain).or_default();
            merge_summary_batch(summaries, batch_summaries, index, min_overlap);
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
            for mapping in &feature.mappings {
                let Some(segment) = index.segment(mapping.segment_id) else {
                    continue;
                };
                record_chains.insert(segment.chain);
                *self
                    .stats
                    .segment_mappings
                    .entry((segment.chain, segment.kind))
                    .or_default() += 1;
            }
            for chain in record_chains {
                *self.stats.chain_records.entry(chain).or_default() += 1;
            }

            fragments.entry(feature.id).or_default().push(feature);
        }

        self.stats.physical_fragments = self
            .stats
            .physical_fragments
            .saturating_add(fragments.len());

        for fragment in fragments.values() {
            let mut all_kinds = HashSet::<SegmentKind>::new();
            let mut by_chain = HashMap::<Chain, (HashSet<SegmentId>, HashSet<SegmentId>)>::new();

            for feature in fragment {
                for mapping in &feature.mappings {
                    let Some(segment) = index.segment(mapping.segment_id) else {
                        continue;
                    };
                    all_kinds.insert(segment.kind);
                    let (receptor, constant) = by_chain.entry(segment.chain).or_default();
                    match segment.kind {
                        SegmentKind::V | SegmentKind::J => {
                            receptor.insert(mapping.segment_id);
                        }
                        SegmentKind::C => {
                            constant.insert(mapping.segment_id);
                        }
                        SegmentKind::D => {}
                    }
                }
            }

            if all_kinds.len() > 1 {
                self.stats.linked_fragments = self.stats.linked_fragments.saturating_add(1);
            }

            for (chain, (receptor, constant)) in by_chain {
                if receptor.is_empty() || constant.is_empty() {
                    continue;
                }
                let mut receptor_segments: Vec<_> = receptor.into_iter().collect();
                let mut constant_segments: Vec<_> = constant.into_iter().collect();
                receptor_segments.sort_unstable();
                constant_segments.sort_unstable();
                let signature = FragmentLinkSignature {
                    chain,
                    receptor_segments,
                    constant_segments,
                };
                let n = self.fragment_link_support.entry(signature).or_default();
                *n = n.saturating_add(1);
            }
        }
    }
}

fn feature_has_chain(feature: &BamFeatureEvidence, index: &VdjIndex, chain: Chain) -> bool {
    feature.mappings.iter().any(|mapping| {
        index
            .segment(mapping.segment_id)
            .is_some_and(|segment| segment.chain == chain)
    })
}
