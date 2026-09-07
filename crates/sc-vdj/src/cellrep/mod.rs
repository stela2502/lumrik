//! BAM-derived per-cell receptor evidence.
//!
//! A `BamFeatureEvidence` is the indivisible input object: sequence and mapper
//! geometry stay together.  We never manufacture a later biological identity
//! from QNAME or UMI strings.

use crate::index::{Chain, SegmentId, VdjIndex};
use scdata::CellHash;
use std::collections::HashSet;

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

#[derive(Debug, Clone, Default)]
pub struct CellEvidence {
    pub features: Vec<BamFeatureEvidence>,
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
    pub fn push(&mut self, cell_id: u64, evidence: BamFeatureEvidence) {
        self.cells
            .entry_cell(cell_id)
            .or_default()
            .features
            .push(evidence)
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
    pub fn chains(&self, index: &VdjIndex) -> Vec<Chain> {
        let mut set = HashSet::new();
        for f in &self.features {
            for m in &f.mappings {
                if let Some(s) = index.segment(m.segment_id) {
                    set.insert(s.chain);
                }
            }
        }
        let mut out: Vec<_> = set.into_iter().collect();
        out.sort();
        out
    }
    pub fn features_for_chain<'a>(
        &'a self,
        index: &'a VdjIndex,
        chain: Chain,
    ) -> Vec<&'a BamFeatureEvidence> {
        self.features
            .iter()
            .filter(|f| {
                f.mappings.iter().any(|m| {
                    index
                        .segment(m.segment_id)
                        .is_some_and(|s| s.chain == chain)
                })
            })
            .collect()
    }
}

mod evidence;
mod summary;

pub use evidence::{RawEvidenceDisplay, SummarizedEvidenceDisplay};
pub use summary::{
    summarize_chain_work, ChainSummaryWork, GermlineAnchor, ReceptorSequenceEvidence,
};
