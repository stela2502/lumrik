use crate::sequence::{merge_sequences_any_orientation, reverse_complement};
use crate::types::{Chain, SegmentKind};
use std::collections::HashMap;

pub(crate) const DEFAULT_SEQUENCE_OVERLAP: usize = 20;

#[derive(Debug, Clone, Default)]
pub(crate) struct BamFeatureSequenceParts {
    /// Stable ID of this provisional sequence entry. Mapper components and
    /// finalized contigs keep this ID so V/D/J/C evidence can be linked back
    /// to the same original BAM feature sequence.
    pub(crate) sequence_id: u64,
    pub(crate) r1: Option<(Vec<u8>, Vec<u8>)>,
    pub(crate) r2: Option<(Vec<u8>, Vec<u8>)>,
}

impl BamFeatureSequenceParts {
    pub(crate) fn sequence(&self, read2: bool) -> Option<(&[u8], &[u8])> {
        let part = if read2 { self.r2.as_ref() } else { self.r1.as_ref() }?;
        Some((&part.0, &part.1))
    }
}

#[derive(Debug, Clone)]
pub(crate) struct MappedSegment {
    pub(crate) segment_index: usize,
    pub(crate) chain: Chain,
    pub(crate) kind: SegmentKind,
}

#[derive(Debug, Clone)]
pub(crate) struct InitialMapperResults {
    pub(crate) segments: Vec<MappedSegment>,
    pub(crate) read2: bool,
    pub(crate) ref_blocks: Vec<(u32, u32)>,
    pub(crate) mapq: u8,
    pub(crate) is_supplementary: bool,
}

#[derive(Debug, Default)]
pub(crate) struct CellEvidence {
    pub(crate) cell: String,
    pub(crate) provisional: Vec<(BamFeatureSequenceParts, Vec<InitialMapperResults>)>,
    pub(crate) final_results: HashMap<Chain, Vec<CellRecombinationResults>>,
}

#[derive(Debug, Clone)]
pub(crate) struct CellRecombinationResults {
    base_counts: [Vec<u16>; 4],
    base_max_qual: [Vec<u8>; 4],
    pub(crate) ref_blocks: Vec<(u32, u32)>,
    pub(crate) mapq: u8,
    pub(crate) has_primary: bool,
    pub(crate) supporting_sequence_ids: Vec<u64>,
    pub(crate) mapped_segment_indices: Vec<usize>,
}

impl Default for CellRecombinationResults {
    fn default() -> Self {
        Self {
            base_counts: std::array::from_fn(|_| Vec::with_capacity(1024)),
            base_max_qual: std::array::from_fn(|_| Vec::with_capacity(1024)),
            ref_blocks: Vec::new(),
            mapq: 0,
            has_primary: false,
            supporting_sequence_ids: Vec::new(),
            mapped_segment_indices: Vec::new(),
        }
    }
}

impl CellRecombinationResults {
    pub(crate) fn add_seq(&mut self, sequence: &[u8], qualities: &[u8]) -> bool {
        self.add_seq_with_overlap(sequence, qualities, DEFAULT_SEQUENCE_OVERLAP)
    }

    fn add_seq_with_overlap(
        &mut self,
        sequence: &[u8],
        qualities: &[u8],
        min_overlap: usize,
    ) -> bool {
        if sequence.is_empty() || sequence.len() != qualities.len() {
            return false;
        }

        if self.is_empty() {
            self.add_at(sequence, qualities, 0);
            return true;
        }

        let (consensus, _) = self.finalize();
        let min_overlap = min_overlap.min(consensus.len()).min(sequence.len()).max(1);
        let Some(merge) = merge_sequences_any_orientation(&consensus, sequence, min_overlap) else {
            return false;
        };

        let (sequence_storage, qualities_storage) = if merge.right_reverse {
            (
                reverse_complement(sequence),
                qualities.iter().rev().copied().collect::<Vec<_>>(),
            )
        } else {
            (sequence.to_vec(), qualities.to_vec())
        };

        self.extend_left_if_needed(merge.right_offset);
        self.add_at(
            &sequence_storage,
            &qualities_storage,
            merge.right_offset.max(0) as usize,
        );
        true
    }

    pub(crate) fn add_mapping(&mut self, mapping: &InitialMapperResults) {
        self.mapq = self.mapq.max(mapping.mapq);
        self.has_primary |= !mapping.is_supplementary;
        merge_reference_blocks(&mut self.ref_blocks, &mapping.ref_blocks);
        for segment in &mapping.segments {
            if !self.mapped_segment_indices.contains(&segment.segment_index) {
                self.mapped_segment_indices.push(segment.segment_index);
            }
        }
    }

    pub(crate) fn finalize(&self) -> (Vec<u8>, Vec<u8>) {
        let len = self.base_counts.iter().map(Vec::len).max().unwrap_or(0);
        let mut sequence = Vec::with_capacity(len);
        let mut qualities = Vec::with_capacity(len);

        for pos in 0..len {
            let mut best: Option<usize> = None;
            for base in 0..4 {
                let count = self.base_counts[base].get(pos).copied().unwrap_or(0);
                if count == 0 {
                    continue;
                }
                let quality = self.base_max_qual[base].get(pos).copied().unwrap_or(0);
                if best.is_none_or(|old| {
                    let old_count = self.base_counts[old].get(pos).copied().unwrap_or(0);
                    let old_quality = self.base_max_qual[old].get(pos).copied().unwrap_or(0);
                    count > old_count || (count == old_count && quality > old_quality)
                }) {
                    best = Some(base);
                }
            }

            if let Some(base) = best {
                sequence.push(b"ACGT"[base]);
                qualities.push(self.base_max_qual[base].get(pos).copied().unwrap_or(0));
            } else {
                sequence.push(b'N');
                qualities.push(0);
            }
        }

        (sequence, qualities)
    }

    fn is_empty(&self) -> bool {
        self.base_counts.iter().all(Vec::is_empty)
    }

    fn add_at(&mut self, sequence: &[u8], qualities: &[u8], offset: usize) {
        let needed = offset + sequence.len();
        self.ensure_len(needed);
        for (idx, (&base, &qual)) in sequence.iter().zip(qualities).enumerate() {
            let Some(base_idx) = dna_base_index(base) else {
                continue;
            };
            let pos = offset + idx;
            self.base_counts[base_idx][pos] = self.base_counts[base_idx][pos].saturating_add(1);
            self.base_max_qual[base_idx][pos] = self.base_max_qual[base_idx][pos].max(qual);
        }
    }

    fn ensure_len(&mut self, len: usize) {
        for base in 0..4 {
            self.base_counts[base].resize(len, 0);
            self.base_max_qual[base].resize(len, 0);
        }
    }

    fn extend_left_if_needed(&mut self, offset: isize) {
        if offset >= 0 {
            return;
        }
        let prefix = (-offset) as usize;
        for base in 0..4 {
            let mut counts = vec![0u16; prefix];
            counts.append(&mut self.base_counts[base]);
            self.base_counts[base] = counts;

            let mut qualities = vec![0u8; prefix];
            qualities.append(&mut self.base_max_qual[base]);
            self.base_max_qual[base] = qualities;
        }
    }
}

pub(crate) fn absorb_provisional(cell: &mut CellEvidence) -> (usize, usize) {
    let provisional = std::mem::take(&mut cell.provisional);
    let mut processed = 0usize;
    let mut merged = 0usize;

    for (parts, mappings) in provisional {
        // A provisional BAM feature sequence is the evidence unit. Multiple
        // mapper records (primary/supplementary or repeated segment hits) only
        // annotate its placement; they must never add the same bases twice.
        for read2 in [false, true] {
            let Some((sequence, qualities)) = parts.sequence(read2) else {
                continue;
            };

            let mut chains = Vec::<Chain>::new();
            for mapping in mappings.iter().filter(|mapping| mapping.read2 == read2) {
                for segment in &mapping.segments {
                    if !chains.contains(&segment.chain) {
                        chains.push(segment.chain);
                    }
                }
            }

            for chain in chains {
                processed += 1;
                let contigs = cell.final_results.entry(chain).or_default();
                let before = contigs.len();

                let mut target = None;
                for (idx, contig) in contigs.iter_mut().enumerate() {
                    if contig.add_seq(sequence, qualities) {
                        target = Some(idx);
                        break;
                    }
                }

                let idx = if let Some(idx) = target {
                    idx
                } else {
                    let mut contig = CellRecombinationResults::default();
                    debug_assert!(contig.add_seq(sequence, qualities));
                    contigs.push(contig);
                    contigs.len() - 1
                };
                if contigs.len() == before {
                    merged += 1;
                }

                if !contigs[idx]
                    .supporting_sequence_ids
                    .contains(&parts.sequence_id)
                {
                    contigs[idx].supporting_sequence_ids.push(parts.sequence_id);
                }
                for mapping in mappings.iter().filter(|mapping| {
                    mapping.read2 == read2 && mapping.segments.iter().any(|segment| segment.chain == chain)
                }) {
                    contigs[idx].add_mapping(mapping);
                }
            }
        }
    }

    (processed, merged)
}

fn dna_base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn merge_reference_blocks(left: &mut Vec<(u32, u32)>, right: &[(u32, u32)]) {
    left.extend_from_slice(right);
    left.sort_unstable();
    let mut merged = Vec::<(u32, u32)>::with_capacity(left.len());
    for (start, end) in left.drain(..) {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    *left = merged;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consensus_prefers_count_before_quality() {
        let mut result = CellRecombinationResults::default();
        for _ in 0..20 {
            assert!(result.add_seq(b"A", &[35]));
        }
        assert!(result.add_seq(b"G", &[40]));
        let (seq, qual) = result.finalize();
        assert_eq!(seq, b"A");
        assert_eq!(qual, vec![35]);
    }

    #[test]
    fn consensus_uses_quality_to_break_equal_count() {
        let mut result = CellRecombinationResults::default();
        assert!(result.add_seq(b"A", &[25]));
        assert!(result.add_seq(b"G", &[38]));
        let (seq, qual) = result.finalize();
        assert_eq!(seq, b"G");
        assert_eq!(qual, vec![38]);
    }
}
