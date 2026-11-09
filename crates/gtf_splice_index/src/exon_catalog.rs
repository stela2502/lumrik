//! Gene-wide exon units and compact, variable-length exon evidence.
//! Each gene defines its own ordered catalogue. Masks are meaningful only
//! within that gene's catalogue.
use crate::model::transcript::Transcript;
use crate::types::{RefBlock, SplicedRead};

pub use scdata::ExonMask;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneExonCatalog {
    /// Non-overlapping genomic intervals, ordered by coordinate.
    pub units: Vec<RefBlock>,
}

impl GeneExonCatalog {
    pub fn from_transcripts<'a>(transcripts: impl IntoIterator<Item = &'a Transcript>) -> Self {
        let mut intervals = Vec::new();
        let mut boundaries = Vec::new();
        for tx in transcripts {
            for &exon in tx.exons() {
                intervals.push(exon);
                boundaries.push(exon.start);
                boundaries.push(exon.end);
            }
        }
        boundaries.sort_unstable();
        boundaries.dedup();
        let units = boundaries.windows(2).filter_map(|pair| {
            let unit = RefBlock { start: pair[0], end: pair[1] };
            intervals.iter().any(|exon| exon.contains(unit)).then_some(unit)
        }).collect();
        Self { units }
    }
    /// Interleaved evidence: even bits are exon units; odd bits are the
    /// genomic gaps between consecutive units (potential intronic sequence).
    /// WARNING: Within a (cell, UMI, gene), masks are OR-merged. Unobserved
    /// holes are ignored; no connectivity or contradiction check is performed.
    pub fn observed_mask(&self, read: &SplicedRead) -> ExonMask {
        self.mask_for_blocks(&read.blocks)
    }
    pub fn transcript_mask(&self, tx: &Transcript) -> ExonMask {
        self.mask_for_blocks(tx.exons())
    }
    pub fn mask_for_blocks(&self, blocks: &[RefBlock]) -> ExonMask {
        let mut mask = ExonMask::empty(self.units.len().saturating_mul(2).saturating_sub(1));
        for (i, unit) in self.units.iter().enumerate() {
            if blocks.iter().any(|block| block.overlaps(*unit)) { mask.set(2 * i); }
            if let Some(next) = self.units.get(i + 1) {
                if unit.end < next.start {
                    let gap = RefBlock { start: unit.end, end: next.start };
                    if blocks.iter().any(|block| block.overlaps(gap)) { mask.set(2 * i + 1); }
                }
            }
        }
        mask
    }
    pub fn compatible(transcript_mask: &ExonMask, observed_mask: &ExonMask) -> bool {
        transcript_mask.contains(observed_mask)
    }
}

/// Classification of a complete (cell, UMI, gene) evidence mask.
/// Novel patterns are evidence masks, not inferred splice-junction models.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UmiTranscriptMatch {
    Intronic,
    Reference(usize),
    Novel(ExonMask),
}

impl GeneExonCatalog {
    /// Resolve one OR-merged UMI. Intronic evidence takes precedence.
    /// For mature molecules, compare reference exon usage only between the
    /// first and last observed exons; missing transcript ends are tolerated.
    /// Ties resolve deterministically to the lowest transcript ID.
    /// WARNING: holes within one UMI are ignored during collection; this
    /// function cannot reconstruct unobserved exon connectivity.
    pub fn resolve_umi<'a>(
        &self,
        observed: &ExonMask,
        transcripts: impl IntoIterator<Item = &'a Transcript>,
    ) -> UmiTranscriptMatch {
        let n = self.units.len();
        if (0..n.saturating_sub(1)).any(|i| mask_bit(observed, 2 * i + 1)) {
            return UmiTranscriptMatch::Intronic;
        }
        let first = (0..n).find(|&i| mask_bit(observed, 2 * i));
        let last = (0..n).rfind(|&i| mask_bit(observed, 2 * i));
        let (Some(first), Some(last)) = (first, last) else {
            return UmiTranscriptMatch::Novel(observed.clone());
        };
        // Missing terminal coverage is not a structural difference. Choose the
        // closest annotated model over the observed span, rather than requiring
        // an exact match or treating unsequenced 5' exons as deletions.
        let mut best: Option<(usize, usize)> = None;
        for tx in transcripts {
            let reference = self.transcript_mask(tx);
            let differences = (first..=last)
                .filter(|&i| mask_bit(observed, 2 * i) != mask_bit(&reference, 2 * i))
                .count();
            if best.is_none_or(|(score, id)| (differences, tx.id) < (score, id)) {
                best = Some((differences, tx.id));
            }
        }
        match best {
            Some((0, id)) => UmiTranscriptMatch::Reference(id),
            Some(_) | None => UmiTranscriptMatch::Novel(observed.clone()),
        }
    }
}

fn mask_bit(mask: &ExonMask, bit: usize) -> bool {
    match mask {
        ExonMask::Small(value) => bit < 64 && (value & (1u64 << bit)) != 0,
        ExonMask::Large(words) => words.get(bit / 64)
            .is_some_and(|value| value & (1u64 << (bit % 64)) != 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn classify_intronic_and_reference_and_novel() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(30, 40), RefBlock::new(50, 60),
        ] };
        let mut tx = Transcript::new(7, 0, "T7", 0, crate::Strand::Plus);
        tx.add_exon(RefBlock::new(10, 20));
        tx.add_exon(RefBlock::new(30, 40));
        tx.add_exon(RefBlock::new(50, 60));
        assert_eq!(catalog.resolve_umi(&ExonMask::Small(0b00100), [&tx]), UmiTranscriptMatch::Reference(7));
        assert_eq!(catalog.resolve_umi(&ExonMask::Small(0b00010), [&tx]), UmiTranscriptMatch::Intronic);
        assert_eq!(catalog.resolve_umi(&ExonMask::Small(0b10001), [&tx]), UmiTranscriptMatch::Novel(ExonMask::Small(0b10001)));
    }
    #[test]
    fn missing_terminal_exons_are_not_novel() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(30, 40),
            RefBlock::new(50, 60), RefBlock::new(70, 80),
        ] };
        let mut tx = Transcript::new(7, 0, "T7", 0, crate::Strand::Minus);
        for unit in &catalog.units { tx.add_exon(*unit); }
        let observed = catalog.mask_for_blocks(&catalog.units[1..]);
        assert_eq!(catalog.resolve_umi(&observed, [&tx]), UmiTranscriptMatch::Reference(7));
    }
    #[test]
    fn internal_missing_exon_is_novel_without_exact_reference() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(30, 40),
            RefBlock::new(50, 60), RefBlock::new(70, 80),
        ] };
        let mut tx = Transcript::new(7, 0, "T7", 0, crate::Strand::Plus);
        for unit in &catalog.units { tx.add_exon(*unit); }
        let observed = catalog.mask_for_blocks(&[
            catalog.units[0], catalog.units[2], catalog.units[3]
        ]);
        assert_eq!(catalog.resolve_umi(&observed, [&tx]), UmiTranscriptMatch::Novel(observed));
    }
    #[test]
    fn combine_disjoint_observations() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(20, 30),
            RefBlock::new(50, 60), RefBlock::new(80, 90),
        ] };
        let mut a = catalog.mask_for_blocks(&[RefBlock::new(10, 19)]);
        let b = catalog.mask_for_blocks(&[RefBlock::new(80, 90)]);
        a.merge(&b);
        assert_eq!(a, ExonMask::Small(0b1000001));
        assert!(GeneExonCatalog::compatible(&ExonMask::Small(0b1010101), &a));
        assert!(!GeneExonCatalog::compatible(&ExonMask::Small(0b0010101), &a));
    }
    #[test]
    fn interleaved_intron_evidence_and_holes() {
        let catalog = GeneExonCatalog { units: vec![
            RefBlock::new(10, 20), RefBlock::new(30, 40), RefBlock::new(50, 60),
        ] };
        let exon_only = catalog.mask_for_blocks(&[RefBlock::new(10, 20), RefBlock::new(50, 60)]);
        assert_eq!(exon_only, ExonMask::Small(0b10001));
        let retained = catalog.mask_for_blocks(&[RefBlock::new(20, 30)]);
        assert_eq!(retained, ExonMask::Small(0b00010));
        let mut merged = exon_only;
        merged.merge(&retained);
        assert_eq!(merged, ExonMask::Small(0b10011));
    }
    #[test]
    fn more_than_64_units_preserved() {
        let catalog = GeneExonCatalog { units: (0..100).map(|i| RefBlock::new(i*2, i*2+1)).collect() };
        let mut a = catalog.mask_for_blocks(&[RefBlock::new(0, 1)]);
        a.merge(&catalog.mask_for_blocks(&[RefBlock::new(198, 199)]));
        assert_eq!(a.count_ones(), 2);
        assert!(a.contains(&catalog.mask_for_blocks(&[RefBlock::new(198, 199)])));
    }
    #[test]
    fn small_and_large_merge_without_losing_bits() {
        let mut a = ExonMask::Small(1);
        let mut b = ExonMask::empty(96);
        b.set(95);
        a.merge(&b);
        assert_eq!(a.count_ones(), 2);
        assert!(a.contains(&ExonMask::Small(1)));
        assert!(a.contains(&b));
    }
}
