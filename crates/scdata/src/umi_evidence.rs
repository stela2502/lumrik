//! Optional per-UMI evidence. Ordinary count matrices use `()`.
pub trait UmiEvidence: Default + Clone + Send + Sync {
    fn merge(&mut self, other: &Self);
}

impl UmiEvidence for () {
    fn merge(&mut self, _other: &Self) {}
}

/// Gene-wide observed exon units for a single (cell, UMI, gene).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ExonMask {
    Small(u64),
    Large(Vec<u64>),
}

impl Default for ExonMask {
    fn default() -> Self { Self::Small(0) }
}

impl ExonMask {
    pub fn empty(units: usize) -> Self {
        if units <= 64 { Self::Small(0) } else { Self::Large(vec![0; units.div_ceil(64)]) }
    }

    pub fn set(&mut self, bit: usize) {
        match self {
            Self::Small(value) => {
                if bit < 64 { *value |= 1u64 << bit; }
                else {
                    let mut words = vec![0; bit / 64 + 1];
                    words[0] = *value;
                    words[bit / 64] |= 1u64 << (bit % 64);
                    *self = Self::Large(words);
                }
            }
            Self::Large(words) => {
                words.resize(words.len().max(bit / 64 + 1), 0);
                words[bit / 64] |= 1u64 << (bit % 64);
            }
        }
    }

    fn word(&self, i: usize) -> u64 {
        match self {
            Self::Small(value) => if i == 0 { *value } else { 0 },
            Self::Large(words) => words.get(i).copied().unwrap_or(0),
        }
    }

    fn words(&self) -> usize {
        match self { Self::Small(_) => 1, Self::Large(words) => words.len() }
    }

    pub fn merge(&mut self, other: &Self) {
        let count = self.words().max(other.words());
        if count == 1 {
            *self = Self::Small(self.word(0) | other.word(0));
        } else {
            *self = Self::Large((0..count).map(|i| self.word(i) | other.word(i)).collect());
        }
    }

    pub fn contains(&self, observed: &Self) -> bool {
        (0..observed.words()).all(|i| self.word(i) & observed.word(i) == observed.word(i))
    }

    pub fn count_ones(&self) -> u32 {
        (0..self.words()).map(|i| self.word(i).count_ones()).sum()
    }
}

impl UmiEvidence for ExonMask {
    fn merge(&mut self, other: &Self) { ExonMask::merge(self, other); }
}

/// Splice and variant observations from the same (cell, UMI, gene).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TranscriptUmiEvidence {
    pub splice: ExonMask,
    pub snps: SnpUmiEvidence,
}

impl UmiEvidence for TranscriptUmiEvidence {
    fn merge(&mut self, other: &Self) {
        self.splice.merge(&other.splice);
        self.snps.merge(&other.snps);
    }
}

/// Gene-wide observed-exon bitmap carried by one feature/UMI observation.
/// Repeated observations combine without counting an exon twice.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ExonMaskEvidence(pub u64);

impl UmiEvidence for ExonMaskEvidence {
    #[inline]
    fn merge(&mut self, other: &Self) {
        self.0 |= other.0;
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SnpSupport {
    pub snp_id: u64,
    pub ref_reads: u8,
    pub alt_reads: u8,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SnpUmiEvidence {
    pub snps: Vec<SnpSupport>,
}

impl UmiEvidence for SnpUmiEvidence {
    fn merge(&mut self, other: &Self) {
        for incoming in &other.snps {
            if let Some(existing) = self.snps.iter_mut().find(|s| s.snp_id == incoming.snp_id) {
                existing.ref_reads = existing.ref_reads.saturating_add(incoming.ref_reads);
                existing.alt_reads = existing.alt_reads.saturating_add(incoming.alt_reads);
            } else {
                self.snps.push(incoming.clone());
            }
        }
    }
}

impl SnpUmiEvidence {
    pub fn add_ref(&mut self, snp_id: u64) {
        self.merge(&Self { snps: vec![SnpSupport { snp_id, ref_reads: 1, alt_reads: 0 }] });
    }
    pub fn add_alt(&mut self, snp_id: u64) {
        self.merge(&Self { snps: vec![SnpSupport { snp_id, ref_reads: 0, alt_reads: 1 }] });
    }
    pub fn conflicts(&self) -> impl Iterator<Item = &SnpSupport> {
        self.snps.iter().filter(|s| s.ref_reads > 0 && s.alt_reads > 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exon_masks_merge_by_or() {
        let mut evidence = ExonMaskEvidence(0b0011);
        evidence.merge(&ExonMaskEvidence(0b0110));
        assert_eq!(evidence.0, 0b0111);
        evidence.merge(&ExonMaskEvidence(0b0110));
        assert_eq!(evidence.0, 0b0111);
    }
    #[test]
    fn shared_exon_mask_merges_across_32_and_64_bits() {
        let mut a = ExonMask::Small(0b11);
        let mut b = ExonMask::empty(100);
        b.set(70);
        a.merge(&b);
        assert_eq!(a.count_ones(), 3);
        assert!(a.contains(&b));
        a.merge(&b);
        assert_eq!(a.count_ones(), 3);
    }
    #[test]
    fn transcript_umi_evidence_merges_splice_and_snp() {
        let mut a = TranscriptUmiEvidence::default();
        a.splice.set(0);
        a.snps.add_ref(17);
        let mut b = TranscriptUmiEvidence::default();
        b.splice.set(70);
        b.snps.add_alt(17);
        a.merge(&b);
        assert_eq!(a.splice.count_ones(), 2);
        assert_eq!(a.snps.snps.len(), 1);
        assert_eq!(a.snps.snps[0].ref_reads, 1);
        assert_eq!(a.snps.snps[0].alt_reads, 1);
    }
    #[test]
    fn saturates_and_detects_conflicts() {
        let mut evidence = SnpUmiEvidence::default();
        for _ in 0..300 { evidence.add_ref(17); }
        evidence.add_alt(17);
        assert_eq!(evidence.snps[0].ref_reads, 255);
        assert_eq!(evidence.snps[0].alt_reads, 1);
        assert_eq!(evidence.conflicts().count(), 1);
    }
}
