//! Optional per-UMI evidence. Ordinary count matrices use `()`.
pub trait UmiEvidence: Default + Clone + Send + Sync {
    fn merge(&mut self, other: &Self);
}

impl UmiEvidence for () {
    fn merge(&mut self, _other: &Self) {}
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
    fn saturates_and_detects_conflicts() {
        let mut evidence = SnpUmiEvidence::default();
        for _ in 0..300 { evidence.add_ref(17); }
        evidence.add_alt(17);
        assert_eq!(evidence.snps[0].ref_reads, 255);
        assert_eq!(evidence.snps[0].alt_reads, 1);
        assert_eq!(evidence.conflicts().count(), 1);
    }
}
