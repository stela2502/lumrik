use std::fmt;

#[derive(Debug, Default, Clone)]
pub struct NormalizeStats {
    pub total_records: u64,
    pub zero_cassette: u64,
    pub one_cassette: u64,
    pub multi_cassette: u64,
    pub emitted_molecules: u64,
    pub forward_molecules: u64,
    pub reverse_molecules: u64,
    pub too_short_after_adapter: u64,
    pub failed_poly_t: u64,
    pub fastq_reads_written: u64,
}

impl fmt::Display for NormalizeStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mean = if self.total_records == 0 {
            0.0
        } else {
            self.emitted_molecules as f64 / self.total_records as f64
        };

        writeln!(f, "bam-ont-normalizer stats")?;
        writeln!(f, "  total BAM records         : {}", self.total_records)?;
        writeln!(f, "  records with 0 cassettes  : {}", self.zero_cassette)?;
        writeln!(f, "  records with 1 cassette   : {}", self.one_cassette)?;
        writeln!(f, "  records with 2+ cassettes : {}", self.multi_cassette)?;
        writeln!(
            f,
            "  total molecules emitted   : {}",
            self.emitted_molecules
        )?;
        writeln!(
            f,
            "  forward molecules         : {}",
            self.forward_molecules
        )?;
        writeln!(
            f,
            "  reverse molecules         : {}",
            self.reverse_molecules
        )?;
        writeln!(
            f,
            "  too_short_after_adapter   : {}",
            self.too_short_after_adapter
        )?;
        writeln!(f, "  failed_polyT              : {}", self.failed_poly_t)?;
        writeln!(f, "  mean molecules per read   : {:.3}", mean)?;
        writeln!(
            f,
            "  FASTQ reads written       : {}",
            self.fastq_reads_written
        )?;

        Ok(())
    }
}
