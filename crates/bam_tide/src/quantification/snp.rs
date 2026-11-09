//snp.rs
use anyhow::{Context, Result};
use run_stats::RunStats;
use scdata::SnpUmiEvidence;
use snp_index::{AlignedRead, SnpIndex, VcfReadOptions};

pub struct SnpSideChannel {
    pub index: SnpIndex,
    min_anchor: u8,
    flank: u8,
}

impl SnpSideChannel {
    pub fn from_vcf_path(
        path: impl AsRef<std::path::Path>,
        chr_names: Vec<String>,
        chr_lengths: Vec<u32>,
        min_anchor: u8,
    flank: u8,
    ) -> Result<Self> {
        let path = path.as_ref();
        let opts = VcfReadOptions::default();

        let index = SnpIndex::from_vcf_path(path, chr_names, chr_lengths, 10_000, &opts)
            .with_context(|| format!("reading SNP VCF {}", path.display()))?;

        eprintln!("{index}");

        Ok(Self { index, min_anchor, flank })
    }

    pub fn add_hits(
        &self,
        aligned: Option<&AlignedRead>,
        report: &mut RunStats,
    ) -> SnpUmiEvidence {
        let mut evidence = SnpUmiEvidence::default();
        let Some(read) = aligned else {
            report.report("snp requested but no aligned read");
            return evidence;
        };

        //panic!("We use the min SNP quality of {} here!\nSNP index info: {}", self.min_anchor, self.index );
        let (ref_ids, alt_ids, _other_ids) = self
            .index
            .get_ref_alt_other_ids_for_read_with_flank(read, self.min_anchor, self.flank);

        if ref_ids.is_empty() && alt_ids.is_empty() {
            report.report("no snp hit");
            return evidence;
        }

        for snp_id in ref_ids {
            evidence.add_ref(snp_id);
        }

        for snp_id in alt_ids {
            evidence.add_alt(snp_id);
        }

        report.report("snp hit");
        evidence
    }

    pub fn min_anchor_for_diagnostic(&self) -> u8 { self.min_anchor }
    pub fn flank_for_diagnostic(&self) -> u8 { self.flank }

    pub fn is_empty(&self) -> bool {
        false
    }
}
