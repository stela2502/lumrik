//snp.rs
use anyhow::{Context, Result};
use run_stats::RunStats;
use scdata::{QuantData, SnpUmiEvidence};
use scdata::cell_data::GeneUmiHash;
use snp_index::{AlignedRead, SnpIndex, VcfReadOptions};

pub struct SnpSideChannel {
    pub index: SnpIndex,
    min_anchor: u8,
}

impl SnpSideChannel {
    pub fn from_vcf_path(
        path: impl AsRef<std::path::Path>,
        chr_names: Vec<String>,
        chr_lengths: Vec<u32>,
        min_anchor: u8,
    ) -> Result<Self> {
        let path = path.as_ref();
        let opts = VcfReadOptions::default();

        let index = SnpIndex::from_vcf_path(path, chr_names, chr_lengths, 10_000, &opts)
            .with_context(|| format!("reading SNP VCF {}", path.display()))?;

        eprintln!("{index}");

        Ok(Self { index, min_anchor })
    }

    pub fn add_hits(
        &self,
        aligned: Option<&AlignedRead>,
        cell: u64,
        umi: u64,
        data: &mut QuantData,
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
            .get_ref_alt_other_ids_for_read(read, self.min_anchor);

        if ref_ids.is_empty() && alt_ids.is_empty() {
            report.report("no snp hit");
            return evidence;
        }

        for snp_id in ref_ids {
            evidence.add_ref(snp_id);
            data.try_insert(
                QuantData::SNP_REF,
                &cell,
                GeneUmiHash(snp_id, umi),
                1.0,
                report,
            );
        }

        for snp_id in alt_ids {
            evidence.add_alt(snp_id);
            data.try_insert(
                QuantData::SNP_ALT,
                &cell,
                GeneUmiHash(snp_id, umi),
                1.0,
                report,
            );
        }

        report.report("snp hit");
        evidence
    }

    pub fn is_empty(&self) -> bool {
        false
    }
}
