//snp.rs
use anyhow::{Context, Result};
use mapping_info::MappingInfo;
use scdata::Scdata;
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
        sc_ref: &mut Scdata,
        sc_alt: &mut Scdata,
        report: &mut MappingInfo,
    ) {
        let Some(read) = aligned else {
            report.report("snp requested but no aligned read");
            return;
        };

        //panic!("We use the min SNP quality of {} here!\nSNP index info: {}", self.min_anchor, self.index );
        let (ref_ids, alt_ids, _other_ids) = self
            .index
            .get_ref_alt_other_ids_for_read(read, self.min_anchor);

        if ref_ids.is_empty() && alt_ids.is_empty() {
            report.report("no snp hit");
            return;
        }

        for snp_id in ref_ids {
            sc_ref.try_insert(&cell, GeneUmiHash(snp_id, umi), 1.0, report);
        }

        for snp_id in alt_ids {
            sc_alt.try_insert(&cell, GeneUmiHash(snp_id, umi), 1.0, report);
        }

        report.report("snp hit");
    }

    pub fn is_empty(&self) -> bool {
        false
    }
}
