//model.rs
use std::fmt::Write;

use crate::encoder::cli::TruthFeatureMode;
use gtf_splice_index::MatchClass;

#[derive(Debug, Clone)]
pub struct EncoderConfig {
    pub seq_len: usize,
    pub reads_per_cell: usize,

    pub antisense_fraction: f64,
    pub unspliced_fraction: f64,
    pub snp_fraction: f64,

    pub body_error_rate: f64,
    pub end_error_rate: f64,
    pub bad_end_bases: usize,

    pub seed: u64,
    pub cell_barcodes: Vec<String>,

    pub truth_path: Option<std::path::PathBuf>,
    pub min_cell_counts: usize,
    pub genes_per_cell: usize,
    pub truth_feature_mode: TruthFeatureMode,
}

impl Default for EncoderConfig {
    fn default() -> Self {
        Self {
            seq_len: 30,
            reads_per_cell: 70,

            antisense_fraction: 10.0 / 70.0,
            unspliced_fraction: 10.0 / 70.0,
            snp_fraction: 0.20,

            body_error_rate: 0.005,
            end_error_rate: 0.15,
            bad_end_bases: 4,

            seed: 1,
            cell_barcodes: vec!["CELL000001-1".to_string()],
            min_cell_counts: 1,
            genes_per_cell: 1,
            truth_path: None,

            truth_feature_mode: TruthFeatureMode::Gene,
        }
    }
}

impl EncoderConfig {
    pub fn validate(&self) -> Result<(), String> {
        if self.seq_len == 0 {
            return Err("seq_len must be > 0".into());
        }

        if self.reads_per_cell == 0 {
            return Err("reads_per_cell must be > 0".into());
        }

        if self.antisense_fraction + self.unspliced_fraction > 1.0 {
            return Err("antisense_fraction + unspliced_fraction must be <= 1".into());
        }

        if self.cell_barcodes.is_empty() {
            return Err("at least one cell barcode is required".into());
        }

        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct SamRead {
    pub qname: String,
    pub flag: u16,
    pub chrom: String,
    pub pos_1based: u64,
    pub mapq: u8,
    pub cigar: String,
    pub seq: String,
    pub qual: String,

    pub cell_barcode: String,
    pub umi: String,

    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub transcript_id: Option<String>,

    pub nm: usize,
    pub expected: MatchClass,
}

impl SamRead {
    pub fn to_sam_line(&self) -> String {
        let raw_cell = self.cell_barcode.trim_end_matches("-1");

        let mut line = format!(
            "{qname}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\
             \tCB:Z:{cb}\tCR:Z:{cr}\tCY:Z:{cy}\
             \tUB:Z:{umi}\tUR:Z:{umi}\tUY:Z:{uy}\
             \tNH:i:1\tNM:i:{nm}\tZE:Z:{expected}",
            qname = self.qname,
            flag = self.flag,
            chrom = self.chrom,
            pos = self.pos_1based,
            mapq = self.mapq,
            cigar = self.cigar,
            seq = self.seq,
            qual = self.qual,
            cb = self.cell_barcode,
            cr = raw_cell,
            cy = "FFFFFFFFFFFFFFFF",
            umi = self.umi,
            uy = "FFFFFFFFFFFF",
            nm = self.nm,
            expected = self.expected,
        );

        if let Some(gene_id) = &self.gene_id {
            write!(line, "\tGX:Z:{gene_id}").unwrap();
        }

        if let Some(gene_name) = &self.gene_name {
            write!(line, "\tGN:Z:{gene_name}").unwrap();
        }

        if let Some(transcript_id) = &self.transcript_id {
            write!(line, "\tTX:Z:{transcript_id}").unwrap();
        }

        line
    }
}
