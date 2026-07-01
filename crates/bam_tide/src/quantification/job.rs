//job.rs
use anyhow::{Context, Result};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{HeaderView, Record};

use crate::core::ref_block::record_to_blocks;
use gtf_splice_index::types::RefBlock;
use gtf_splice_index::{SplicedRead, Strand};

use int_to_str::int_to_str::IntToStr;

use mapping_info::MappingInfo;
use snp_index::{AlignedRead, Genome, RefineOptions, SnpIndex};

use read_tag_table::ReadTagTable;

#[derive(Clone)]
pub struct Job {
    pub cell: u64,
    pub umi: u64,
    pub spliced: SplicedRead,

    /// Refined/cleaned read used for SNP matching.
    ///
    /// Currently created when genome/SNP support is active.
    pub aligned: Option<AlignedRead>,
}

pub struct JobBuilder<'a> {
    header: &'a HeaderView,
    chr_map: &'a std::collections::HashMap<String, usize>,
    genome: Option<&'a Genome>,
    snp: Option<&'a SnpIndex>,
    read_tag_table: Option<&'a ReadTagTable>,

    cell_tag: [u8; 2],
    umi_tag: [u8; 2],

    min_mapq: u8,
    read1_only: bool,
    refine_against_genome: bool,
}

impl<'a> JobBuilder<'a> {
    pub fn new(
        header: &'a HeaderView,
        chr_map: &'a std::collections::HashMap<String, usize>,
        cell_tag: [u8; 2],
        umi_tag: [u8; 2],
    ) -> Self {
        Self {
            header,
            chr_map,
            genome: None,
            snp: None,
            read_tag_table: None,
            cell_tag,
            umi_tag,
            min_mapq: 0,
            read1_only: false,
            refine_against_genome: false,
        }
    }

    pub fn with_genome(mut self, genome: Option<&'a Genome>, refine_against_genome: bool) -> Self {
        self.genome = genome;
        self.refine_against_genome = refine_against_genome;
        self
    }

    pub fn with_snp_index(mut self, snp: Option<&'a SnpIndex>) -> Self {
        self.snp = snp;
        self
    }

    pub fn with_read_tag_table(mut self, table: Option<&'a ReadTagTable>) -> Self {
        self.read_tag_table = table;
        self
    }

    pub fn with_min_mapq(mut self, min_mapq: u8) -> Self {
        self.min_mapq = min_mapq;
        self
    }

    pub fn read1_only(mut self, read1_only: bool) -> Self {
        self.read1_only = read1_only;
        self
    }



    pub fn build(&self, rec: &Record, report: &mut MappingInfo) -> Result<Option<Job>> {
        if rec.is_unmapped() {
            report.report("unmapped");
            return Ok(None);
        }

        if rec.mapq() < self.min_mapq {
            report.report("mapq failed");
            return Ok(None);
        }

        if rec.is_secondary() || rec.is_supplementary() {
            report.report("secondary or supplementary");
            return Ok(None);
        }

        if self.read1_only && !rec.is_first_in_template() {
            report.report("read!=1");
            return Ok(None);
        }

        let (cb_raw, ub) = if let Some(read_tag_table) = self.read_tag_table {
            let read_id = std::str::from_utf8(rec.qname())
                .context("Invalid BAM read name / qname")?;

            match read_tag_table.cell_umi_for_read(read_id) {
                Some((cell, umi)) => {
                    report.report("cell_umi from read-tag table");
                    (cell.to_string(), umi.to_string())
                }
                None => {
                    report.report("read not in read-tag table");
                    return Ok(None);
                }
            }
        } else {
            let cb_raw = match Self::aux_tag_str(rec, self.cell_tag ) {
                Some(v) => v,
                None => {
                    report.report(
                        &format!("no {} tag", String::from_utf8_lossy(&self.cell_tag))
                    );
                    return Ok(None);
                }
            };

            let ub = match Self::aux_tag_str(rec, self.umi_tag ) {
                Some(v) => v,
                None => {
                    report.report(
                        &format!("no {} tag",  
                        String::from_utf8_lossy(&self.umi_tag))
                    );
                    return Ok(None);
                }
            };

            report.report("cell_umi from BAM tags");
            (cb_raw.to_string(), ub.to_string())
        };

        let cb = Self::normalize_10x_barcode(&cb_raw);

        let cell = match Self::dna_to_u64(cb) {
            Some(v) => v,
            None => {
                report.report("invalid CB");
                return Ok(None);
            }
        };

        let umi = match Self::dna_to_u64(&ub) {
            Some(v) => v,
            None => {
                report.report("invalid UB");
                return Ok(None);
            }
        };

        let tid = rec.tid();

        if tid < 0 {
            report.report("tid below 0");
            return Ok(None);
        }

        let chr_name = std::str::from_utf8(self.header.tid2name(tid as u32))
            .context("Invalid chromosome name in BAM header")?;

        let chr_id = match self.fuzzy_chr_id(chr_name) {
            Some(id) => id,
            None => {
                report.report(&format!("contig {chr_name} not in index - checked with and without chr"));
                return Ok(None);
            }
        };

        let spliced = match Self::record_to_spliced_read(rec, chr_id) {
            Some(v) => v,
            None => {
                report.report("no ref blocks");
                return Ok(None);
            }
        };

        let aligned = self.build_aligned_read(rec, chr_name);

        Ok(Some(Job {
            cell,
            umi,
            spliced,
            aligned,
        }))
    }

    fn fuzzy_chr_id(&self, chr_name: &str) -> Option<usize> {
        if let Some(&id) = self.chr_map.get(chr_name) {
            return Some(id);
        }

        if let Some(stripped) = chr_name.strip_prefix("chr") {
            if let Some(&id) = self.chr_map.get(stripped) {
                return Some(id);
            }

            if stripped == "M" {
                if let Some(&id) = self.chr_map.get("MT") {
                    return Some(id);
                }
            }
        } else {
            let with_chr = format!("chr{chr_name}");

            if let Some(&id) = self.chr_map.get(&with_chr) {
                return Some(id);
            }

            match chr_name {
                "MT" => {
                    if let Some(&id) = self.chr_map.get("chrM") {
                        return Some(id);
                    }
                    if let Some(&id) = self.chr_map.get("M") {
                        return Some(id);
                    }
                }
                "M" => {
                    if let Some(&id) = self.chr_map.get("chrM") {
                        return Some(id);
                    }
                    if let Some(&id) = self.chr_map.get("MT") {
                        return Some(id);
                    }
                }
                _ => {}
            }
        }

        None
    }

    fn build_aligned_read(&self, rec: &Record, chr_name: &str) -> Option<AlignedRead> {
        if self.genome.is_none() && self.snp.is_none() {
            return None;
        }

        let chr_id = if let Some(snp) = self.snp {
            snp.chr_id(chr_name)?
        } else {
            *self.chr_map.get(chr_name)?
        };

        let mut read = AlignedRead::from_record(rec, chr_id);

        if let Some(genome) = self.genome
            && self.refine_against_genome
        {
            read.refine_against_genome(genome, RefineOptions::default());
        }

        Some(read)
    }

    fn aux_tag_str(rec: &Record, tag: [u8; 2]) -> Option<&str> {
        match rec.aux(&tag).ok()? {
            Aux::String(s) => Some(s),
            _ => None,
        }
    }

    fn normalize_10x_barcode(cb: &str) -> &str {
        match cb.split_once('-') {
            Some((core, _)) => core,
            None => cb,
        }
    }

    /// a temporary bugfid for the https://github.com/imallona/rock_roi_method
    /// pile of crap - sorry
    fn normalize_rock_roi(raw: &str) -> String {
        raw.chars()
            .filter(|c| *c != '_' && *c != '-')
            .collect::<String>()
    }

    fn dna_to_u64(seq: &str) -> Option<u64> {

        let seq = Self::normalize_rock_roi( seq );

        if !seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            return None;
        }

        let tool = IntToStr::new(seq.as_bytes());
        Some(tool.into_u64())
    }

    fn record_to_spliced_read(rec: &Record, chr_id: usize) -> Option<SplicedRead> {
        let blocks: Vec<RefBlock> = record_to_blocks(rec);

        if blocks.is_empty() {
            return None;
        }

        let strand = if rec.is_reverse() {
            Strand::Minus
        } else {
            Strand::Plus
        };

        let mut spliced = SplicedRead::new(chr_id, strand, blocks);
        spliced.finalize();

        Some(spliced)
    }
}
