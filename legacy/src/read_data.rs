use rust_htslib::bam::record::{Record};
use rust_htslib::errors::Result;
use std::collections::HashMap;
use int_to_str::int_to_str::IntToStr;
use crate::bam_flag::BamFlag;

use core::fmt;

/// This class converts the start from bam (0 based) to Gtf (1 based)
#[derive(Debug, Clone)]
pub struct ReadData {
    pub cell_id: String,
    pub umi: u64,
    pub start: i32,
    pub flag: BamFlag,
    pub cigar: String,
    pub chromosome: String,
    //pub is_reverse: bool,
    pub sequence: Vec<u8>,
    pub qualities: Vec<u8>,
}


// Implementing Display trait for SecondSeq
impl fmt::Display for ReadData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}\t{}\t{}\t{}", 
            self.chromosome, self.start, self.cigar, self.flag, self.cell_id, 
            String::from_utf8_lossy(&self.sequence) )
    }
}

/// Returns the reverse complement of a DNA sequence.
/// Supports A, T, C, G, N (case-insensitive).
pub fn rev_compl(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N', // fallback for unknown bases
        })
        .collect()
}

impl ReadData {
    /// Extracts all required info from `bam_feature` and makes it thread-safe.
    pub fn new(
        cell_id: &str,
        umi: u64,
        start: i32,
        flag: u16,
        cigar: &str,
        chromosome: &str,
        //is_reverse: bool,
        sequence: Vec<u8>,
        qualities: Vec<u8>,
    ) -> Self {
        Self {
            cell_id: cell_id.to_string(),
            umi,
            start: start, // Convert BAM to GTF notation
            flag: BamFlag::new( flag ),
            cigar: cigar.to_string(),
            chromosome: chromosome.to_string(),
            //is_reverse,
            sequence,
            qualities,
        }
    }


    // Replace bam::record::tags::{StringType, TagValue}; with rust_htslib's tag handling

    /// The DataTuple stores the data from the Bam files:
    /// ( 0,      1,   2,     3,     4,   5,                 6,                7,               )
    /// (cell_id, umi, start, cigar, chr, is_reverse_strand, sequence as u8's, quality as u8's  )
    /// In this case, we are using ReadData in the return values.
    pub fn get_tag_value(record: &Record, tag: &[u8; 2]) -> Option<String> {
        if let Some(value) = record.aux(tag).ok() {
            match value {
                rust_htslib::bam::record::Aux::String(s) => Some(s.to_string()),
                _ => None,
            }
        } else {
            None
        }
    }

    /// Get the values from a bowtie2-mapped single cell experiment
    pub fn from_singlecell_bowtie2<'a>(
        bam_feature: &'a Record,
        chromosome_mappings: &'a HashMap<i32, String>,
        pseudo_umi: u64
    ) -> Result<(String, Self), &'a str> {
        // Get the read ID (qname)
        let id = bam_feature.qname();
        let id_str = std::str::from_utf8(id).map_err(|_| "invalid UTF-8 in qname")?;

        // Split by colon and extract the last part
        let parts: Vec<&str> = id_str.split(':').collect();
        let cell_seq = parts.last().ok_or("missing read ID component")?;

        if cell_seq.len() < 16 {
            return Err("sequence too short to contain the cell id");
        }

        // Check if the last part is a valid DNA string and split it
        if !cell_seq.chars().all(|c| matches!(c, 'A' | 'T' | 'C' | 'G' | 'N')) {
            //panic!("Problem: this is not a cell id: {}",  cell_seq);
            return Err("last part of ID is not a valid DNA sequence");
        }

        //panic!("We have the cell id {} and the umi {} from the seq {}", cell_barcode, umi, cell_seq);

        // Now call your existing constructor for the record
        
        Self::from_bulk_bowtie2(bam_feature, chromosome_mappings, cell_seq , pseudo_umi )
    }

    /// get_values_bulk extracts the values and creates a ReadData object.
    /// the CIGAR here actually contains the MD tag!
    pub fn from_bulk_bowtie2<'a>(
        bam_feature: &'a Record,
        chromosmome_mappings: &'a HashMap<i32, String>,
        cell_id: &str,
        umi_u64: u64
    ) -> Result<(String, Self), &'a str> {
        // Extract the chromosome (reference name)
        let chr = match chromosmome_mappings.get(&bam_feature.tid()) {
            Some(name) => name,
            None => return Err("missing_Chromosome"),
        };

        // Try to extract the UMI (UB), and report if missing
        //let umi_u64 = IntToStr::new(umi).into_u64();

        // Extract the start position (1-based for your GTF comparison)
        let start = bam_feature.pos(); // This is 0-based, needs adjustment for GTF

        let cigar = Self::get_tag_value(bam_feature, b"MD").unwrap_or_else(|| bam_feature.seq().len().to_string() ); 

        // Create a new ReadData object
        let res = Self::new(
            cell_id,
            umi_u64,
            (start + 1).try_into().unwrap(), // Adjust for GTF 1-based start
            bam_feature.flags(),
            &cigar,
            chr,
            bam_feature.seq().as_bytes(),
            bam_feature.qual().to_vec(),
        );

        let id = std::str::from_utf8(bam_feature.qname()).unwrap().to_string(); // Convert read name to string

        Ok((id, res))
    }

    /// get_values will extract the required fields (CellID, UMI, etc.) from a BAM record.
    pub fn from_single_cell<'a>(
        bam_feature: &'a Record,
        chromosmome_mappings: &'a HashMap<i32, String>,
        bam_cell_tag: &[u8; 2],
        bam_umi_tag: &[u8; 2]    ) -> Result<(String, Self), &'a str> {
        
        let cell_id = match Self::get_tag_value(bam_feature, bam_cell_tag) {
            Some(id) => id, // Convert to string, or use empty string
            None => return Err("missing_CellID"),
        };

        Self::from_bulk(bam_feature, chromosmome_mappings, bam_umi_tag, 1_u64, &cell_id)
    }

    /// get_values_bulk extracts the values and creates a ReadData object.
    pub fn from_bulk<'a>(
        bam_feature: &'a Record,
        chromosmome_mappings: &'a HashMap<i32, String>,
        bam_umi_tag: &[u8; 2],
        pseudo_umi: u64,
        cell_id: &str
    ) -> Result<(String, Self), &'a str> {
        // Extract the chromosome (reference name)
        let chr = match chromosmome_mappings.get(&bam_feature.tid()) {
            Some(name) => name,
            None => return Err("missing_Chromosome"),
        };

        // Try to extract the UMI (UB), and report if missing
        let umi_u64 = if bam_umi_tag == b"No" {
            pseudo_umi
        } else {
            let umi = match Self::get_tag_value(bam_feature, bam_umi_tag) {
                Some(u) => u,
                None => return Err("missing_UMI"),
            };
            IntToStr::new(umi).into_u64()
        };

        // Extract the start position (1-based for your GTF comparison)
        let start = bam_feature.pos(); // This is 0-based, needs adjustment for GTF

        // Create a new ReadData object
        let res = Self::new(
            cell_id,
            umi_u64,
            (start).try_into().unwrap(), // Adjust for GTF 1-based start
            bam_feature.flags(),
            &bam_feature.cigar().to_string(),
            chr,
            bam_feature.seq().as_bytes(),
            bam_feature.qual().to_vec(),
        );

        let id = std::str::from_utf8(bam_feature.qname()).unwrap().to_string(); // Convert read name to string

        Ok((id, res))
    }


    pub fn is(&self, what:&str ) -> bool {
        self.flag.is( what )
    }
}
