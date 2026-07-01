#![allow(clippy::too_many_arguments)]
// Core infrastructure
pub mod cli;
pub mod core;
pub mod encoder;
pub mod index;
pub mod results;

// Data structures / IO
pub mod bed_data;
pub mod data_iter;
pub mod fastq;
pub mod read_tag_table;
//pub mod tags;

// Quantification
pub mod compare_report;
pub mod quantification;

// Normalizers
pub mod illumina_normalizer;
pub mod ngs_normalizer;
pub mod ont_normalizer;
pub mod primer_restore;

// Converters / utilities
pub mod multi_subset_bam;
pub mod transcriptome_to_genome;

pub use results::QuantData;
