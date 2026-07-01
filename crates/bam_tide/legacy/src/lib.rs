
pub mod atac_to_rna_mapper;
pub mod bed_data;
pub mod data_iter;
pub mod bam_flag;

pub mod gtf;
pub mod gtf_logics;
pub mod mutation_processor;

pub mod read_data;

pub mod feature_matcher;
pub use crate::feature_matcher::QueryErrors; 
pub use crate::atac_to_rna_mapper::ATACtoRNAMapper;

use directories::ProjectDirs;
use std::path::PathBuf;


pub fn get_data_dir() -> PathBuf {
    ProjectDirs::from("org", "MyOrg", "bam_tide")
        .expect("Could not determine user data directory")
        .data_dir()
        .to_path_buf()
}
