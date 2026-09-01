// src/lib.rs

pub mod ambient_rna_detect;
pub mod cell_data;
pub mod feature_index;
pub mod mtx;
pub mod sparse_matrix;

pub use crate::cell_data::GeneUmiHash;
pub use crate::feature_index::FeatureIndex;
pub use crate::mtx::{MexFeature, MexFeatureIndex, load_mtx_feature_matrix, read_mtx_cell_ids};
pub use crate::sparse_matrix::{MatrixValueType, Scdata};
