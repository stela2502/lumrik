//! Lightweight, opinionated single-cell analysis primitives used by Lumrik tools.
mod bioinformatics;
mod cluster;
mod data;
mod embedding;
mod mex;
mod normalize;
mod norn;
mod partition;
mod pca;
mod qc;
mod report;
mod stats;
pub use bioinformatics::{AnalysisConfig, AnalysisSummary, analyze_exon_matrix};
pub use data::{CellAnnotations, SingleCellData};
