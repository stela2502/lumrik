pub mod split;
pub mod annotate;
mod bioinformatics;
pub mod cli;
pub mod module_score;
pub mod cell_cycle;
mod cluster;
mod data;
mod embedding;
mod leiden;
mod mex;
mod multiplet;
mod normalize;
mod norn;
mod partition;
mod pca;
mod qc;
mod report;
pub mod snp;
mod stats;
pub use bioinformatics::{
    AnalysisConfig, AnalysisSummary, analyze_exon_matrix, analyze_sample_matrices,
    merge_sample_inputs,
};
pub use data::{CellAnnotations, SingleCellData};
pub use report::write_group_stats;

/// Describe the analysis operations currently exposed to an AI investigator.
/// This is a capability declaration, not a claim about any loaded dataset.
pub fn contains() -> String {
    "sc-analysis capabilities: load pre-analysed expression.mtx.gz/genes.tsv.gz/cells.tsv.gz; group cells by an annotation column; retrieve per-gene detection fractions; analyze raw exon MEX using normalization, PCA, clustering and reporting. H5AD ingestion is provided by sc_anndata through an adapter. Do not infer mutation status, causal trajectories, or biological replication from clusters."
        .to_string()
}

#[cfg(test)]
mod capability_tests {
    #[test]
    fn describes_real_capabilities_and_limits() {
        let text = super::contains();
        assert!(text.contains("group cells"));
        assert!(text.contains("H5AD"));
        assert!(text.contains("Do not infer"));
    }
}
