//! Plotting entry points for Ommverse reference-model visualizations.
//!
//! Keep plot families behind this module so the CLI can grow beyond genomic
//! gene/transcript views without coupling plotting to the command parser.

mod genes_transcripts;

pub use genes_transcripts::{
    GeneTranscriptPlotOptions, plot_gene_transcripts, plot_gene_transcripts_with_options,
};

use std::path::{Path, PathBuf};

/// Resolve a multi-gene plot output prefix deterministically.
///
/// Unlike the former implementation, the meaning never depends on whether
/// the path happens to exist as a directory.
pub fn gene_plot_output_path(out: &Path, gene: &str) -> PathBuf {
    let mut path = out.as_os_str().to_os_string();
    path.push(gene);
    path.push(".svg");
    PathBuf::from(path)
}

#[cfg(test)]
mod tests {
    use super::gene_plot_output_path;
    use std::path::Path;

    #[test]
    fn genome_plot_output_is_always_a_prefix() {
        assert_eq!(
            gene_plot_output_path(Path::new("results/TP53_"), "TP53"),
            Path::new("results/TP53_TP53.svg")
        );
        assert_eq!(
            gene_plot_output_path(Path::new("results/GENE_"), "TP53"),
            Path::new("results/GENE_TP53.svg")
        );
    }
}
