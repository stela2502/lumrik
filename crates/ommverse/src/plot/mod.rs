//! Plotting entry points for Ommverse reference-model visualizations.
//!
//! Keep plot families behind this module so the CLI can grow beyond genomic
//! gene/transcript views without coupling plotting to the command parser.

mod genes_transcripts;

pub use genes_transcripts::{
    GeneTranscriptPlotOptions, plot_gene_transcripts, plot_gene_transcripts_with_options,
};

use std::path::{Path, PathBuf};

/// Resolve a multi-gene plot output target. Existing directories and paths
/// ending in a separator produce <gene>.svg inside that directory; every
/// other path is treated as a filename prefix.
pub fn gene_plot_output_path(out: &Path, gene: &str) -> PathBuf {
    let text = out.to_string_lossy();
    if out.is_dir() || text.ends_with(std::path::MAIN_SEPARATOR) || text.ends_with('/') {
        out.join(format!("{gene}.svg"))
    } else {
        PathBuf::from(format!("{text}{gene}.svg"))
    }
}
