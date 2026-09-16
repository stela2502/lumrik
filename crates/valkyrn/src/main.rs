use anyhow::Result;
use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "valkyrn",
    about = "Turn reconstructed immune receptors into biological hypotheses"
)]
struct Cli {
    /// nelrune-vdj output directory containing vdj_calls.tsv and airr_rearrangements.tsv.
    #[arg(long)]
    vdj_dir: PathBuf,
    /// Valkyrn output directory.
    #[arg(long, default_value = "valkyrn")]
    out: PathBuf,
    /// Retained for CLI compatibility. Primary clone identity now comes from Lumrik HC:/LC: structural recombination IDs.
    #[arg(long, default_value_t = 1, hide = true)]
    max_cdr3_distance: usize,
    /// Rayon worker threads used for CPU-heavy receptor alignment.
    #[arg(long, default_value_t = 8)]
    threads: usize,
    /// Minimum IGH family size considered for structural prioritization.
    #[arg(long, default_value_t = 3)]
    min_structure_family: usize,
    /// Minimum IGH family size to send through ClonoMap PCA/MST plotting.
    #[arg(long, default_value_t = 100)]
    min_clonomap_family: usize,
    /// Minimum cells sharing one HC+LC structural combination to receive its own rooted ClonoMap.
    #[arg(long, default_value_t = 20)]
    min_clonomap_paired_family: usize,
    /// PCA dimensions retained by ClonoMap for large-family geometry.
    #[arg(long, default_value_t = 30)]
    clonomap_k: usize,
}
fn main() -> Result<()> {
    let c = Cli::parse();
    valkyrn::analyze(
        &c.vdj_dir,
        &c.out,
        c.max_cdr3_distance,
        c.min_structure_family,
        c.threads,
        c.min_clonomap_family,
        c.min_clonomap_paired_family,
        c.clonomap_k,
    )
}
