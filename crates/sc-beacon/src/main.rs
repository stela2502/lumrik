use std::fs;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use scdata::{load_mtx_feature_matrix, read_mtx_cell_ids, MexFeatureIndex};

mod cli;

use cli::GuideModelCli;
use sc_beacon::run_from_scdata;

#[derive(Debug, Parser)]
#[command(
    name = "lumrik-guides",
    about = "Ambient-aware multi-guide caller for 10x CRISPR feature-barcode matrices"
)]
struct Cli {
    /// 10x raw_feature_bc_matrix directory.
    #[arg(long)]
    raw: PathBuf,

    /// 10x filtered_feature_bc_matrix directory.
    #[arg(long)]
    filtered: PathBuf,

    /// Output directory.
    #[arg(long)]
    out: PathBuf,

    /// 10x feature type containing the guide counts.
    #[arg(long, default_value = "CRISPR Guide Capture")]
    feature_type: String,

    /// Number of worker threads.
    #[arg(long)]
    threads: Option<usize>,

    #[command(flatten)]
    model: GuideModelCli,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    fs::create_dir_all(&cli.out)
        .with_context(|| format!("creating {}", cli.out.display()))?;

    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads.max(1))
            .build_global()
            .context("failed to initialize Rayon thread pool")?;
    }

    let threads = cli.threads.unwrap_or_else(rayon::current_num_threads);

    let (raw, raw_index, cell_barcode_len) =
        load_mtx_feature_matrix(&cli.raw, &cli.feature_type, threads)?;
    let filtered_index = MexFeatureIndex::from_dir(&cli.filtered, &cli.feature_type)?;
    raw_index.validate_compatible(&filtered_index)?;
    let filtered_cells = read_mtx_cell_ids(&cli.filtered)?;

    eprintln!("Found {} guide features.", raw_index.features().len());
    eprintln!(
        "Running Beacon on {} retained cells using raw droplets for ambient estimation...",
        filtered_cells.len(),
    );

    let (filtered_counts, background_counts) = raw.split_by_cells(&filtered_cells);

    let mut result = run_from_scdata(
        &filtered_counts,
        &background_counts,
        &filtered_cells,
        cell_barcode_len,
        &raw_index,
        &cli.model.background_config(),
        &cli.model.fit_config(),
        &cli.model.call_config(),
    )?;

    result.write(&cli.out, &raw_index, cell_barcode_len)?;

    Ok(())
}
