use std::path::PathBuf;

use anyhow::{Context, Result};
use bam_tide::{AdditionalFeatureSource, FeatureTagCounts};
use clap::Parser;

/// Convert the legacy observation-only feature file into the self-contained
/// feature_observations.bin format consumed by current nelrune quant.
#[derive(Debug, Parser)]
#[command(
    name = "convert-feature-observations",
    about = "Convert legacy Nelrune feature_observations.bin to the self-contained format"
)]
struct Cli {
    /// Legacy feature_observations.bin produced by an older prepare-fastqs run.
    #[arg(long)]
    input: PathBuf,

    /// Feature definition(s) used by that original prepare-fastqs run.
    /// Built-ins are bd_sample_mouse and bd_sample_human; FASTA files are also accepted.
    #[arg(long, num_args = 1.., default_value = "bd_sample_mouse")]
    additional_features: Vec<AdditionalFeatureSource>,

    /// FastTagMapper minimum hits used by the original prepare-fastqs run.
    /// This reconstructs the legacy feature dictionary only; it does not remap reads.
    #[arg(long, default_value_t = 4)]
    additional_feature_min_hits: u32,

    /// New self-contained feature observation file. Must differ from --input.
    #[arg(long, short)]
    out: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    anyhow::ensure!(
        cli.input.is_file(),
        "legacy feature observation file does not exist: {}",
        cli.input.display()
    );

    let input = std::fs::canonicalize(&cli.input)
        .with_context(|| format!("resolving {}", cli.input.display()))?;
    if cli.out.exists() {
        let out = std::fs::canonicalize(&cli.out)
            .with_context(|| format!("resolving {}", cli.out.display()))?;
        anyhow::ensure!(
            input != out,
            "refusing to overwrite the legacy input in place"
        );
    }

    eprintln!(
        "[convert-feature-observations] input: {}",
        cli.input.display()
    );
    eprintln!(
        "[convert-feature-observations] feature sources: {:?}",
        cli.additional_features
    );

    let counts = FeatureTagCounts::load_legacy_observations(
        &cli.input,
        &cli.additional_features,
        cli.additional_feature_min_hits,
    )
    .context("loading legacy feature observations")?;

    anyhow::ensure!(
        !counts.is_empty(),
        "legacy feature observation file contained no observations"
    );

    eprintln!(
        "[convert-feature-observations] feature-bearing cells: {}",
        counts.len()
    );

    counts
        .save_observations(&cli.out)
        .with_context(|| format!("writing converted observations to {}", cli.out.display()))?;

    eprintln!(
        "[convert-feature-observations] output: {}",
        cli.out.display()
    );
    Ok(())
}
