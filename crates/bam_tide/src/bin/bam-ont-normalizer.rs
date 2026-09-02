use anyhow::Result;
use clap::Parser;
use bam_tide::ont_normalizer::{OntNormalizer, cli::Cli};

fn main() -> Result<()> {
    let cli = Cli::parse();
    let mut normalizer = OntNormalizer::from_cli(cli)?;
    normalizer.run()?;

    eprintln!("{}", normalizer.stats_report());
    let log_path = normalizer.config().out.with_extension("log.tsv");
    normalizer.stats().report_to_csv(log_path.to_str().unwrap());

    Ok(())
}
