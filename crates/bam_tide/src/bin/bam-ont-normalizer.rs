use anyhow::Result;
use bam_tide::ont_normalizer::{OntNormalizer, cli::Cli};

fn main() -> Result<()> {
    let cli = Cli::parse_args();
    let mut normalizer = OntNormalizer::from_cli(cli)?;
    normalizer.run()?;

    normalizer.write_feature_tag_table_if_present()?;

    eprintln!("{}", normalizer.stats_report());
    let log_path = normalizer.config().out.with_extension("log.tsv");
    normalizer.stats().report_to_csv(log_path.to_str().unwrap());

    Ok(())
}
