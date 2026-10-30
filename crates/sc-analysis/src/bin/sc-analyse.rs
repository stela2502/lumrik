use anyhow::Result;
use clap::{Parser, Subcommand};
use sc_analysis::annotate::{AnnotateCli, run_annotate};
use sc_analysis::cli::{AnalysisCli, run_analysis};
use sc_analysis::snp::{SnpCli, run};
use sc_analysis::split::{SplitCli, run_split};

#[derive(Parser)]
#[command(name = "sc-analyse", version, about = "Lumrik single-cell analysis")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Analyze single-cell gene expression data.
    Gex(AnalysisCli),
    /// Analyze single-cell SNP evidence.
    Snp(SnpCli),
    /// Import annotation table columns into H5AD obs.
    Annotate(AnnotateCli),
    /// Split an H5AD by values in an obs column.
    Split(SplitCli),
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Gex(cli) => run_analysis(cli),
        Command::Snp(cli) => run(cli),
        Command::Annotate(cli) => run_annotate(cli),
        Command::Split(cli) => run_split(cli),
    }
}
