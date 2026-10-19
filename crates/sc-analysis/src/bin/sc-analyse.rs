use anyhow::Result;
use clap::{Parser, Subcommand};
use sc_analysis::cli::{AnalysisCli, run_analysis};
use sc_analysis::snp::{SnpCli, run};

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
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Gex(cli) => run_analysis(cli),
        Command::Snp(cli) => run(cli),
    }
}
