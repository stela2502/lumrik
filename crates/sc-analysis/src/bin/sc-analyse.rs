use anyhow::Result;
use clap::{Parser, Subcommand};
use sc_analysis::cli::{AnalysisCli,run_analysis};
use sc_analysis::snp::{SnpCli,run};

#[derive(Parser)]
#[command(name="sc-analyse",version)]
struct Cli {
    #[command(subcommand)] command: Command,
}
#[derive(Subcommand)]
enum Command { Snp(SnpCli) }
fn main() -> Result<()> {
    // Preserve the legacy `sc-analyse INPUT --out ...` interface.
    if std::env::args_os().nth(1).is_some_and(|x|x!="snp" && x!="--help" && x!="-h" && x!="--version") {
        return run_analysis(AnalysisCli::parse());
    }
    match Cli::parse().command {Command::Snp(cli)=>run(cli)}
}
