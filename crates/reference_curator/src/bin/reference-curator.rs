use anyhow::Result;
use clap::{Parser, Subcommand};
use reference_curator::{CandidateFilter, ReferenceCurator};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "reference-curator", about = "Inspect and export Lumrik's persistent reference curation store")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    Inspect { #[arg(long)] store: PathBuf },
    ExportFasta {
        #[arg(long)] store: PathBuf,
        #[arg(long)] out: PathBuf,
        #[arg(long, conflicts_with = "unresolved_only")]
        resolved_only: bool,
        #[arg(long, conflicts_with = "resolved_only")]
        unresolved_only: bool,
    },
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Inspect { store } => {
            let curator = ReferenceCurator::open(store)?;
            println!("candidates\t{}", curator.len());
            println!("resolved\t{}", curator.resolved().count());
            println!("unresolved\t{}", curator.unresolved().count());
        }
        Command::ExportFasta { store, out, resolved_only, unresolved_only } => {
            let filter = if unresolved_only {
                CandidateFilter::UnresolvedOnly
            } else if resolved_only {
                CandidateFilter::ResolvedOnly
            } else {
                CandidateFilter::All
            };
            ReferenceCurator::open(store)?.export_fasta(out, filter)?;
        }
    }
    Ok(())
}
