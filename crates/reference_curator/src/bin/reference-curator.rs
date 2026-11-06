use anyhow::Result;
use lumrik_clap::{ExistingOrNewFile, OutputFile, Parser, Subcommand};
use reference_curator::{CandidateFilter, ReferenceCurator};

#[derive(Parser)]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(
    name = "reference-curator",
    about = "Inspect and export Lumrik's persistent reference curation store"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    Inspect {
        #[arg(long)]
        store: ExistingOrNewFile,
    },
    ExportFasta {
        #[arg(long)]
        store: ExistingOrNewFile,
        #[arg(long)]
        out: OutputFile,
        #[arg(long, conflicts_with = "unresolved_only")]
        resolved_only: bool,
        #[arg(long, conflicts_with = "resolved_only")]
        unresolved_only: bool,
    },
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Inspect { store } => {
            let curator = ReferenceCurator::open(store.as_path())?;
            println!("candidates\t{}", curator.len());
            println!("resolved\t{}", curator.resolved().count());
            println!("unresolved\t{}", curator.unresolved().count());
        }
        Command::ExportFasta {
            store,
            out,
            resolved_only,
            unresolved_only,
        } => {
            let filter = if unresolved_only {
                CandidateFilter::UnresolvedOnly
            } else if resolved_only {
                CandidateFilter::ResolvedOnly
            } else {
                CandidateFilter::All
            };
            out.prepare()?;
            ReferenceCurator::open(store.as_path())?.export_fasta(out.as_path(), filter)?;
        }
    }
    Ok(())
}
