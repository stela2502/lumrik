use anyhow::Result;
use bam_tide::primer_restore::cli::Cli;
use bam_tide::primer_restore::restore::PrimerRestore;
use clap::Parser;

fn main() -> Result<()> {
    let cli = Cli::parse();

    let mut restore = PrimerRestore::from_cli(cli)?;
    restore.run()
}
