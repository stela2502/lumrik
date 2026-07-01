use anyhow::Result;
use bam_tide::primer_restore::cli::Cli;
use bam_tide::primer_restore::restore::PrimerRestore;

fn main() -> Result<()> {
    let cli = Cli::parse_args();

    let mut restore = PrimerRestore::from_cli(cli)?;
    restore.run()
}