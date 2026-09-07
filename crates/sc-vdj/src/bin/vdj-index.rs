use anyhow::{Context, Result};
use clap::Parser;
use sc_vdj::VdjIndexBuilder;
use std::path::PathBuf;
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "vdj-index",
    about = "Build a compact V/D/J/C index for sc-vdj"
)]
struct Cli {
    #[arg(long)]
    gtf: PathBuf,
    #[arg(long)]
    genome: PathBuf,
    #[arg(long)]
    out: PathBuf,
}
fn main() -> Result<()> {
    let c = Cli::parse();
    let idx = VdjIndexBuilder::default()
        .build(&c.gtf, &c.genome)
        .context("building VDJ index")?;
    idx.save(&c.out)?;
    println!("VDJ index: {} segments -> {}", idx.len(), c.out.display());
    for ((chain, kind), n) in idx.counts() {
        println!("  {chain} {kind:?}: {n}")
    }
    Ok(())
}
