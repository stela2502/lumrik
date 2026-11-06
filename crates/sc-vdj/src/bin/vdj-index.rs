use anyhow::{Context, Result};
use clap::Parser;
use sc_vdj::VdjIndexBuilder;
use lumrik_clap::{ReadableFile, OutputFile};
#[derive(Debug, Parser)]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(
    author,
    version,
    name = "vdj-index",
    about = "Build a compact V/D/J/C index for sc-vdj"
)]
struct Cli {
    #[arg(long)]
    gtf: ReadableFile,
    #[arg(long)]
    genome: ReadableFile,
    #[arg(long)]
    out: OutputFile,
}
fn main() -> Result<()> {
    let c = Cli::parse();
    c.out.prepare()?;
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
