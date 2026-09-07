use anyhow::{bail, Result};
use clap::Parser;
use sc_vdj::{VdjIndex, VdjIndexBuilder};
use std::collections::HashSet;
use std::path::PathBuf;
#[derive(Debug, Parser)]
#[command(author, version, name = "create-test-vdj-index")]
struct Cli {
    #[arg(long)]
    gtf: PathBuf,
    #[arg(long)]
    genome: PathBuf,
    #[arg(long, value_delimiter = ',')]
    segments: Vec<String>,
    #[arg(long)]
    out: PathBuf,
}
fn main() -> Result<()> {
    let c = Cli::parse();
    let idx = VdjIndexBuilder::default().build(c.gtf, c.genome)?;
    let wanted: HashSet<_> = c.segments.iter().map(String::as_str).collect();
    let kept: Vec<_> = idx
        .segments
        .into_iter()
        .filter(|s| wanted.contains(s.name.as_str()))
        .collect();
    if kept.is_empty() {
        bail!("no requested segments found")
    };
    let out = VdjIndex::from_segments(kept)?;
    out.save(&c.out)?;
    println!("{} segments -> {}", out.len(), c.out.display());
    Ok(())
}
