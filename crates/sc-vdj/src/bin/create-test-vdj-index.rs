use anyhow::Result;
use clap::Parser;
use sc_vdj::VdjIndexBuilder;
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
    let idx = VdjIndexBuilder::default()
        .with_gene_names(c.segments)
        .build(c.gtf, c.genome)?;
    idx.save(&c.out)?;
    println!("{} segments -> {}", idx.len(), c.out.display());
    Ok(())
}
