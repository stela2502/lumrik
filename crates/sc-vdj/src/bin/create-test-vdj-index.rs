use anyhow::Result;
use clap::Parser;
use sc_vdj::VdjIndexBuilder;
use lumrik_clap::{ReadableFile, OutputFile};

#[derive(Debug, Parser)]
#[command(author, version, name = "create-test-vdj-index")]
struct Cli {
    #[arg(long)]
    gtf: ReadableFile,
    #[arg(long)]
    genome: ReadableFile,
    #[arg(long, value_delimiter = ',')]
    segments: Vec<String>,
    #[arg(long)]
    out: OutputFile,
}

fn main() -> Result<()> {
    let c = Cli::parse();
    c.out.prepare()?;
    let idx = VdjIndexBuilder::default()
        .with_gene_names(c.segments)
        .build(c.gtf.as_path(), c.genome.as_path())?;
    idx.save(&c.out)?;
    println!("{} segments -> {}", idx.len(), c.out.display());
    Ok(())
}
