use anyhow::{Context, Result};
use clap::Parser;
use sc_vdj::{VdjMapper, VdjMapperConfig, VdjReferenceBuilder};
use std::path::PathBuf;


#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "vdj-index",
    about = "Compile a reusable Lumrik V(D)J germline + seed index",
    after_help = "The resulting .vdjidx is the preferred reference input for nelrune-vdj and
vdj-summary; it avoids rebuilding the germline and seed indices for every run."
)]
struct Cli {
    /// Genome annotation containing antigen-receptor segments.
    #[arg(long, value_name = "GTF")]
    gtf: PathBuf,

    /// Genome FASTA matching the annotation.
    #[arg(long, value_name = "FASTA")]
    genome: PathBuf,

    /// Output compiled V(D)J index.
    #[arg(long, value_name = "FILE.vdjidx")]
    out: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let gtf = cli.gtf;
    let genome = cli.genome;
    let out = cli.out;
    eprintln!("building V(D)J germline reference...");
    let reference = VdjReferenceBuilder::default().build(&gtf, &genome)?;
    eprintln!("building ambiguity-aware coverage + identity seed indices for {} segments...", reference.len());
    let mapper = VdjMapper::new(reference, VdjMapperConfig::default());
    mapper.save_index(&out).with_context(|| format!("writing {}", out.display()))?;
    println!("VDJ index v4 written to {}", out.display());
    Ok(())
}
