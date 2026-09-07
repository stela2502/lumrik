use anyhow::{bail, Context, Result};
use clap::Parser;
use sc_vdj::{VdjMapper, VdjMapperConfig, VdjReference, VdjReferenceBuilder};
use std::collections::{BTreeMap, HashSet};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "create-test-vdj-index",
    about = "Create a compact Lumrik .vdjidx containing only selected germline segments",
    after_help = "Example:\n  create-test-vdj-index \\\n    --gtf annotation.gtf \\\n    --genome genome.fa \\\n    --segments Ighv1-64,Ighd1-1,Ighj3,Igha,Igkv6-15,Igkj2,Igkc,Iglv3,Iglj2 \\\n    --out tests/data/vdj-integration-cell/reference.vdjidx"
)]
struct Cli {
    /// Genome annotation containing antigen-receptor segments.
    #[arg(long, value_name = "GTF")]
    gtf: PathBuf,

    /// Genome FASTA matching the annotation.
    #[arg(long, value_name = "FASTA")]
    genome: PathBuf,

    /// Exact V/D/J/C gene names to retain, comma-separated.
    #[arg(long, value_delimiter = ',', num_args = 1.., required = true, value_name = "GENE,...")]
    segments: Vec<String>,

    /// Output compact V(D)J index.
    #[arg(long, value_name = "FILE.vdjidx")]
    out: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let requested: HashSet<String> = cli.segments.into_iter().collect();
    if requested.is_empty() {
        bail!("no V(D)J segment names requested");
    }

    eprintln!("building V(D)J germline reference from annotation...");
    let full_reference = VdjReferenceBuilder::default()
        .build(&cli.gtf, &cli.genome)
        .with_context(|| {
            format!(
                "building V(D)J reference from {} and {}",
                cli.gtf.display(),
                cli.genome.display()
            )
        })?;

    let mut found = HashSet::new();
    let segments = full_reference
        .segments
        .into_iter()
        .filter(|segment| {
            let keep = requested.contains(&segment.name);
            if keep {
                found.insert(segment.name.clone());
            }
            keep
        })
        .collect::<Vec<_>>();

    let mut missing = requested
        .difference(&found)
        .cloned()
        .collect::<Vec<_>>();
    missing.sort();
    if !missing.is_empty() {
        bail!(
            "requested V(D)J segment(s) not found in annotation: {}",
            missing.join(", ")
        );
    }
    if segments.is_empty() {
        bail!("selected V(D)J reference is empty");
    }

    let selected = VdjReference { segments };
    let counts: BTreeMap<_, _> = selected.counts();

    eprintln!(
        "building compact V(D)J index for {} selected reference segments:",
        selected.len()
    );
    for ((chain, kind), count) in counts {
        eprintln!("  {chain} {kind:?}: {count}");
    }
    for segment in &selected.segments {
        eprintln!(
            "  {}\t{}\t{:?}\t{}:{}-{}",
            segment.chain,
            segment.name,
            segment.kind,
            segment.chr,
            segment.start,
            segment.end
        );
    }

    let mapper = VdjMapper::new(selected, VdjMapperConfig::default());
    mapper
        .save_index(&cli.out)
        .with_context(|| format!("writing {}", cli.out.display()))?;

    println!("compact VDJ index v4 written to {}", cli.out.display());
    Ok(())
}
