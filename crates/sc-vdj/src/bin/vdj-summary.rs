use anyhow::{bail, Context, Result};
use clap::Parser;
use sc_vdj::{NelruneIdentityResolver, VdjIndex, VdjIndexBuilder, VdjRunner, VdjRunnerConfig};
use std::path::PathBuf;
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "vdj-summary",
    about = "Print per-cell sc-vdj recombination summaries"
)]
struct Cli {
    #[arg(long)]
    bam: PathBuf,
    #[arg(long)]
    index: Option<PathBuf>,
    #[arg(long)]
    gtf: Option<PathBuf>,
    #[arg(long)]
    genome: Option<PathBuf>,
    #[arg(long)]
    out: Option<PathBuf>,
    #[arg(long)]
    exonic: Option<PathBuf>,
    #[arg(long, default_value_t = 12)]
    min_sequence_overlap: usize,
}
fn main() -> Result<()> {
    let c = Cli::parse();
    let index = match (&c.index, &c.gtf, &c.genome) {
        (Some(p), None, None) => {
            VdjIndex::load(p).with_context(|| format!("loading {}", p.display()))?
        }
        (None, Some(g), Some(f)) => VdjIndexBuilder::default().build(g, f)?,
        _ => bail!("provide either --index or --gtf + --genome"),
    };
    let mut runner = VdjRunner::new(
        index,
        VdjRunnerConfig {
            min_sequence_overlap: c.min_sequence_overlap,
        },
    );
    runner.read_bam(&c.bam, &NelruneIdentityResolver)?;
    for (cell, rs) in runner.identify() {
        let name = runner
            .cell_names
            .get(&cell)
            .map(String::as_str)
            .unwrap_or("unknown");
        for r in rs {
            let v = &runner.index.segment(r.v).unwrap().name;
            let d =
                r.d.and_then(|x| runner.index.segment(x))
                    .map(|x| x.name.as_str())
                    .unwrap_or("");
            let j = &runner.index.segment(r.j).unwrap().name;
            let cc = r
                .constant
                .as_ref()
                .and_then(|x| runner.index.segment(x.segment))
                .map(|x| x.name.as_str())
                .unwrap_or("");
            println!(
                "{name}\t{}\tV={v}\tD={d}\tJ={j}\tC={cc}\t{}",
                r.chain, r.stable_id
            )
        }
    }
    Ok(())
}
