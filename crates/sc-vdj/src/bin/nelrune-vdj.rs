use anyhow::{bail, Context, Result};
use clap::Parser;
use sc_vdj::output::{write_mapping_info, ReportWriter};
use sc_vdj::{NelruneIdentityResolver, VdjIndex, VdjIndexBuilder, VdjRunner, VdjRunnerConfig};
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "nelrune-vdj",
    about = "Reconstruct single-cell V(D)J receptors from a retained Nelrune BAM"
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
    out: PathBuf,
    #[arg(long, default_value_t = 12)]
    min_sequence_overlap: usize,
    #[arg(long)]
    exonic: Option<PathBuf>,
    #[arg(long)]
    write_sequences: bool,
    // Kept for CLI compatibility; live status wiring will be restored separately.
    #[arg(long, default_value_t = 8787)]
    health_port: u16,
    #[arg(long)]
    health_hostname: Option<String>,
    #[arg(long, default_value_t = false)]
    no_health_server: bool,
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
    let n = runner.read_bam(&c.bam, &NelruneIdentityResolver)?;
    let calls = runner.identify();

    let mut writer = ReportWriter::create(&c.out, c.write_sequences)?;
    let mut by_chain = HashMap::<String, usize>::new();
    let mut nr = 0usize;
    for (cell_id, rs) in &calls {
        let name = runner
            .cell_names
            .get(cell_id)
            .map(String::as_str)
            .unwrap_or("unknown");
        writer.write_cell(name, rs, &runner.index)?;
        for r in rs {
            *by_chain.entry(r.chain.to_string()).or_default() += 1;
            nr += 1;
        }
    }
    writer.finish()?;
    write_mapping_info(
        c.out.join("vdj-mapping-info.txt"),
        n,
        calls.len(),
        nr,
        &by_chain,
    )?;

    if !c.no_health_server {
        eprintln!("nelrune-vdj: health-server CLI retained, but live status wiring is not restored yet (requested port {})", c.health_port);
    }
    let _ = (&c.exonic, &c.health_hostname);
    eprintln!("nelrune-vdj: {n} receptor-overlapping BAM records; {} cell(s); {nr} recombination(s); outputs in {}", calls.len(), c.out.display());
    Ok(())
}
