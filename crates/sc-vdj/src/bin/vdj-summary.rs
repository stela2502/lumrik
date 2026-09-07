#[path = "vdj-summary/evidence_report.rs"]
mod evidence_report;

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

    /// Print raw, pre-merge receptor evidence as an ASCII alignment canvas.
    /// Output goes to stderr so the normal TSV-like summary on stdout stays clean.
    #[arg(long)]
    visual_evidence: bool,

    /// Print a human-readable audit of BAM evidence, direct fragment links, and final calls.
    #[arg(long)]
    evidence_report: bool,

    /// Re-scan the BAM with cell-specific CDR3/J baits and constant-region baits,
    /// then use direct same-read or same-fragment evidence to fill missing C calls.
    #[arg(long)]
    rescan_recombination_evidence: bool,

    /// Restrict --visual-evidence to one cell barcode/name.
    #[arg(long, requires = "visual_evidence")]
    visual_evidence_cell: Option<String>,
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

    if c.visual_evidence {
        for (cell_id, evidence) in runner.evidence.cells() {
            let name = runner
                .cell_names
                .get(&cell_id)
                .map(String::as_str)
                .unwrap_or("unknown");
            if c.visual_evidence_cell.as_deref().is_some_and(|wanted| wanted != name) {
                continue;
            }
            eprintln!("\nCELL {name} ({cell_id})");
            for chain in evidence.chains(&runner.index) {
                eprintln!(
                    "{}",
                    evidence.display_raw_chain(
                        &runner.index,
                        chain,
                        c.min_sequence_overlap,
                    )
                );
                eprintln!(
                    "{}",
                    evidence.display_summarized_chain(
                        &runner.index,
                        chain,
                        c.min_sequence_overlap,
                    )
                );
            }
        }
    }

    let mut calls = runner.identify();
    let initial_calls = c
        .rescan_recombination_evidence
        .then(|| calls.clone());
    let rescan_report = if c.rescan_recombination_evidence {
        Some(runner.rescue_missing_constants_from_bam_with_report(
            &c.bam,
            &NelruneIdentityResolver,
            &mut calls,
        )?)
    } else {
        None
    };

    if c.evidence_report {
        eprintln!("{}", evidence_report::describe_vdj_run(&c.bam, &runner, &calls)?);
    }
    if let (Some(before), Some(report)) = (initial_calls.as_ref(), rescan_report.as_ref()) {
        eprintln!(
            "{}",
            evidence_report::describe_recombination_rescan(&runner, before, &calls, report)?
        );
    }

    for (cell, rs) in calls {
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
