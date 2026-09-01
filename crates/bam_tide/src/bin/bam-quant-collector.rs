//! Quantify one or more BAM files using bam-tide's `BamCollector`.
//!
//! The binary is intentionally orchestration-only.
//! BAM parsing, reference loading, SNP handling, job construction,
//! chunk processing, and quantification are owned by `BamCollector`.

use anyhow::{Context, Result};

use clap::Parser;

use bam_tide::index::{GeneFeatureIndex, TranscriptFeatureIndex};

use bam_tide::quantification::{
    bam_collector::{BamCollector, BamCollectorConfig},
    cli::{QuantCli, QuantMode},
};
use gtf_splice_index::SpliceIndex;

fn main() -> Result<()> {
    let args = QuantCli::parse();

    run(args)
}

fn run(args: QuantCli) -> Result<()> {
    configure_rayon(args.threads);

    let config = BamCollectorConfig {
        index: args.index.clone(),

        genome: args.genome.clone(),

        vcf: args.vcf.clone(),

        quant_mode: args.quant_mode,

        min_mapq: args.min_mapq,

        max_reads: args.max_reads,

        read1_only: args.read1_only,

        no_genome_refine: args.no_genome_refine,

        require_strand: args.require_strand,

        require_exact_junction_chain: args.require_exact_junction_chain,

        max_5p_overhang_bp: args.max_5p_overhang_bp,

        max_3p_overhang_bp: args.max_3p_overhang_bp,

        allowed_intronic_gap_size: args.allowed_intronic_gap_size,

        snp_min_anchor: args.snp_min_anchor,

        /*
         * bam-quant consumes existing BAM files.
         * It does not need to emit another BAM.
         */
        bam_out: None,
    };

    let collector = BamCollector::from_cli(config)?;

    let result = collector
        .run_paths(&args.bam)
        .context("collecting BAM quantification")?;

    let mut data = result.data;

    let snp = result.snp;

    /*
     * Export still needs the splice index in order to create
     * the appropriate FeatureIndex.
     *
     * BamCollector already loaded this index internally; we can
     * remove this second load later when its result exposes the
     * export index cleanly.
     */
    let idx = SpliceIndex::load(&args.index)
        .with_context(|| format!("reading splice index {} for export", args.index.display()))?;

    println!("Writing outfiles");

    match args.quant_mode {
        QuantMode::Gene => {
            let features = GeneFeatureIndex::new(&idx);

            data.write(
                &args.outpath,
                args.min_cell_counts,
                &features,
                snp.as_ref().map(|s| &s.index),
            )
            .map_err(anyhow::Error::msg)
            .context("writing gene quantification")?;
        }

        QuantMode::Transcript => {
            let features = TranscriptFeatureIndex::new(&idx);

            data.write(
                &args.outpath,
                args.min_cell_counts,
                &features,
                snp.as_ref().map(|s| &s.index),
            )
            .map_err(anyhow::Error::msg)
            .context("writing transcript quantification")?;
        }
    }

    data.report.stop_file_io_time();

    println!("{}", data.report);

    Ok(())
}

fn configure_rayon(threads: usize) {
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
}
