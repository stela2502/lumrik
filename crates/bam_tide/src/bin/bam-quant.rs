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
    cli::{CellCallingMode, QuantCli, QuantMode},
};
use gtf_splice_index::SpliceIndex;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Quantify BAM files through bam-tide BamCollector"
)]
struct Cli {
    #[command(flatten)]
    quant: QuantCli,
}

fn main() -> Result<()> {
    let args = Cli::parse().quant;

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

        read_tags: args.read_tags.clone(),
        analysis_type: args.analysis_type,
        cell_tag: args.cell_tag,
        umi_tag: args.umi_tag,
        grammar_type: args.grammar_type,

        /*
         * bam-quant consumes existing BAM files.
         * It does not need to emit another BAM.
         */
        bam_out: None,
    };

    let mut collector = BamCollector::from_cli(config)?;
    if let Some(structure) = args.primer_structure.as_deref() {
        let grammar = sc_primer::Grammar::parse("bam-quant", structure)
            .map_err(anyhow::Error::msg)
            .context("parsing --primer-structure")?;
        collector = collector.with_grammar(grammar);
    }

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
            write_quantification(&mut data, &args, &features, snp.as_ref().map(|s| &s.index))
                .context("writing gene quantification")?;
        }

        QuantMode::Transcript => {
            let features = TranscriptFeatureIndex::new(&idx);
            write_quantification(&mut data, &args, &features, snp.as_ref().map(|s| &s.index))
                .context("writing transcript quantification")?;
        }
    }

    data.report.stop_file_io_time();

    println!("{}", data.report);

    Ok(())
}

fn write_quantification<T, F>(
    data: &mut bam_tide::QuantData,
    args: &QuantCli,
    features: &T,
    snp_index: Option<&F>,
) -> Result<()>
where
    T: scdata::FeatureIndex,
    F: scdata::FeatureIndex,
{
    match args.cell_calling {
        CellCallingMode::Fixed => {
            let (retained, accounting) = data
                .write_with_unfiltered(
                    &args.outpath,
                    args.min_cell_counts,
                    features,
                    snp_index,
                    None,
                )
                .map_err(anyhow::Error::msg)?;
            println!(
                "{accounting}Cell calling\n------------\nmethod: fixed\nminimum UMIs: {}\nretained cells: {}\nremoved by exonic cutoff: {}",
                args.min_cell_counts,
                retained.len(),
                accounting.exonic_cells.saturating_sub(retained.len())
            );
        }
        CellCallingMode::Beacon => {
            println!("Running sc-beacon barcode-rank knee cell identification and QC...");
            let calling = data.beacon_cell_calling().map_err(anyhow::Error::msg)?;
            let retained_count = calling.retained.len();

            std::fs::create_dir_all(&args.outpath)
                .with_context(|| format!("creating {}", args.outpath.display()))?;
            calling.write_tsv(args.outpath.join("cell_calling.tsv")).map_err(anyhow::Error::msg)?;
            calling.write_qc(args.outpath.join("qc")).map_err(anyhow::Error::msg)?;

            let accounting = data.write_with_unfiltered_for_cells(
                &args.outpath, &calling.retained, features, snp_index, None,
            ).map_err(anyhow::Error::msg)?;

            println!(
                "sc-beacon cell identification complete\n{accounting}Cell calling\n------------\nmethod: sc-beacon barcode-rank knee\ncandidate barcodes: {}\ninformative barcodes (>1 UMI): {}\nknee rank: {}\nUMI cutoff: {}\nknee score: {:.6}\nretained cells: {}\nnot called: {}\ncell diagnostics: {}\ncell QC plots: {}",
                calling.fit.candidate_barcodes,
                calling.fit.informative_barcodes,
                calling.fit.knee_rank,
                calling.fit.umi_cutoff,
                calling.fit.score,
                retained_count,
                accounting.exonic_cells.saturating_sub(retained_count),
                args.outpath.join("cell_calling.tsv").display(),
                args.outpath.join("qc").display(),
            );
        }
    }
    Ok(())
}

fn configure_rayon(threads: usize) {
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
}
