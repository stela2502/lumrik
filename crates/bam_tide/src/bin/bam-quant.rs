//! Quantify one or more BAM files using bam-tide's `BamCollector`.
//!
//! The binary is intentionally orchestration-only.
//! BAM parsing, reference loading, SNP handling, job construction,
//! chunk processing, and quantification are owned by `BamCollector`.

use anyhow::{Context, Result};
use std::collections::HashMap;

use clap::Parser;

use bam_tide::quantification::{
    bam_collector::{BamCollector, BamCollectorConfig},
    cli::{CellCallingMode, QuantCli, QuantMode},
};
use gtf_splice_index::{SpliceIndex, SpliceMatchMode};

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
    if args.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

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
    let mut report = result.report;
    let snp = result.snp;

    /*
     * Export still needs the splice index in order to create
     * the appropriate FeatureIndex.
     *
     * BamCollector already loaded this index internally; we can
     * remove this second load later when its result exposes the
     * export index cleanly.
     */
    let match_mode = match args.quant_mode {
        QuantMode::Gene => SpliceMatchMode::Gene,
        QuantMode::Transcript => SpliceMatchMode::Transcript,
    };
    let idx = SpliceIndex::load(&args.index)
        .with_context(|| format!("reading splice index {} for export", args.index.display()))?
        .with_match_mode(match_mode);
    let features = idx.feature_index();

    println!("Writing outfiles");
    write_quantification(&mut data, &args, &features, snp.as_ref().map(|s| &s.index))
        .context("writing quantification")?;

    report.stop_file_io_time();

    println!("{}", report);

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
    let accounting = data.cell_accounting();

    let (retained, calling) = match args.cell_calling {
        CellCallingMode::Fixed => (
            data.cells_with_min_exonic_umis(args.min_umi_count),
            None,
        ),
        CellCallingMode::Beacon => {
            println!("Running sc-beacon barcode-rank knee cell identification and QC...");
            let calling = bam_tide::results::beacon_cell_calling(&data).map_err(anyhow::Error::msg)?;

            std::fs::create_dir_all(&args.outpath)
                .with_context(|| format!("creating {}", args.outpath.display()))?;
            calling
                .write_tsv(args.outpath.join("cell_calling.tsv"))
                .map_err(anyhow::Error::msg)?;
            calling
                .write_qc(args.outpath.join("qc"))
                .map_err(anyhow::Error::msg)?;

            (calling.retained.clone(), Some(calling))
        }
    };

    let mut indexes: HashMap<String, &dyn scdata::FeatureIndex> = HashMap::new();
    indexes.insert(gtf_splice_index::QuantClass::Exonic.as_str().to_string(), features);
    indexes.insert(gtf_splice_index::QuantClass::Intronic.as_str().to_string(), features);
    if let Some(snp_index) = snp_index {
        indexes.insert(scdata::QuantData::SNP_REF.to_string(), snp_index);
        indexes.insert(scdata::QuantData::SNP_ALT.to_string(), snp_index);
    }

    data.write_raw_and_filtered_for_cells(&args.outpath, &retained, &indexes, None)
        .map_err(anyhow::Error::msg)?;

    match calling {
        None => println!(
            "{accounting}Cell calling\n------------\nmethod: fixed\nminimum UMIs: {}\nretained cells: {}\nremoved by exonic cutoff: {}",
            args.min_umi_count,
            retained.len(),
            accounting.exonic_cells.saturating_sub(retained.len())
        ),
        Some(calling) => println!(
            "sc-beacon cell identification complete\n{accounting}Cell calling\n------------\nmethod: sc-beacon barcode-rank knee\ncandidate barcodes: {}\ninformative barcodes (>1 UMI): {}\nknee rank: {}\nUMI cutoff: {}\nknee score: {:.6}\nretained cells: {}\nnot called: {}\ncell diagnostics: {}\ncell QC plots: {}",
            calling.fit.candidate_barcodes,
            calling.fit.informative_barcodes,
            calling.fit.knee_rank,
            calling.fit.umi_cutoff,
            calling.fit.score,
            retained.len(),
            accounting.exonic_cells.saturating_sub(retained.len()),
            args.outpath.join("cell_calling.tsv").display(),
            args.outpath.join("qc").display(),
        ),
    }

    Ok(())
}
