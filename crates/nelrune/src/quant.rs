use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;

use anyhow::{Context, Result};
use bam_tide::{AdditionalFeatureSource, FeatureTagCounts};
use bam_tide::index::{GeneFeatureIndex, TranscriptFeatureIndex};
use bam_tide::quantification::bam_collector::{BamCollector, BamCollectorConfig};
use bam_tide::quantification::cli::QuantMode;
use clap::Parser;
use gtf_splice_index::SpliceIndex;
use sc_primer::PrimerCli;

#[derive(Debug, Parser)]
#[command(about = "Quantify a STAR BAM and integrate FASTQ-derived Nelrune observations")]
pub struct QuantCli {
    #[arg(long, num_args = 1..)]
    bam: Vec<PathBuf>,
    #[arg(long, help = "Directory produced by nelrune prepare-fastqs")]
    prepare: PathBuf,
    #[command(flatten)]
    primer: PrimerCli,
    #[command(flatten)]
    bam_collector: BamCollectorConfig,
    #[arg(long, num_args = 1..)]
    additional_features: Vec<AdditionalFeatureSource>,
    #[arg(long, default_value_t = 4)]
    additional_feature_min_hits: u32,
    #[arg(long)]
    min_cell_counts: Option<usize>,
    #[arg(long, short)]
    outpath: PathBuf,
}

pub fn run() -> Result<()> {
    let args = QuantCli::parse_from(std::env::args().skip(1));
    fs::create_dir_all(&args.outpath)?;

    let primer = args.primer.detector().map_err(anyhow::Error::msg)?;
    let cell_barcode_len = primer.cell_len();
    let manifest = fs::read_to_string(args.prepare.join("prepare-manifest.tsv"))
        .context("reading prepare-manifest.tsv")?;
    if let Some(line) = manifest.lines().find(|line| line.starts_with("cell_barcode_len\t")) {
        let prepared_len: usize = line.split('\t').nth(1).unwrap_or("0").parse()?;
        anyhow::ensure!(prepared_len == cell_barcode_len,
            "prepare/quant cell barcode length mismatch: prepared {prepared_len}, quant {cell_barcode_len}");
    }

    let collector = BamCollector::from_cli(args.bam_collector.clone())?
        .with_grammar(primer.grammar().clone());
    let result = collector.run_paths(&args.bam)?;
    let mut data = result.data;
    let index = SpliceIndex::load(&args.bam_collector.index)?;

    let retained_cells: HashSet<u64> = if let Some(min) = args.min_cell_counts {
        data.cells_with_min_exonic_umis(min)
    } else {
        let calling = data.beacon_cell_calling().map_err(anyhow::Error::msg)?;
        calling.write_tsv(args.outpath.join("cell_calling.tsv")).map_err(anyhow::Error::msg)?;
        calling.write_qc(args.outpath.join("qc")).map_err(anyhow::Error::msg)?;
        calling.retained
    };

    match args.bam_collector.quant_mode {
        QuantMode::Gene => {
            let features = GeneFeatureIndex::new(&index);
            data.write_with_unfiltered_for_cells(&args.outpath, &retained_cells, &features,
                result.snp.as_ref().map(|s| &s.index), Some(cell_barcode_len))
                .map_err(anyhow::Error::msg)?;
        }
        QuantMode::Transcript => {
            let features = TranscriptFeatureIndex::new(&index);
            data.write_with_unfiltered_for_cells(&args.outpath, &retained_cells, &features,
                result.snp.as_ref().map(|s| &s.index), Some(cell_barcode_len))
                .map_err(anyhow::Error::msg)?;
        }
    }

    let observation_path = args.prepare.join("feature_observations.bin");
    if observation_path.is_file() {
        let mut feature_counts = FeatureTagCounts::load_observations(
            &observation_path, &args.additional_features, args.additional_feature_min_hits)?;
        // Preserve the existing Nelrune integration point: canonical GEX cells
        // are decided from BAM quantification, then applied to FASTQ observations.
        feature_counts.finalize_and_write(&retained_cells, cell_barcode_len, &args.outpath)?;
    }

    fs::write(args.outpath.join("nelrune-report.txt"), data.report.to_string())?;
    eprintln!("[nelrune quant] complete: {} canonical cells", retained_cells.len());
    Ok(())
}
