use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bam_tide::FeatureTagCounts;
use bam_tide::index::{GeneFeatureIndex, TranscriptFeatureIndex};
use bam_tide::quantification::bam_collector::{BamCollector, BamCollectorConfig};
use bam_tide::quantification::cli::QuantMode;
use clap::Parser;
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use gtf_splice_index::SpliceIndex;
use lumrik_status::{public_hostname, spawn_status_server};
use nelrune::progress::RunProgress;
use sc_primer::{BdCellVersion, Chemistry, PrimerCli, RhapsodyWhitelist};

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
    /// Minimum unique exonic gene-associated UMIs required to call a barcode a cell.
    /// When omitted, sc-beacon barcode-rank cell calling is used.
    #[arg(long)]
    min_umi_count: Option<usize>,
    /// Include intronic GEX molecules in the canonical expression matrix and cell calling.
    /// The separate intronic matrix is still retained for QC/auditing.
    #[arg(long, default_value_t = false)]
    include_intronic: bool,
    #[arg(long, default_value_t = 8787)]
    health_port: u16,
    #[arg(long)]
    health_hostname: Option<String>,
    #[arg(long, default_value_t = false)]
    no_health_server: bool,
    #[arg(long, short)]
    outpath: PathBuf,
}

pub fn run() -> Result<()> {
    let args = QuantCli::parse_from(std::env::args().skip(1));
    fs::create_dir_all(&args.outpath)?;

    let mut progress = RunProgress::new();
    progress.open_log(args.outpath.join("nelrune.log"))?;
    let _health_server = if args.no_health_server {
        progress.stage("health server disabled");
        None
    } else {
        let health_addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), args.health_port);
        let server = spawn_status_server(progress.state_handle(), health_addr)?;
        let hostname = public_hostname(args.health_hostname.as_deref());
        let url = format!("http://{}:{}", hostname, server.addr().port());
        progress.set_public_url(url.clone());
        eprintln!("[nelrune quant] health server: {url}");
        Some(server)
    };

    progress.stage("loading prepare manifest / chemistry");
    progress.update_memory();
    let primer = args.primer.detector().map_err(anyhow::Error::msg)?;
    let manifest = fs::read_to_string(args.prepare.join("prepare-manifest.tsv"))
        .context("reading prepare-manifest.tsv")?;
    let cell_barcode_len: usize = manifest
        .lines()
        .find_map(|line| line.strip_prefix("cell_barcode_len\t"))
        .context("prepare-manifest.tsv is missing cell_barcode_len")?
        .parse()
        .context("invalid cell_barcode_len in prepare-manifest.tsv")?;

    progress.stage("loading splice index / quantifier");
    progress.update_memory();
    let collector =
        BamCollector::from_cli(args.bam_collector.clone())?.with_grammar(primer.grammar().clone());

    progress.stage("quantifying BAM");
    let result = collector.run_paths_with_progress(&args.bam, |data| {
        progress.update_quantification_live(data);
    })?;
    progress.update_quantification_live(&result.data);
    let mut data = result.data;

    if args.include_intronic {
        // Merge at (cell, feature, UMI) level, so a molecule represented in both layers
        // remains one molecule. Keep `intron` untouched as an auditable evidence layer.
        data.gene.merge(&data.intron);
        progress.stage("including intronic GEX evidence in canonical expression");
    }

    progress.stage("loading output feature index");
    progress.update_memory();
    let index = SpliceIndex::load(&args.bam_collector.index)?;

    // IMPORTANT OUTPUT CONTRACT:
    //
    // Quantification must preserve BOTH views of the GEX evidence:
    //   <outpath>/unfiltered/{exonic,intronic,...}  all observed barcodes (>=1 UMI)
    //   <outpath>/{exonic,intronic,...}             canonical called cells only
    //
    // The unfiltered matrices are not disposable debug output. They are required
    // for auditing cell calling and recovering evidence when a caller is too
    // stringent. Do not replace this with a filtered-only writer.
    progress.stage("calling canonical GEX cells");
    progress.update_memory();
    let (retained_cells, calling) = if let Some(min) = args.min_umi_count {
        eprintln!(
            "[nelrune quant] user-defined UMI cutoff: >= {min}; sc-beacon cell calling disabled"
        );
        (data.cells_with_min_exonic_umis(min), None)
    } else {
        let calling = data.beacon_cell_calling().map_err(anyhow::Error::msg)?;
        calling
            .write_tsv(args.outpath.join("cell_calling.tsv"))
            .map_err(anyhow::Error::msg)?;
        calling
            .write_qc(args.outpath.join("qc"))
            .map_err(anyhow::Error::msg)?;
        let retained = calling.retained.clone();
        (retained, Some(calling))
    };

    progress.stage("writing filtered + unfiltered GEX matrices");
    progress.update_memory();
    let cell_accounting = match args.bam_collector.quant_mode {
        QuantMode::Gene => {
            let features = GeneFeatureIndex::new(&index);
            data.write_with_unfiltered_for_cells(
                &args.outpath,
                &retained_cells,
                &features,
                result.snp.as_ref().map(|s| &s.index),
                Some(cell_barcode_len),
            )
            .map_err(anyhow::Error::msg)?
        }
        QuantMode::Transcript => {
            let features = TranscriptFeatureIndex::new(&index);
            data.write_with_unfiltered_for_cells(
                &args.outpath,
                &retained_cells,
                &features,
                result.snp.as_ref().map(|s| &s.index),
                Some(cell_barcode_len),
            )
            .map_err(anyhow::Error::msg)?
        }
    };

    let matrix_feature_counts = MatrixFeatureCounts::from_output(&args.outpath)?;

    progress.stage("loading feature observations");
    progress.update_memory();
    let observation_path = args.prepare.join("feature_observations.bin");
    if observation_path.is_file() {
        let mut feature_counts = FeatureTagCounts::load_observations(&observation_path)?;
        progress.stage("integrating feature tags / Beacon");
        progress.update_memory();
        // Preserve the existing Nelrune integration point: canonical GEX cells
        // are decided from BAM quantification, then applied to FASTQ observations.
        feature_counts.finalize_and_write(&retained_cells, cell_barcode_len, &args.outpath)?;
    }

    progress.stage("writing numeric cell identifiers");
    progress.update_memory();
    if let Some(whitelist) = rhapsody_whitelist(&args.primer.chemistry) {
        for dir in [
            args.outpath.join("exonic"),
            args.outpath.join("intronic"),
            args.outpath.join("unfiltered/exonic"),
            args.outpath.join("unfiltered/intronic"),
        ] {
            write_numeric_barcodes(&dir, &whitelist)?;
        }
    }

    // Keep the cell-calling decision auditable in the permanent report. In
    // particular, report how many barcodes existed before filtering and how
    // many survive several UMI thresholds; otherwise a surprisingly small
    // canonical set cannot be diagnosed after the run.
    let calling_report = if calling.is_none() {
        let cutoff = args
            .min_umi_count
            .expect("fixed cell calling requires --min-umi-count");

        format!(
            "Cell calling\n\
         ------------\n\
         method: user-defined UMI cutoff\n\
         UMI cutoff: {}\n\
         retained cells: {}\n\
         not called: {}\n\n{}",
            cutoff,
            retained_cells.len(),
            cell_accounting
                .exonic_cells
                .saturating_sub(retained_cells.len()),
            cell_accounting,
        )
    } else {
        let calling = calling
            .as_ref()
            .expect("sc-beacon cell calling requires a calling result");

        format!(
            "Cell calling\n\
         ------------\n\
         method: sc-beacon barcode-rank knee\n\
         candidate barcodes: {}\n\
         informative barcodes (>1 UMI): {}\n\
         knee rank: {}\n\
         UMI cutoff: {}\n\
         knee score: {:.6}\n\
         retained cells: {}\n\
         not called: {}\n\n{}",
            calling.fit.candidate_barcodes,
            calling.fit.informative_barcodes,
            calling.fit.knee_rank,
            calling.fit.umi_cutoff,
            calling.fit.score,
            retained_cells.len(),
            cell_accounting
                .exonic_cells
                .saturating_sub(retained_cells.len()),
            cell_accounting,
        )
    };
    progress.stage("writing final quantification report");
    progress.set_quantification_summary(
        data.report.get_issue_count("bam_records_seen"),
        data.report.get_issue_count("quantified_bam_records"),
        data.report.get_issue_count("compatible"),
        data.report.get_issue_count("unmapped"),
        retained_cells.len(),
    );
    progress.update_quantification_live(&data);
    let report = format!(
        "{}\n\n{}\n\n{}",
        data.report, calling_report, matrix_feature_counts
    );
    fs::write(args.outpath.join("nelrune-report.txt"), report)?;
    eprintln!(
        "[nelrune quant] unfiltered exonic cells: {}",
        cell_accounting.exonic_cells
    );
    eprintln!(
        "[nelrune quant] canonical GEX cells: {}",
        retained_cells.len()
    );
    eprintln!(
        "[nelrune quant] filtered out: {}",
        cell_accounting
            .exonic_cells
            .saturating_sub(retained_cells.len())
    );
    eprintln!(
        "[nelrune quant] unfiltered exonic features written: {}",
        matrix_feature_counts.unfiltered_exonic
    );
    eprintln!(
        "[nelrune quant] filtered exonic features written: {}",
        matrix_feature_counts.filtered_exonic
    );
    eprintln!(
        "[nelrune quant] unfiltered intronic features written: {}",
        matrix_feature_counts.unfiltered_intronic
    );
    eprintln!(
        "[nelrune quant] filtered intronic features written: {}",
        matrix_feature_counts.filtered_intronic
    );
    eprintln!(
        "[nelrune quant] unfiltered matrices: {}",
        args.outpath.join("unfiltered").display()
    );
    eprintln!(
        "[nelrune quant] filtered matrices: {}",
        args.outpath.display()
    );
    progress.finish();
    progress.write_final_status(&args.outpath)?;
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct MatrixFeatureCounts {
    unfiltered_exonic: usize,
    filtered_exonic: usize,
    unfiltered_intronic: usize,
    filtered_intronic: usize,
}

impl MatrixFeatureCounts {
    fn from_output(outpath: &Path) -> Result<Self> {
        Ok(Self {
            unfiltered_exonic: count_gzip_lines(
                &outpath.join("unfiltered/exonic/features.tsv.gz"),
            )?,
            filtered_exonic: count_gzip_lines(&outpath.join("exonic/features.tsv.gz"))?,
            unfiltered_intronic: count_gzip_lines(
                &outpath.join("unfiltered/intronic/features.tsv.gz"),
            )?,
            filtered_intronic: count_gzip_lines(&outpath.join("intronic/features.tsv.gz"))?,
        })
    }
}

impl std::fmt::Display for MatrixFeatureCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Matrix feature counts")?;
        writeln!(f, "---------------------")?;
        writeln!(
            f,
            "unfiltered exonic features written: {}",
            self.unfiltered_exonic
        )?;
        writeln!(
            f,
            "filtered exonic features written: {}",
            self.filtered_exonic
        )?;
        writeln!(
            f,
            "unfiltered intronic features written: {}",
            self.unfiltered_intronic
        )?;
        writeln!(
            f,
            "filtered intronic features written: {}",
            self.filtered_intronic
        )?;
        Ok(())
    }
}

fn count_gzip_lines(path: &Path) -> Result<usize> {
    let file = fs::File::open(path)
        .with_context(|| format!("opening matrix feature file {}", path.display()))?;
    let reader = BufReader::new(GzDecoder::new(file));
    let mut count = 0usize;
    for line in reader.lines() {
        line.with_context(|| format!("reading matrix feature file {}", path.display()))?;
        count += 1;
    }
    Ok(count)
}

fn rhapsody_whitelist(chemistries: &[Chemistry]) -> Option<RhapsodyWhitelist> {
    let version = if chemistries.contains(&Chemistry::BdV2_384) {
        BdCellVersion::V2_384
    } else if chemistries.contains(&Chemistry::BdV2_96) {
        BdCellVersion::V2_96
    } else if chemistries.contains(&Chemistry::BdV1) {
        BdCellVersion::V1
    } else if chemistries.contains(&Chemistry::BdV2_384Vdj) {
        BdCellVersion::V2_384Vdj
    } else {
        return None;
    };
    Some(RhapsodyWhitelist::builtin(version))
}

fn write_numeric_barcodes(dir: &Path, whitelist: &RhapsodyWhitelist) -> Result<()> {
    let input = dir.join("barcodes.tsv.gz");
    if !input.is_file() {
        return Ok(());
    }

    let output = dir.join("barcodes.numeric.tsv.gz");
    let reader = BufReader::new(GzDecoder::new(fs::File::open(&input)?));
    let mut writer = BufWriter::new(GzEncoder::new(
        fs::File::create(&output)?,
        Compression::default(),
    ));

    for line in reader.lines() {
        let barcode = line?;
        let numeric = whitelist
            .cell_id_for_seq(barcode.as_bytes())
            .with_context(|| {
                format!(
                    "BD barcode '{}' from {} has no numeric Rhapsody id",
                    barcode,
                    input.display()
                )
            })?;
        writeln!(writer, "{numeric}")?;
    }
    Ok(())
}
