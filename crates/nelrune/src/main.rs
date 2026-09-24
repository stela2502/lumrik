mod prepare_fastqs;
mod quant;
use std::fs;
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{CommandFactory, FromArgMatches};
use rust_htslib::bam;

use bam_tide::FeatureTagCounts;
use bam_tide::illumina_normalizer::cli::{InsertRead, PrimerRead};
use bam_tide::illumina_normalizer::{IlluminaNormalizer, IlluminaNormalizerConfig};
use bam_tide::index::{GeneFeatureIndex, TranscriptFeatureIndex};
use bam_tide::ont_normalizer::OntNormalizer;
use bam_tide::ont_normalizer::normalizer::OntNormalizerConfig;
use bam_tide::quantification::bam_collector::BamCollector;
use bam_tide::quantification::cli::QuantMode;

use gtf_splice_index::SpliceIndex;
use sc_mapper::{MappingCall, StreamingMapper};

use lumrik_status::{public_hostname, spawn_status_server};
use nelrune::cli::Cli;
use nelrune::progress::RunProgress;

fn main() -> Result<()> {
    match std::env::args().nth(1).as_deref() {
        Some("prepare-fastqs") => return prepare_fastqs::run(),
        Some("quant") => return quant::run(),
        _ => {}
    }

    let command = Cli::command();
    let matches = command.clone().get_matches();

    let args = Cli::from_arg_matches(&matches)?;

    args.validate()?;
    Cli::print_overrides(&command, &matches);

    run(args)
}

fn run(args: Cli) -> Result<()> {
    fs::create_dir_all(&args.outpath)
        .with_context(|| format!("creating output directory {}", args.outpath.display()))?;

    configure_rayon(args.threads);

    let mut progress = RunProgress::new();
    progress.start_timer("nelrune/startup");

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
        eprintln!("[nelrune] health server: {url}");
        Some(server)
    };

    progress.stage("loading chemistry");
    let primer = args
        .primer
        .detector()
        .map_err(anyhow::Error::msg)
        .context("configuring sc_primer chemistry")?;

    let cell_barcode_len = primer.cell_len();

    progress.stage("starting mapper");
    let mut mapper = args
        .mapper
        .from_cli()
        .context("starting streaming mapper")?;
    let mapper_bam = args
        .bam_collector
        .bam_out
        .clone()
        .unwrap_or_else(|| args.outpath.join("nelrune.mapper.bam"));

    if let Some(parent) = mapper_bam.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("creating mapper BAM directory {}", parent.display()))?;
    }

    let mut sink = MapperBamSink::new(mapper_bam.clone());

    // Spawning STAR is cheap; loading its index is not.  The first progress
    // sample below is therefore the useful boundary between mapper warmup and
    // steady processing.
    progress.stop_timer("nelrune/startup");
    progress.start_timer("nelrune/input_and_mapper");

    progress.stage("waiting for mapper / first batch");
    let (submitted, mut feature_counts) = if let Some(bam_input) = &args.bam {
        run_ont_input(
            &args,
            &primer,
            bam_input,
            &mut mapper,
            &mut sink,
            &mut progress,
        )?
    } else {
        run_illumina_input(&args, &primer, &mut mapper, &mut sink, &mut progress)?
    };

    progress.clear_input_file();

    progress.stop_timer("nelrune/input_and_mapper");

    // Additional-feature reads are filtered out by bam-tide before the mapper
    // Tey are also handles by the normalizers.

    if submitted == 0 {
        bail!("normalization emitted no mapper-ready reads");
    }

    progress.stage("finishing mapper");
    progress.start_timer("nelrune/mapper_finish");

    sink.ensure_open_blocking(&mut mapper)?;
    sink.drain_ready(&mut mapper)?;

    let remaining = mapper.finish().context("finishing streaming mapper")?;
    sink.write_calls(remaining)?;
    sink.finish()?;
    progress.stop_timer("nelrune/mapper_finish");

    progress.stage("quantifying mapper BAM");
    progress.start_timer("nelrune/quantification");

    let collector = BamCollector::from_cli(args.bam_collector.clone())
        .context("configuring BAM collector")?
        .with_grammar(primer.grammar().clone());
    let result = collector
        .run_paths(std::slice::from_ref(&mapper_bam))
        .context("collecting BAM quantification")?;
    let mut data = result.data;
    progress.stop_timer("nelrune/quantification");

    let index = SpliceIndex::load(&args.bam_collector.index).with_context(|| {
        format!(
            "reading splice index {} for export",
            args.bam_collector.index.display()
        )
    })?;

    progress.stage("writing quantification");
    progress.start_timer("nelrune/writing");

    let (retained_cells, calling) = if let Some(min_umi_count) = args.min_umi_count {
        println!("Applying user-defined cell UMI cutoff: {}", min_umi_count);
        (data.cells_with_min_exonic_umis(min_umi_count), None)
    } else {
        println!("Running sc-beacon barcode-rank knee cell identification and QC...");

        let calling = data.beacon_cell_calling().map_err(anyhow::Error::msg)?;

        std::fs::create_dir_all(&args.outpath)
            .with_context(|| format!("creating {}", args.outpath.display()))?;

        calling
            .write_tsv(args.outpath.join("cell_calling.tsv"))
            .map_err(anyhow::Error::msg)?;

        calling
            .write_qc(args.outpath.join("qc"))
            .map_err(anyhow::Error::msg)?;

        let retained = calling.retained.clone();
        (retained, Some(calling))
    };
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
            .map_err(anyhow::Error::msg)
            .context("writing gene quantification")?
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
            .map_err(anyhow::Error::msg)
            .context("writing transcript quantification")?
        }
    };

    // Additional-feature tables have their own feature-type-aware writer.
    // Preserve every observed feature-tag cell under unfiltered/, then retain
    // only the canonical GEX cells in the normal output.
    feature_counts
        .finalize_and_write_all(cell_barcode_len, &args.outpath.join("unfiltered"))
        .context("writing unfiltered additional feature tables")?;
    feature_counts
        .finalize_and_write(&retained_cells, cell_barcode_len, &args.outpath)
        .context("writing additional feature tables")?;
    progress.stop_timer("nelrune/writing");

    let bam_records_seen = data.report.get_issue_count("bam_records_seen");
    let quantified_bam_records = data.report.get_issue_count("cell_umi from BAM tags");
    let compatible_bam_records = data.report.get_issue_count("Compatible");
    let unmapped_bam_records = data.report.get_issue_count("unmapped");
    progress.set_quantification_summary(
        bam_records_seen,
        quantified_bam_records,
        compatible_bam_records,
        unmapped_bam_records,
        retained_cells.len(),
    );

    if let Some(calling) = &calling {
        progress.report_block(
            "Cell calling",
            &format!(
                "method: sc-beacon barcode-rank knee\ncandidate barcodes: {}\ninformative barcodes (>1 UMI): {}\nknee rank: {}\nUMI cutoff: {}\nknee score: {:.6}\nretained cells: {}\nnot called: {}\ncell diagnostics: {}\ncell QC plots: {}\n\n{}",
                calling.fit.candidate_barcodes,
                calling.fit.informative_barcodes,
                calling.fit.knee_rank,
                calling.fit.umi_cutoff,
                calling.fit.score,
                retained_cells.len(),
                cell_accounting.exonic_cells.saturating_sub(retained_cells.len()),
                args.outpath.join("cell_calling.tsv").display(),
                args.outpath.join("qc").display(),
                cell_accounting,
            ),
        );
    } else {
        let cutoff = args
            .min_umi_count
            .expect("fixed cell calling requires --min-umi-count");
        progress.report_block(
            "Cell calling",
            &format!(
                "method: user-defined UMI cutoff\nUMI cutoff: {}\nretained cells: {}\nnot called: {}\n\n{}",
                cutoff,
                retained_cells.len(),
                cell_accounting.exonic_cells.saturating_sub(retained_cells.len()),
                cell_accounting,
            ),
        );
    }

    // Preserve Nelrune's broad orchestration timings in the final report too.
    data.report.merge(progress.mapping_info());

    progress.report_timings();
    progress.report_block("Quantification report", &data.report);

    fs::write(
        args.outpath.join("nelrune-report.txt"),
        data.report.to_string(),
    )
    .context("writing nelrune-report.txt")?;

    progress.stage(format!("mapper BAM retained at {}", mapper_bam.display()));

    progress.finish();
    progress.write_final_status(&args.outpath)?;

    eprintln!(
        "[nelrune] complete: {} reads in {:.1}s ({:.0} reads/s overall, {:.0} reads/s steady after first progress)",
        progress.reads_seen(),
        progress.elapsed().as_secs_f64(),
        progress.average_reads_per_second(),
        progress.processing_reads_per_second(),
    );

    Ok(())
}

fn run_illumina_input(
    args: &Cli,
    primer: &sc_primer::PrimerDetector,
    mapper: &mut StreamingMapper,
    sink: &mut MapperBamSink,
    progress: &mut RunProgress,
) -> Result<(usize, FeatureTagCounts)> {
    let mut submitted = 0usize;

    let config = IlluminaNormalizerConfig {
        out: args.outpath.join("normalized.fastq"),
        read_tags: args.outpath.join("read_tags.tsv"),
        primer_read: PrimerRead::R1,
        insert_read: InsertRead::R2,
        primer: primer.clone(),
        additional_features: args.additional_features.clone(),
        additional_feature_min_hits: args.additional_feature_min_hits,
        min_insert_len: args.min_insert_len,
        max_reads: args.bam_collector.max_reads,
        threads: args.threads,
        gzip_level: 1,
        gzip: false,
    };

    let mut normalizer = IlluminaNormalizer::new(config)?;
    let mut first_progress_seen = false;

    for (r1_path, r2_path) in args.r1.iter().zip(args.r2.iter()) {
        /*
         * Local to this FASTQ pair.
         *
         * --max-reads therefore applies independently to every
         * R1/R2 input pair.
         */

        if first_progress_seen {
            progress.stage(format!("normalizing {}", display_name(r1_path)));
        }

        progress.input_file(r1_path.to_string_lossy());

        normalizer.nelrune_run(
            r1_path,
            r2_path,
            |reads| {
                for (r1, r2) in reads {
                    mapper.submit(r2, r1.as_ref())?;
                }

                submitted += 1;

                sink.drain_ready(mapper)?;

                Ok(true)
            },
            |stats| {
                if !first_progress_seen {
                    first_progress_seen = true;
                    progress.stage(format!("normalizing {}", display_name(r1_path)));
                }
                progress.update_from_mapping_info(stats);
            },
        )?;

        progress.report_block(
            &format!("Illumina normalization: {}", r1_path.display()),
            normalizer.stats(),
        );
    }

    let feature_counts = normalizer.take_feature_tag_counts();

    Ok((submitted, feature_counts))
}

fn run_ont_input(
    args: &Cli,
    primer: &sc_primer::PrimerDetector,
    bam_input: &Path,
    mapper: &mut StreamingMapper,
    sink: &mut MapperBamSink,
    progress: &mut RunProgress,
) -> Result<(usize, FeatureTagCounts)> {
    progress.stage(format!("preparing {}", display_name(bam_input)));

    progress.input_file(bam_input.to_string_lossy());

    let config = OntNormalizerConfig {
        bam: bam_input.to_path_buf(),
        out: args.outpath.join("normalized.fastq"),
        read_tags: args.outpath.join("read_tags.tsv"),
        min_transcript_len: args.min_transcript_len,
        max_reads: args.bam_collector.max_reads,
        primer: primer.clone(),
        additional_features: args.additional_features.clone(),
        additional_feature_min_hits: args.additional_feature_min_hits,
        threads: args.threads,
        gzip_level: 1,
        gzip: false,
    };

    let mut normalizer = OntNormalizer::new(config)?;

    let mut submitted = 0usize;
    let mut first_progress_seen = false;

    normalizer.nelrune_run(
        |reads| {
            for (r1, r2) in reads {
                mapper.submit(r2, r1.as_ref())?;
            }

            submitted += 1;

            sink.drain_ready(mapper)?;

            Ok(true)
        },
        |stats| {
            if !first_progress_seen {
                first_progress_seen = true;
                progress.stage(format!("normalizing {}", display_name(bam_input)));
            }
            progress.update_from_mapping_info(stats);
        },
    )?;

    progress.report_block(
        &format!("ONT normalization: {}", bam_input.display()),
        normalizer.stats(),
    );

    let feature_counts = normalizer.take_feature_tag_counts();

    Ok((submitted, feature_counts))
}

fn display_name(path: &Path) -> String {
    path.file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_else(|| path.to_str().unwrap_or("-"))
        .to_string()
}

struct MapperBamSink {
    path: PathBuf,
    writer: Option<bam::Writer>,
}

impl MapperBamSink {
    fn new(path: PathBuf) -> Self {
        Self { path, writer: None }
    }

    fn drain_ready(&mut self, mapper: &mut StreamingMapper) -> Result<()> {
        if self.writer.is_none() && mapper.header_loaded() {
            self.open_from_mapper(mapper)?;
        }

        if self.writer.is_none() {
            return Ok(());
        }

        while let Some(call) = mapper.try_next()? {
            self.write_call(call)?;
        }

        Ok(())
    }

    fn ensure_open_blocking(&mut self, mapper: &mut StreamingMapper) -> Result<()> {
        if self.writer.is_none() {
            self.open_from_mapper(mapper)?;
        }
        Ok(())
    }

    fn open_from_mapper(&mut self, mapper: &mut StreamingMapper) -> Result<()> {
        let header = mapper
            .header()
            .context("waiting for mapper BAM header")?
            .clone();

        let writer = bam::Writer::from_path(&self.path, &header, bam::Format::Bam)
            .with_context(|| format!("creating mapper BAM {}", self.path.display()))?;

        self.writer = Some(writer);
        Ok(())
    }

    fn write_call(&mut self, call: MappingCall) -> Result<()> {
        let writer = self
            .writer
            .as_mut()
            .context("mapper BAM writer not initialized")?;

        for record in call.records.records {
            let record = record.into_inner();
            writer.write(&record).context("writing mapper BAM record")?;
        }

        Ok(())
    }

    fn write_calls(&mut self, calls: Vec<MappingCall>) -> Result<()> {
        for call in calls {
            self.write_call(call)?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.writer.take();
        Ok(())
    }
}

fn configure_rayon(threads: usize) {
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }
}
