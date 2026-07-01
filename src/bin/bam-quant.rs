//! Quantify a 10x-style BAM against a splice index into `scdata`.
//!
//! Pipeline:
//! 1. Load splice index
//! 2. Optionally load genome FASTA
//! 3. Optionally load SNP VCF
//! 4. Stream BAM into `Job`s
//! 5. Process jobs in rayon chunks
//! 6. Export gene/transcript matrix and optional SNP side-channel matrices
//!
//! Important:
//! - Genome refinement is an addition to the normal quantification path.
//! - SNPs are collected as a side-channel when `--vcf` is supplied.
//! - The binary should stay orchestration-only; quantification logic lives in
//!   `src/quantification/`.

use anyhow::{Context, Result, anyhow};
use clap::Parser;
use rust_htslib::bam::{Read, Reader};
use std::fs::File;
use std::io::Write;

use bam_tide::index::{GeneFeatureIndex, TranscriptFeatureIndex};
use read_tag_table::ReadTagTable;

use bam_tide::quantification::chunk_processor::ChunkProcessor;
use bam_tide::quantification::cli::{QuantCli, QuantMode};
use bam_tide::quantification::job::{Job, JobBuilder};
use bam_tide::quantification::snp::SnpSideChannel;
use bam_tide::quantification::processor_options::ProcessorOptions;
use bam_tide::results::QuantData;

use gtf_splice_index::{MatchOptions, SpliceIndex};
use snp_index::Genome;


const CHUNK: usize = 2_000_000;

fn main() -> Result<()> {
    let args = QuantCli::parse();
    run(args)
}

/// Application runner for `bam-quant`.
///
/// This function intentionally only coordinates the pipeline:
///
/// - parse/load input resources
/// - configure worker objects
/// - stream BAM records into `Job`s
/// - flush chunks into `QuantData`
/// - write final matrices and log
fn run(args: QuantCli) -> Result<()> {
    configure_rayon(args.threads);

    let mut data = QuantData::new();
    data.report = mapping_info::MappingInfo::new(
        None,
        args.min_mapq as f32,
        args.max_reads.unwrap_or(usize::MAX),
    );
    data.report.start_counter();

    let idx = SpliceIndex::load(&args.index)
        .with_context(|| format!("reading index {}", args.index.display()))?;

    println!("{idx}");

    let genome = load_genome(&args)?;

    if args.vcf.is_some() && genome.is_none() {
        return Err(anyhow!("--vcf requires --genome"));
    }

    if !args.read_tags.read_tag_table.is_empty() && args.read_tags.read_tag_table.len() != args.bam.len() {
        return Err(anyhow!(
            "Number of --read-tag-table files ({}) must match number of -b/--bam files ({})",
            args.read_tags.read_tag_table.len(),
            args.bam.len()
        ));
    }

    let first_reader = Reader::from_path(&args.bam[0])
        .with_context(|| format!("bam file could not be read: {}", args.bam[0].display()))?;

    let header = first_reader.header().clone();
    drop(first_reader);

    //let chr_map = idx.chr_map();
    let snp = load_snp_side_channel(&args, &header)?;

    print_run_info(genome.as_ref(), snp.as_ref());

    let match_opts = MatchOptions {
        require_strand: args.require_strand,
        require_exact_junction_chain: args.require_exact_junction_chain,
        max_5p_overhang_bp: args.max_5p_overhang_bp,
        max_3p_overhang_bp: args.max_3p_overhang_bp,
        allowed_intronic_gap_size: args.allowed_intronic_gap_size,
    };
    data.report.stop_file_io_time();

    let processor = ChunkProcessor::new(&idx, snp.as_ref(), match_opts,  ProcessorOptions::from(&args) );

    let mut n_seen = 0usize;

    for (id, bam_path) in args.bam.iter().enumerate() {
        if let Some(max_reads) = args.max_reads
            && n_seen >= max_reads
        {
            break;
        }

        let read_tag_table = if args.read_tags.read_tag_table.is_empty() {
            None
        } else {
            let tag_path = &args.read_tags.read_tag_table[id];

            let config = args.read_tags.to_config_for_id(id)?;

            let tab = ReadTagTable::from_config(&config)
                .with_context(|| {
                    format!(
                        "reading read tag table {} for BAM {}",
                        tag_path.display(),
                        bam_path.display()
                    )
                })?;

            println!("{tab}");
            Some(tab)
        };

        process_bam_file(
            bam_path,
            &header,
            idx.chr_map(),
            genome.as_ref(),
            snp.as_ref(),
            read_tag_table.as_ref(),
            &processor,
            &args,
            &mut data,
            &mut n_seen,
        )?;
    }

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
    write_log(&args, &data)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_bam_file(
    bam_path: &std::path::Path,
    expected_header: &rust_htslib::bam::HeaderView,
    chr_map: &std::collections::HashMap<String, usize>,
    genome: Option<&Genome>,
    snp: Option<&SnpSideChannel>,
    read_tag_table: Option<&ReadTagTable>,
    processor: &ChunkProcessor<'_>,
    args: &QuantCli,
    data: &mut QuantData,
    n_seen: &mut usize,
) -> Result<()> {
    println!("Reading BAM {}", bam_path.display());

    let mut reader = Reader::from_path(bam_path)
        .with_context(|| format!("bam file could not be read: {}", bam_path.display()))?;
    reader.set_threads(3)?; // e.g. n = 4

    validate_bam_header_compatible(expected_header, reader.header())
        .with_context(|| format!("BAM header mismatch for {}", bam_path.display()))?;
    let header = reader.header().clone();
    let job_builder = JobBuilder::new(&header, chr_map, args.cell_tag.0, args.umi_tag.0 )
        .with_genome(genome, !args.no_genome_refine)
        .with_snp_index(snp.map(|s| &s.index))
        .with_read_tag_table( read_tag_table )
        .with_min_mapq(args.min_mapq)
        .read1_only(args.read1_only);

    let mut jobs: Vec<Job> = Vec::with_capacity(CHUNK);

    #[cfg(debug_assertions)]
    let mut i = 0usize;

    for record_result in reader.records() {
        #[cfg(debug_assertions)]
        {
            i += 1;
            if i > 10_000_000 && jobs.is_empty() {
                panic!(
                    "More than 10 million reads failed initial filters in {}\nprocessed total {}\n{}",
                    bam_path.display(),
                    n_seen,
                    data.report
                );
            }
        }

        let rec = record_result.context("BAM read error")?;

        if let Some(job) = job_builder.build(&rec, &mut data.report)? {
            jobs.push(job);
            *n_seen += 1;
        }

        if let Some(max_reads) = args.max_reads
            && *n_seen >= max_reads
        {
            break;
        }

        if jobs.len() >= CHUNK {
            flush_jobs(processor, args.quant_mode, &mut jobs, data)?;
        }
    }

    flush_jobs(processor, args.quant_mode, &mut jobs, data)?;

    Ok(())
}

fn validate_bam_header_compatible(
    expected: &rust_htslib::bam::HeaderView,
    observed: &rust_htslib::bam::HeaderView,
) -> Result<()> {
    if expected.target_count() != observed.target_count() {
        return Err(anyhow!(
            "different number of reference sequences: expected {}, got {}",
            expected.target_count(),
            observed.target_count()
        ));
    }

    for tid in 0..expected.target_count() {
        let exp_name = std::str::from_utf8(expected.tid2name(tid))?;
        let obs_name = std::str::from_utf8(observed.tid2name(tid))?;

        let exp_len = expected.target_len(tid);
        let obs_len = observed.target_len(tid);

        if exp_name != obs_name || exp_len != obs_len {
            return Err(anyhow!(
                "reference mismatch at tid {}: expected {} {:?}, got {} {:?}",
                tid,
                exp_name,
                exp_len,
                obs_name,
                obs_len
            ));
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

fn load_genome(args: &QuantCli) -> Result<Option<Genome>> {
    match &args.genome {
        Some(path) => Genome::from_fasta(path)
            .with_context(|| format!("reading genome FASTA {}", path.display()))
            .map(Some),
        None => Ok(None),
    }
}

fn load_snp_side_channel(
    args: &QuantCli,
    header: &rust_htslib::bam::HeaderView,
) -> Result<Option<SnpSideChannel>> {
    let Some(vcf) = &args.vcf else {
        return Ok(None);
    };

    println!("Loading SNP index");

    let chr_names: Vec<String> = (0..header.target_count())
        .map(|i| std::str::from_utf8(header.tid2name(i)).unwrap().to_string())
        .collect();

    let chr_lengths: Vec<u32> = (0..header.target_count())
        .map(|i| header.target_len(i).unwrap_or(0) as u32)
        .collect();

    SnpSideChannel::from_vcf_path(vcf, chr_names, chr_lengths, args.snp_min_anchor).map(Some)
}

fn print_run_info(genome: Option<&Genome>, snp: Option<&SnpSideChannel>) {
    if genome.is_some() {
        eprintln!(
            "[INFO] genome refinement enabled: reads will be represented as AlignedRead and refined before SNP matching"
        );
    }

    if snp.is_some() {
        eprintln!("[INFO] SNP side-channel enabled: writing additional ref/alt matrices");
    }
}

fn flush_jobs(
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    if jobs.is_empty() {
        return Ok(());
    }

    println!("Processing chunk of size {}", jobs.len());
    data.report.stop_file_io_time();

    processor.process_into(quant_mode, jobs, data)?;

    data.report.stop_single_processor_time();

    jobs.clear();

    Ok(())
}

fn write_log(args: &QuantCli, data: &QuantData) -> Result<()> {
    let log_path = args.outpath.with_extension("log");
    let log_str = format!("{}", data.report);

    let mut file = File::create(&log_path)
        .with_context(|| format!("failed to create log file {}", log_path.display()))?;

    file.write_all(log_str.as_bytes())
        .with_context(|| format!("failed to write log file {}", log_path.display()))?;

    Ok(())
}
