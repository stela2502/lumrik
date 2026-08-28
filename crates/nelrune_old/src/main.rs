use std::collections::HashMap;
use std::fs;

use std::net::{
    IpAddr,
    Ipv4Addr,
    SocketAddr,
};

use nelrune::server::spawn_health_server;

use anyhow::{
    anyhow,
    bail,
    Context,
    Result,
};

use clap::Parser;

use rust_htslib::bam::HeaderView;

use mapping_info::MappingInfo;


use gtf_splice_index::{
    MatchOptions,
    SpliceIndex,
};

use sc_primer::PrimerDetector;

use sc_mapper::StreamingMapper;

use fast_tag_mapper::FastTagFeatureIndex;


use bam_tide::fastq::FastqPairReader;

use bam_tide::index::{
    GeneFeatureIndex,
    TranscriptFeatureIndex,
};

use bam_tide::quantification::chunk_processor::ChunkProcessor;

use bam_tide::quantification::cli::QuantMode;

use bam_tide::quantification::job::{
    Job,
    JobBuilder,
};

use bam_tide::quantification::processor_options::ProcessorOptions;

use bam_tide::results::QuantData;

use nelrune::cli::Cli;
use nelrune::mapper::flush_jobs;

use nelrune::fast_features::{
    build_fast_features,
    process_fast_feature_read,
    FastFeatures,
};

use nelrune::fastq::next_parsed_pair;

use nelrune::mapper::{
    consume_mapping_call,
    drain_mapper,
    remember_tags,
    ReadTags,
    validate_reference_compatibility,
};

use nelrune::progress::RunProgress;

use nelrune::quant::{
    configure_rayon,
    load_genome,
    load_snp_side_channel,
    write_quantification,
    CHUNK,
};

use nelrune::summary::RunSummary;
fn main() -> Result<()> {
    run(Cli::parse())
}

fn run(args: Cli) -> Result<()> {

    let mut progress =
        RunProgress::new();

    let health_addr =
        SocketAddr::new(
            IpAddr::V4(
                Ipv4Addr::UNSPECIFIED
            ),
            args.health_port,
        );

    let _health_server =
        spawn_health_server(
            progress.state_handle(),
            health_addr,
        )?;

    let hostname =
        std::env::var("SLURMD_NODENAME")
            .or_else(|_| std::env::var("HOSTNAME"))
            .unwrap_or_else(|_| "localhost".to_string());

    eprintln!(
        "[nelrune] health server listening on http://{}:{}",
        hostname,
        args.health_port
    );

    progress.stage(
        "loading chemistry"
    );

    args.validate()?;

    configure_rayon(args.threads);

    fs::create_dir_all(&args.outpath)
        .with_context(|| {
            format!(
                "creating output directory {}",
                args.outpath.display()
            )
        })?;

    // --------------------------------------------------------
    // Chemistry
    // --------------------------------------------------------

    let primer = args
        .primer
        .detector()
        .map_err(|e| {
            anyhow!("failed to configure sc_primer: {e}")
        })?;

    // --------------------------------------------------------
    // Fast feature mapper
    // --------------------------------------------------------
    progress.stage(
        "loading additional fasta index files"
    );


    let mut fast_features =
        build_fast_features(&args)?;

    // --------------------------------------------------------
    // Quantification resources
    // --------------------------------------------------------

    progress.stage(
        "loading splice index"
    );


    let idx =
        SpliceIndex::load(&args.index)
            .with_context(|| {
                format!(
                    "reading splice index {}",
                    args.index.display()
                )
            })?;

    println!("{idx}");

    let genome =
        load_genome(&args)?;

    if args.vcf.is_some()
        && genome.is_none()
    {
        bail!("--vcf requires --genome");
    }

    // --------------------------------------------------------
    // FASTQ
    // --------------------------------------------------------

    let mut input_iter =
        args.r1
            .iter()
            .zip(args.r2.iter());

    let Some((first_r1, first_r2)) =
        input_iter.next()
    else {
        bail!("no FASTQ inputs supplied");
    };

    progress.stage(
        format!(
            "processing {}",
            first_r1.display(),
        )
    );
    progress.input_file( first_r1.to_string_lossy(), );

    let mut fastq =
        FastqPairReader::from_paths(
            first_r1,
            first_r2,
        )?;

    // --------------------------------------------------------
    // Mapper
    // --------------------------------------------------------

    progress.stage(
        "starting mapper"
    );

    let mut mapper =
        args.mapper
            .from_cli()
            .context(
                "starting streaming genomic mapper"
            )?;

    // --------------------------------------------------------
    // Quantification output
    // --------------------------------------------------------

    let mut data =
        QuantData::new();

    data.report =
        MappingInfo::new(
            None,
            args.min_mapq as f32,
            args.max_reads
                .unwrap_or(usize::MAX),
        );

    data.report.start_counter();

    // --------------------------------------------------------
    // Primer metadata waiting for asynchronous mappings
    // --------------------------------------------------------

    let mut pending_tags:
        HashMap<String, ReadTags> =
        HashMap::new();

    // --------------------------------------------------------
    // Bootstrap mapper
    // --------------------------------------------------------
    //
    // Submit one usable read before waiting for the mapper header.
    //
    // Some external mappers only produce stdout/header once input
    // has begun arriving.

    progress.stage(
        "feeding initial FASTQ to mapper"
    );

    let mut submitted = 0usize;
    const BOOTSTRAP_READS: usize = 10_000;

    while !mapper.header_loaded() {
        let mut exhausted = false;

        for _ in 0..BOOTSTRAP_READS {
            let Some(read) =
                next_parsed_pair(
                    &mut fastq,
                    &primer,
                )?
            else {
                exhausted = true;
                break;
            };

            process_fast_feature_read(
                &read,
                fast_features.as_mut(),
            );

            remember_tags(
                &read,
                &mut pending_tags,
            )?;

            mapper.submit(
                &read.r2,
                None,
            )?;

            submitted += 1;
            progress.read_processed();

            if args
                .max_reads
                .is_some_and(|max| submitted >= max)
            {
                exhausted = true;
                break;
            }
        }

        if exhausted {
            break;
        }
    }

    progress.stage(
        "waiting for mapper header"
    );

    let mapper_header =
        mapper
            .header()
            .context(
                "waiting for mapper BAM header"
            )?
            .clone();

    // Synchronization point:
    // wait until the mapper has started, loaded its reference,
    // and emitted its SAM/BAM header.

    let header =
        HeaderView::from_header(
            &mapper_header,
        );
    progress.stage(
        "mapper ready"
    );

    progress.stage(
        "checking reference compatibility"
    );
    validate_reference_compatibility(&header, &idx )?;


    // --------------------------------------------------------
    // SNP support
    // --------------------------------------------------------

    let snp =
        load_snp_side_channel(
            &args,
            &header,
        )?;

    // --------------------------------------------------------
    // bam-tide set up processors and collectors
    // --------------------------------------------------------

    let match_opts =
        MatchOptions {
            require_strand:
                args.require_strand,

            require_exact_junction_chain:
                args.require_exact_junction_chain,

            max_5p_overhang_bp:
                args.max_5p_overhang_bp,

            max_3p_overhang_bp:
                args.max_3p_overhang_bp,

            allowed_intronic_gap_size:
                args.allowed_intronic_gap_size,
        };

    let processor_options =
        ProcessorOptions {
            min_mapq:
                args.min_mapq,

            read1_only:
                false,

            require_strand:
                args.require_strand,

            quant_mode:
                args.quant_mode,

            ..ProcessorOptions::default()
        };

    let processor =
        ChunkProcessor::new(
            &idx,
            snp.as_ref(),
            match_opts,
            processor_options,
        );

    let job_builder =
        JobBuilder::new(
            &header,
            idx.chr_map(),
            *b"CB",
            *b"UB",
        )
        .with_genome(
            genome.as_ref(),
            !args.no_genome_refine,
        )
        .with_snp_index(
            snp.as_ref()
                .map(|s| &s.index),
        )
        .with_min_mapq(
            args.min_mapq,
        );

    data.report.stop_file_io_time();

    // --------------------------------------------------------
    // Main streaming loop
    // --------------------------------------------------------

    let mut jobs =
        Vec::<Job>::with_capacity(CHUNK);

    let mut submitted =
        1usize;

    process_fastq_reader(
        &mut fastq,
        &primer,
        &mut fast_features,
        &mut mapper,
        &mut pending_tags,
        &job_builder,
        &processor,
        args.quant_mode,
        &mut jobs,
        &mut data,
        &mut progress,
        &mut submitted,
        args.max_reads,
    )?;

    for (r1, r2) in input_iter {
        progress.stage(
            format!(
                "processing {}",
                r1.display(),
            )
        );

        progress.input_file(
            r1.to_string_lossy(),
        );

        let mut fastq =
            FastqPairReader::from_paths(
                r1,
                r2,
            )?;

        if !process_fastq_reader(
            &mut fastq,
            &primer,
            &mut fast_features,
            &mut mapper,
            &mut pending_tags,
            &job_builder,
            &processor,
            args.quant_mode,
            &mut jobs,
            &mut data,
            &mut progress,
            &mut submitted,
            args.max_reads,
        )? {
            break;
        }
    }

    // --------------------------------------------------------
    // Finish mapper
    // --------------------------------------------------------

    progress.stage(
        "finishing mapper"
    );

    for call in mapper.finish()? {
        consume_mapping_call(
            call,
            &mut pending_tags,
            &job_builder,
            &processor,
            args.quant_mode,
            &mut jobs,
            &mut data,
        )?;
    }

    flush_jobs(
        &processor,
        args.quant_mode,
        &mut jobs,
        &mut data,
    )?;

    if !pending_tags.is_empty() {
        eprintln!(
            "[WARN] {} submitted reads never produced mapper output",
            pending_tags.len()
        );
    }

    progress.clear_input_file();


    // --------------------------------------------------------
    // Output
    // --------------------------------------------------------
    let cells = match args.quant_mode {
        QuantMode::Gene => {
            let gene_index =
                GeneFeatureIndex::new(&idx);

            data.finalize_for_export(
                args.min_cell_counts,
                &gene_index,
                snp.as_ref().map(|s| &s.index),
            )
        }

        QuantMode::Transcript => {
            let transcript_index =
                TranscriptFeatureIndex::new(&idx);

            data.finalize_for_export(
                args.min_cell_counts,
                &transcript_index,
                snp.as_ref().map(|s| &s.index),
            )
        }
    };

    eprintln!(
        "[nelrune] detected {} cells",
        cells.len(),
    );

    let mut beacon = None;

    if let Some(features) = fast_features.as_mut() {
        let feature_index =
            FastTagFeatureIndex::new(
                &features.mapper,
            );

        let raw =
            std::mem::take(
                &mut features.data,
            );

        progress.stage(
            "running Beacon"
        );

        let (
            beacon_result,
            mut filtered_features,
        ) =
            sc_beacon::runner::run_from_scdata(
                raw,
                &cells,
                primer.grammar().cell_len(),
                &feature_index,
                &args.beacon.background_config(),
                &args.beacon.fit_config(),
                &args.beacon.call_config(),
            )?;

        /*
         * The split already selected the correct cells,
         * but finalize_for_cells() rebuilds export metadata.
         */
        filtered_features.finalize_for_cells(
            &cells,
            &feature_index,
        );

        filtered_features
            .write_sparse(
                &args.outpath.join("fast_features"),
                &feature_index,
            )
            .map_err(anyhow::Error::msg)?;

        beacon_result.write(
            args.outpath.join("fast_features_stats"),
            &feature_index,
            primer.grammar().cell_len(),
        )?;

        /*
         * Keep it for the final Nelrune summary.
         */
        beacon = Some(beacon_result);

        /*
         * Optionally retain the filtered matrix in FastFeatures too,
         * if later code expects it there.
         */
        features.data = filtered_features;
    }

    progress.stage(
        "writing output"
    );

    write_quantification(
        &args,
        &idx,
        snp.as_ref(),
        &mut data,
    )?;


    data.report.stop_file_io_time();

    println!("{}", data.report);

    if let Some(features) =
        fast_features.as_ref()
    {
        println!(
            "Fast feature mapping:\n{}",
            features.report
        );
    }


    let summary =
        RunSummary::from_run(
            &progress,
            &data,
            cells.len(),
            fast_features.as_ref(),
            beacon.as_ref(),
        );

    println!();
    println!("{summary}");

    summary.write(
        args.outpath.join(
            "nelrune.pretty.log"
        ),
    )?;

    Ok(())
}


pub fn process_fastq_reader(
    fastq: &mut FastqPairReader,
    primer: &PrimerDetector,
    fast_features: &mut Option<FastFeatures>,
    mapper: &mut StreamingMapper,
    pending_tags: &mut HashMap<String, ReadTags>,
    job_builder: &JobBuilder<'_>,
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
    progress: &mut RunProgress,
    submitted: &mut usize,
    max_reads: Option<usize>,
) -> Result<bool> {
    loop {
        if max_reads.is_some_and(
            |max| *submitted >= max
        ) {
            return Ok(false);
        }

        let Some(read) =
            next_parsed_pair(
                fastq,
                primer,
            )?
        else {
            return Ok(true);
        };

        process_fast_feature_read(
            &read,
            fast_features.as_mut(),
        );

        remember_tags(
            &read,
            pending_tags,
        )?;

        mapper.submit(
            &read.r2,
            None,
        )?;

        *submitted += 1;

        progress.read_processed();

        drain_mapper(
            mapper,
            pending_tags,
            job_builder,
            processor,
            quant_mode,
            jobs,
            data,
        )?;
    }
}