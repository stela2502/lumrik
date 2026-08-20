use std::collections::HashMap;

use anyhow::{bail, Result};
use clap::Parser;

use bam_tide::fastq::FastqPairReader;
use bam_tide::quantification::job::Job;
use bam_tide::results::QuantData;

use nelrune::beacon::run_beacon;
use nelrune::cli::Cli;
use nelrune::fast_features::{
    build_fast_features,
    process_fast_feature_read,
};
use nelrune::fastq::next_parsed_pair;
use nelrune::mapper::{
    drain_mapper,
    remember_tags,
    ReadTags,
};
use nelrune::quant::{
    configure_rayon,
    load_genome,
    load_snp_side_channel,
    CHUNK,
};

use nelrune::progress::RunProgress;
use nelrune::summary::RunSummary;

fn main() -> Result<()> {
    run(Cli::parse())
}

fn run(args: Cli) -> Result<()> {

    let mut progress =
        RunProgress::new();

    progress.stage(
        "loading chemistry"
    );

    validate_args(&args)?;

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

    let mut fastq =
        FastqPairReader::from_paths(
            &args.r1,
            &args.r2,
        )?;

    // --------------------------------------------------------
    // Mapper
    // --------------------------------------------------------

    progress.stage(
        "processing FASTQ"
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

    let Some(first) =
        next_parsed_pair(
            &mut fastq,
            &primer,
        )?
    else {
        bail!(
            "no primer-compatible FASTQ pairs found"
        );
    };

    process_fast_feature_read(
        &first,
        fast_features.as_mut(),
    );

    remember_tags(
        &first,
        &mut pending_tags,
    )?;

    mapper.submit(
        &first.r2,
        None,
    )?;

    let mapper_header =
        mapper
            .header() // blocks!
            .context(
                "waiting for mapper BAM header"
            )?
            .clone();

    let header =
        HeaderView::from_header(
            &mapper_header,
        );

    // --------------------------------------------------------
    // SNP support
    // --------------------------------------------------------

    let snp =
        load_snp_side_channel(
            &args,
            &header,
        )?;

    // --------------------------------------------------------
    // bam-tide processor
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

    loop {
        if args
            .max_reads
            .is_some_and(
                |max| submitted >= max
            )
        {
            break;
        }

        let Some(read) =
            next_parsed_pair(
                &mut fastq,
                &primer,
            )?
        else {
            break;
        };

        // -----------------------------------------------
        // Fast feature side-channel
        // -----------------------------------------------

        process_fast_feature_read(
            &read,
            fast_features.as_mut(),
        );

        // -----------------------------------------------
        // Genomic mapping
        // -----------------------------------------------

        remember_tags(
            &read,
            &mut pending_tags,
        )?;

        mapper.submit(
            &read.r2,
            None,
        )?;

        submitted += 1;

        // -----------------------------------------------
        // Opportunistically drain mapper output
        // -----------------------------------------------
        progress.read_processed();
        
        drain_mapper(
            &mut mapper,
            &mut pending_tags,
            &job_builder,
            &processor,
            args.quant_mode,
            &mut jobs,
            &mut data,
        )?;
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

    progress.stage(
        "selecting cells"
    );    

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

    if let Some(features) = fast_features.as_mut() {
        let feature_index =
            FastTagFeatureIndex::new(
                &features.mapper,
            );
        let raw = std::mem::take(&mut features.data);

        progress.stage(
            "running Beacon"
        );

        let mut beacon =
            sc_beacon::runner::run_from_scdata(
                raw,
                &cells,
                primer.grammar().cell_len(),
                &feature_index,
                &args.beacon.background_config(),
                &args.beacon.fit_config(),
                &args.beacon.call_config(),
            )?;

        let _ = beacon.1.write_sparse(
            &args.outpath.join("fast_features"),
            &feature_index,
        );
        beacon.0.write(
            args.outpath.join("fast_features_stats"),
            &feature_index,
            primer.grammar().cell_len(),
        )?;

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

    write_log(
        &args,
        &data,
        fast_features.as_ref(),
    )?;

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
