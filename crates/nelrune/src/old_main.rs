//! Nelrune
//!
//! Thin orchestration binary:
//!
//! FASTQ
//!   -> sc_primer
//!   -> FastTagMapper side-channel
//!   -> sc_mapper
//!   -> bam_tide quantification
//!   -> scdata output
//!
//! Current initial assumptions:
//!
//! - R1 contains the single-cell chemistry / cell barcode / UMI.
//! - R2 contains the biological insert.
//! - R2 is searched against all FastTagMapper features.
//! - R2 is also submitted to the genomic mapper.
//!
//! These assumptions can later become chemistry/collection-driven without
//! changing the downstream libraries.

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::{ PathBuf};
use std::str::FromStr;

use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;

use rust_htslib::bam;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{HeaderView, Record};

use int_to_str::IntToStr;

use bam_tide::fastq::{
    FastqPairReader,
    FastqRecord,
};
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
use bam_tide::quantification::snp::SnpSideChannel;
use bam_tide::results::QuantData;

use fast_tag_mapper::{
    BuiltinTagSet,
    FastTagFeatureIndex,
    FastTagMapper,
};

use gtf_splice_index::{
    MatchOptions,
    SpliceIndex,
};

use mapping_info::MappingInfo;

use sc_beacon::cli::GuideModelCli;
use sc_beacon::runner::run_from_scdata;

use sc_mapper::{
    MappingCall,
    StreamingMapper,
    StreamingMapperCli,
};

use sc_primer::{
    PrimerCli,
    PrimerDetector,
};

use scdata::{
    GeneUmiHash,
    MatrixValueType,
    Scdata,
};

use snp_index::Genome;


const CHUNK: usize = 2_000_000;


// ============================================================
// Fast feature source
// ============================================================

#[derive(Debug, Clone)]
enum FastFeatureSource {
    BdSampleHuman,
    BdSampleMouse,
    Fasta(PathBuf),
}

impl FromStr for FastFeatureSource {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "bd_sample_human" => {
                Ok(Self::BdSampleHuman)
            }

            "bd_sample_mouse" => {
                Ok(Self::BdSampleMouse)
            }

            value if value.is_empty() => {
                bail!("empty --fast-features entry")
            }

            path => {
                Ok(Self::Fasta(PathBuf::from(path)))
            }
        }
    }
}


// ============================================================
// CLI
// ============================================================

#[derive(Parser, Debug)]
#[command(
    name = "nelrune",
    about = "Single-cell FASTQ processing using the Lumrik libraries"
)]
struct Cli {
    /// Barcode / UMI read.
    #[arg(long)]
    r1: PathBuf,

    /// Biological insert read.
    #[arg(long)]
    r2: PathBuf,

    // --------------------------------------------------------
    // Chemistry
    // --------------------------------------------------------

    #[command(flatten)]
    primer: PrimerCli,

    // --------------------------------------------------------
    // Genomic mapping
    // --------------------------------------------------------

    #[command(flatten)]
    mapper: StreamingMapperCli,

    /// Gene stats model
    #[command(flatten)]
    beacon: GuideModelCli,

    /// bam-tide splice index.
    #[arg(long, short)]
    index: PathBuf,


    // --------------------------------------------------------
    // Fast feature mapping
    // --------------------------------------------------------

    /// FastTagMapper references.
    ///
    /// Built-ins:
    ///
    ///   bd_sample_human
    ///   bd_sample_mouse
    ///
    /// Anything else is interpreted as a FASTA path.
    ///
    /// Example:
    ///
    ///   --fast-features \
    ///       'bd_sample_human;guides.fa;hto.fa'
    #[arg(
        long,
        value_delimiter = ';'
    )]
    fast_features: Vec<FastFeatureSource>,

    /// Minimum supporting 8-mer hits for fast feature mapping.
    #[arg(long, default_value_t = 4)]
    fast_feature_min_hits: u32,

    // --------------------------------------------------------
    // Quantification
    // --------------------------------------------------------

    #[arg(long, default_value_t = 0)]
    min_mapq: u8,

    #[arg(
        long,
        value_enum,
        default_value_t = QuantMode::Gene
    )]
    quant_mode: QuantMode,

    #[arg(long, default_value_t = 400)]
    min_cell_counts: usize,

    #[arg(long)]
    max_reads: Option<usize>,

    #[arg(long, default_value_t = 0)]
    threads: usize,

    // --------------------------------------------------------
    // Optional genome / SNP support
    // --------------------------------------------------------

    #[arg(long)]
    genome: Option<PathBuf>,

    #[arg(long)]
    vcf: Option<PathBuf>,

    #[arg(long, default_value_t = 20)]
    snp_min_anchor: u8,

    #[arg(long, default_value_t = false)]
    no_genome_refine: bool,

    // --------------------------------------------------------
    // Splice matching
    // --------------------------------------------------------

    #[arg(long, default_value_t = false)]
    require_strand: bool,

    #[arg(long, default_value_t = false)]
    require_exact_junction_chain: bool,

    #[arg(long, default_value_t = 100)]
    max_5p_overhang_bp: u32,

    #[arg(long, default_value_t = 100)]
    max_3p_overhang_bp: u32,

    #[arg(long, default_value_t = 5)]
    allowed_intronic_gap_size: u32,

    // --------------------------------------------------------
    // Output
    // --------------------------------------------------------

    #[arg(long, short)]
    outpath: PathBuf,
}


// ============================================================
// Per-read metadata
// ============================================================

#[derive(Debug)]
struct ReadTags {
    cell_seq: Vec<u8>,
    umi_seq: Vec<u8>,
}


// ============================================================
// Fast feature matrix
// ============================================================

struct FastFeatures {
    mapper: FastTagMapper,
    data: Scdata,
    report: MappingInfo,
}

impl FastFeatures {
    fn new(
        mapper: FastTagMapper,
        threads: usize,
    ) -> Self {
        Self {
            mapper,
            data: Scdata::new(
                threads.max(1),
                MatrixValueType::Real,
            ),
            report: MappingInfo::new(
                None,
                0.0,
                usize::MAX,
            ),
        }
    }
}


// ============================================================
// main
// ============================================================

fn main() -> Result<()> {
    run(Cli::parse())
}


// ============================================================
// Runner
// ============================================================

fn run(args: Cli) -> Result<()> {
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

    let mut fast_features =
        build_fast_features(&args)?;

    // --------------------------------------------------------
    // Quantification resources
    // --------------------------------------------------------

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

    if let Some(features) = fast_features.as_mut() {
        let feature_index =
            FastTagFeatureIndex::new(
                &features.mapper,
            );
        let raw = std::mem::take(&mut features.data);

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

    Ok(())
}


// ============================================================
// Parsed FASTQ pair
// ============================================================

struct ParsedPair {
    read_id: String,

    cell_seq: Vec<u8>,
    umi_seq: Vec<u8>,

    cell_id: u64,
    umi_id: u64,

    r2: FastqRecord,
}


// ============================================================
// FASTQ + sc_primer
// ============================================================

fn next_parsed_pair(
    reader: &mut FastqPairReader,
    primer: &PrimerDetector,
) -> Result<Option<ParsedPair>> {
    loop {
        let Some((r1, r2)) =
            reader.next_pair()?
        else {
            return Ok(None);
        };

        let r1_id =
            r1.clean_id();

        let r2_id =
            r2.clean_id();

        if r1_id != r2_id {
            bail!(
                "FASTQ pair ID mismatch: R1='{r1_id}', R2='{r2_id}'"
            );
        }

        let Some(hit) =
            primer
                .detect_first(
                    &r1.seq,
                    &r1.qual,
                )
                .map_err(|e| {
                    anyhow!(
                        "primer detection failed for read '{}': {e}",
                        r1_id
                    )
                })?
        else {
            continue;
        };

        let cell =
            hit.get_cell(
                &r1.seq,
                &r1.qual,
            )
            .map_err(|e| {
                anyhow!(
                    "failed to extract cell barcode from '{}': {e}",
                    r1_id
                )
            })?;

        let umi =
            hit.get_umi(
                &r1.seq,
                &r1.qual,
            )
            .map_err(|e| {
                anyhow!(
                    "failed to extract UMI from '{}': {e}",
                    r1_id
                )
            })?;

        let cell_id =
            dna_to_u64(
                &cell.seq,
                "cell barcode",
            )?;

        let umi_id =
            dna_to_u64(
                &umi.seq,
                "UMI",
            )?;

        return Ok(
            Some(
                ParsedPair {
                    read_id: r2_id,

                    cell_seq:
                        cell.seq,

                    umi_seq:
                        umi.seq,

                    cell_id,
                    umi_id,

                    r2,
                }
            )
        );
    }
}


// ============================================================
// Fast features
// ============================================================

fn build_fast_features(
    args: &Cli,
) -> Result<Option<FastFeatures>> {
    if args.fast_features.is_empty() {
        return Ok(None);
    }

    let mut mapper = FastTagMapper::new();

    for source in
        &args.fast_features
    {
        match source {
            FastFeatureSource::BdSampleHuman => {
                let n =
                    mapper.add_builtin(
                        BuiltinTagSet::Human,
                    );

                println!(
                    "Loaded {n} built-in BD human sample tags"
                );
            }

            FastFeatureSource::BdSampleMouse => {
                let n =
                    mapper.add_builtin(
                        BuiltinTagSet::Mouse,
                    );

                println!(
                    "Loaded {n} built-in BD mouse sample tags"
                );
            }

            FastFeatureSource::Fasta(path) => {
                let n =
                    mapper
                        .load_fasta(path)
                        .with_context(|| {
                            format!(
                                "loading fast feature FASTA {}",
                                path.display()
                            )
                        })?;

                println!(
                    "Loaded {n} fast features from {}",
                    path.display()
                );
            }
        }
    }

    mapper =
        mapper.with_min_hits(
            args.fast_feature_min_hits,
        );

    println!(
        "Fast feature mapper contains {} features",
        mapper.feature_count()
    );

    Ok(Some(
        FastFeatures::new(
            mapper,
            args.threads,
        )
    ))
}


fn process_fast_feature_read(
    read: &ParsedPair,
    features: Option<&mut FastFeatures>,
) {
    let Some(features) =
        features
    else {
        return;
    };

    /*
     * Initial Nelrune model:
     *
     * R2 is searched against every requested short-feature
     * reference.
     *
     * No distinction exists here between guide / HTO /
     * sample tag / feature barcode / extra synthetic gene.
     */

    let Some(feature_id) =
        features
            .mapper
            .map_feature_id(
                &read.r2.seq,
                &mut features.report,
            )
    else {
        return;
    };

    features.data.try_insert(
        &read.cell_id,
        GeneUmiHash(
            feature_id,
            read.umi_id,
        ),
        1.0,
        &mut features.report,
    );
}


// ============================================================
// Pending CB / UMI metadata
// ============================================================

fn remember_tags(
    read: &ParsedPair,
    pending: &mut HashMap<String, ReadTags>,
) -> Result<()> {
    if pending
        .insert(
            read.read_id.clone(),
            ReadTags {
                cell_seq:
                    read.cell_seq.clone(),

                umi_seq:
                    read.umi_seq.clone(),
            },
        )
        .is_some()
    {
        bail!(
            "duplicate outstanding read id '{}'",
            read.read_id
        );
    }

    Ok(())
}


// ============================================================
// sc_mapper -> bam_tide
// ============================================================

fn drain_mapper(
    mapper: &mut StreamingMapper,
    pending: &mut HashMap<String, ReadTags>,
    job_builder: &JobBuilder<'_>,
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    while let Some(call) =
        mapper.try_next()?
    {
        consume_mapping_call(
            call,
            pending,
            job_builder,
            processor,
            quant_mode,
            jobs,
            data,
        )?;
    }

    Ok(())
}


fn consume_mapping_call(
    call: MappingCall,
    pending: &mut HashMap<String, ReadTags>,
    job_builder: &JobBuilder<'_>,
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    let tags =
        pending
            .remove(&call.read_id)
            .ok_or_else(|| {
                anyhow!(
                    "mapper returned read '{}' without matching primer metadata",
                    call.read_id
                )
            })?;

    let cell =
        std::str::from_utf8(
            &tags.cell_seq,
        )
        .context(
            "cell barcode is not valid ASCII"
        )?;

    let umi =
        std::str::from_utf8(
            &tags.umi_seq,
        )
        .context(
            "UMI is not valid ASCII"
        )?;

    /*
     * sc_mapper may return multiple alignments.
     *
     * For now these are passed through the existing JobBuilder
     * semantics. JobBuilder currently ignores secondary /
     * supplementary records.
     *
     * Proper multimapper-aware quantification belongs in
     * bam_tide, not here.
     */

    for mapper_record in
        call.records.records
    {
        let mut rec =
            mapper_record.into_inner();

        set_aux_string(
            &mut rec,
            b"CB",
            cell,
        )?;

        set_aux_string(
            &mut rec,
            b"UB",
            umi,
        )?;

        if let Some(job) =
            job_builder.build(
                &rec,
                &mut data.report,
            )?
        {
            jobs.push(job);
        }
    }

    if jobs.len() >= CHUNK {
        flush_jobs(
            processor,
            quant_mode,
            jobs,
            data,
        )?;
    }

    Ok(())
}


fn set_aux_string(
    record: &mut Record,
    tag: &[u8; 2],
    value: &str,
) -> Result<()> {
    /*
     * Mapper output should not already carry Nelrune's
     * barcode tags, but remove first so this is deterministic.
     */

    let _ =
        record.remove_aux(tag);

    record
        .push_aux(
            tag,
            Aux::String(value),
        )
        .with_context(|| {
            format!(
                "adding BAM tag {}",
                String::from_utf8_lossy(tag)
            )
        })?;

    Ok(())
}


// ============================================================
// bam-tide chunk handling
// ============================================================

fn flush_jobs(
    processor: &ChunkProcessor<'_>,
    quant_mode: QuantMode,
    jobs: &mut Vec<Job>,
    data: &mut QuantData,
) -> Result<()> {
    if jobs.is_empty() {
        return Ok(());
    }

    println!(
        "Processing {} mapping jobs",
        jobs.len()
    );

    data.report.stop_file_io_time();

    processor.process_into(
        quant_mode,
        jobs,
        data,
    )?;

    jobs.clear();

    Ok(())
}


// ============================================================
// Quantification output
// ============================================================

fn write_quantification(
    args: &Cli,
    idx: &SpliceIndex,
    snp: Option<&SnpSideChannel>,
    data: &mut QuantData,
) -> Result<()> {
    println!(
        "Writing genomic quantification"
    );

    match args.quant_mode {
        QuantMode::Gene => {
            let features =
                GeneFeatureIndex::new(idx);

            data.write(
                &args.outpath,
                args.min_cell_counts,
                &features,
                snp.map(|s| &s.index),
            )
            .map_err(anyhow::Error::msg)
            .context(
                "writing gene quantification"
            )?;
        }

        QuantMode::Transcript => {
            let features =
                TranscriptFeatureIndex::new(idx);

            data.write(
                &args.outpath,
                args.min_cell_counts,
                &features,
                snp.map(|s| &s.index),
            )
            .map_err(anyhow::Error::msg)
            .context(
                "writing transcript quantification"
            )?;
        }
    }

    Ok(())
}


// ============================================================
// Fast feature output
// ============================================================

fn write_fast_features(
    args: &Cli,
    features: &mut FastFeatures,
    data: &QuantData,
) -> Result<()> {
    let cell_set:
        HashSet<u64> =
        data
            .gene
            .export_cell_ids()
            .iter()
            .copied()
            .collect();

    let index =
        FastTagFeatureIndex::new(
            &features.mapper,
        );

    features
        .data
        .finalize_for_cells(
            &cell_set,
            &index,
        );

    let out =
        args.outpath.join(
            "fast_features",
        );

    fs::create_dir_all(&out)
        .with_context(|| {
            format!(
                "creating {}",
                out.display()
            )
        })?;

    features
        .data
        .write_sparse(
            &out,
            &index,
        )
        .map_err(anyhow::Error::msg)
        .context(
            "writing fast feature matrix"
        )?;

    Ok(())
}


// ============================================================
// Genome / SNP
// ============================================================

fn load_genome(
    args: &Cli,
) -> Result<Option<Genome>> {
    match &args.genome {
        Some(path) => {
            Genome::from_fasta(path)
                .with_context(|| {
                    format!(
                        "reading genome FASTA {}",
                        path.display()
                    )
                })
                .map(Some)
        }

        None => Ok(None),
    }
}


fn load_snp_side_channel(
    args: &Cli,
    header: &HeaderView,
) -> Result<Option<SnpSideChannel>> {
    let Some(vcf) =
        &args.vcf
    else {
        return Ok(None);
    };

    println!("Loading SNP index");

    let chr_names:
        Vec<String> =
        (0..header.target_count())
            .map(|tid| {
                String::from_utf8_lossy(
                    header.tid2name(tid),
                )
                .to_string()
            })
            .collect();

    let chr_lengths:
        Vec<u32> =
        (0..header.target_count())
            .map(|tid| {
                header
                    .target_len(tid)
                    .unwrap_or(0)
                    as u32
            })
            .collect();

    SnpSideChannel::from_vcf_path(
        vcf,
        chr_names,
        chr_lengths,
        args.snp_min_anchor,
    )
    .map(Some)
}


// ============================================================
// Misc
// ============================================================

fn dna_to_u64(
    seq: &[u8],
    label: &str,
) -> Result<u64> {
    IntToStr::try_new(seq)
        .map(|x| x.into_u64())
        .map_err(|e| anyhow!("{label}: {e}"))
}


fn configure_rayon(
    threads: usize,
) {
    if threads == 0 {
        return;
    }

    let _ =
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
}


fn validate_args(
    args: &Cli,
) -> Result<()> {
    /*
     * At the moment R1 is chemistry and R2 is biological
     * insert, therefore mapper-level paired mode does not
     * describe the input correctly yet.
     */

    if args.mapper.mapper_paired {
        bail!(
            "Nelrune currently maps R2 as a single biological insert; \
             --mapper-paired is not supported yet"
        );
    }

    Ok(())
}


fn write_log(
    args: &Cli,
    data: &QuantData,
    fast_features: Option<&FastFeatures>,
) -> Result<()> {
    let path =
        args.outpath.join(
            "nelrune.log"
        );

    let mut file =
        File::create(&path)
            .with_context(|| {
                format!(
                    "creating {}",
                    path.display()
                )
            })?;

    writeln!(
        file,
        "Genomic mapping / quantification"
    )?;

    writeln!(
        file,
        "=============================="
    )?;

    writeln!(
        file,
        "{}",
        data.report
    )?;

    if let Some(features) =
        fast_features
    {
        writeln!(file)?;
        writeln!(
            file,
            "Fast feature mapping"
        )?;

        writeln!(
            file,
            "===================="
        )?;

        writeln!(
            file,
            "{}",
            features.report
        )?;
    }

    Ok(())
}