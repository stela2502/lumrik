use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};

use clap::{Args, Parser, Subcommand};

use gtf_splice_index::{IdNameKeys, SpliceIndex}; // <-- adjust crate path/module as needed
use bam_tide::core::fasta::FastaRecord;
use snp_index::Genome;

use flate2::write::GzEncoder;
use flate2::Compression;

/// Build, inspect, or serialize a splice index.
#[derive(Parser, Debug)]
#[command(name = "splice-index")]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build an index from a GTF/GFF annotation and write it to disk
    Build(BuildArgs),

    /// Load an index from disk and print summary stats
    Stats(StatsArgs),

/*
    /// Query transcripts covering a genomic position
    Query(QueryArgs),

*/

    /// Build transcriptome FASTA from genomic annotation + genome FASTA
    Transcriptome(TranscriptomeArgs),


}

#[derive(Args, Debug)]
struct StatsArgs {
    /// Serialized index file
    #[arg(long, short)]
    index: PathBuf,
}


#[derive(Args, Debug)]
struct TranscriptomeArgs {
    /// Genomic GTF/GFF or binary gtf-splice-index used to define transcripts and exon structure.
    ///
    /// The annotation remains in genomic coordinates. This command only emits
    /// transcript sequences extracted from the genome FASTA.
    #[arg(long, short = 'a')]
    index: PathBuf,

    /// Genome FASTA used as sequence source.
    ///
    /// Plain `.fa/.fasta` and gzip-compressed `.fa.gz/.fasta.gz` are supported
    /// if the underlying Genome loader supports them.
    #[arg(long, short = 'g')]
    genome: PathBuf,

    /// Output transcriptome FASTA.
    ///
    /// If the path ends in `.gz`, gzip-compressed FASTA is written.
    #[arg(long, short = 'o')]
    out: PathBuf,

    /// Optional genomic region restriction, e.g. chr17:7668402-7687550.
    ///
    /// Coordinates are 1-based inclusive.
    #[arg(long, short = 'r')]
    region: Option<String>,

    /// Bin width used when building the temporary splice index from annotation.
    #[arg(long, default_value_t = 1_000_000)]
    bin_width: u32,

    /// Number of bases per FASTA sequence line.
    #[arg(long, default_value_t = 80)]
    line_width: usize,
}

/*
#[derive(Args, Debug)]
struct QueryArgs {
    /// Build the splice index directly from this GTF/GFF annotation file.
    ///
    /// Use this for quick checks without first creating a `.gtf.dat` index.
    /// Mutually exclusive with `--index`.
    #[arg(long, short = 'a', conflicts_with = "index", required_unless_present = "index")]
    annotation: Option<PathBuf>,

    /// Load a prebuilt binary splice index.
    ///
    /// This is faster than rebuilding from annotation and is the preferred
    /// mode inside pipelines.
    /// Mutually exclusive with `--annotation`.
    #[arg(long, short = 'i', conflicts_with = "annotation", required_unless_present = "annotation")]
    index: Option<PathBuf>,

    /// Bin width used when building from `--annotation`.
    ///
    /// Ignored when loading an existing `--index`.
    #[arg(long, default_value_t = 1_000_000)]
    bin_width: u32,

    /// Query region in genomic coordinates.
    ///
    /// Accepted forms:
    ///   chr17:7674220
    ///   chr17:7674220-7674300
    ///   chr17:7,674,220-7,674,300
    ///
    /// Coordinates are 1-based and inclusive, like VCF/GTF-style positions.
    #[arg(long, short = 'r', conflicts_with_all = ["chr", "pos"])]
    region: Option<String>,

    /// Chromosome/contig name for point queries.
    ///
    /// Example: `chr17`, `17`, `chrM`, `MT`.
    /// Use together with `--pos`.
    #[arg(long, short = 'c', required_unless_present = "region")]
    chr: Option<String>,

    /// 1-based genomic position for point queries.
    ///
    /// Example: `7674220`.
    /// Use together with `--chr`.
    #[arg(long, short = 'p', required_unless_present = "region")]
    pos: Option<u32>,
}
*/

#[derive(Args, Debug)]
struct BuildArgs {
    /// Input annotation file (.gtf/.gff/.gff3)
    #[arg(long, short)]
    annotation: PathBuf,

    /// Bin width in base pairs
    #[arg(long, short, default_value_t = 1_000_000)]
    bin_width: u32,

    /// Output serialized index file
    #[arg(long, short)]
    index: PathBuf,

    // -------------------------
    // Attribute key options
    // -------------------------
    /// Attribute keys to use for gene ID (repeatable).
    /// Default (GTF-safe): gene_id
    #[arg(
        long = "gene-id-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["gene_id".to_string()]
    )]
    gene_id_keys: Vec<String>,

    /// Attribute keys to use for gene name (repeatable).
    /// Default (GTF-safe): gene_name
    #[arg(
        long = "gene-name-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["gene_name".to_string()]
    )]
    gene_name_keys: Vec<String>,

    /// Attribute keys to use for transcript ID (repeatable).
    /// Default (GTF-safe): transcript_id
    #[arg(
        long = "transcript-id-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["transcript_id".to_string()]
    )]
    transcript_id_keys: Vec<String>,

    /// Attribute keys to use for transcript name (repeatable).
    /// Default (GTF-safe): transcript_name
    #[arg(
        long = "transcript-name-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["transcript_name".to_string()]
    )]
    transcript_name_keys: Vec<String>,

    /// GFF3 exon->transcript linkage keys (repeatable).
    /// Default: Parent
    #[arg(
        long = "parent-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["Parent".to_string()]
    )]
    parent_keys: Vec<String>,

    /// Feature types that count as exon blocks (repeatable).
    /// Default: exon
    #[arg(
        long = "exon-feature-type",
        value_name = "TYPE",
        num_args = 1..,
        default_values_t = vec!["exon".to_string()]
    )]
    exon_feature_types: Vec<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.cmd {
        Command::Build(args) => {
            let keys = IdNameKeys {
                gene_id_keys: args.gene_id_keys,
                gene_name_keys: args.gene_name_keys,
                transcript_id_keys: args.transcript_id_keys,
                transcript_name_keys: args.transcript_name_keys,
                parent_keys: args.parent_keys,
                exon_feature_types: args.exon_feature_types,
            };

            let idx = SpliceIndex::from_path(&args.annotation, args.bin_width, keys)
                .with_context(|| format!("building index from {}", args.annotation.display()))?;

            println!("{idx}");

            idx.save(&args.index)
                .with_context(|| format!("writing index to {}", args.index.display()))?;

            eprintln!("Index written to {}", args.index.display());
        }

        Command::Transcriptome(args) => {
            build_transcriptome(&args)?
        }

        Command::Stats(args) => {
            let idx = SpliceIndex::load(&args.index)
                .with_context(|| format!("reading index {}", args.index.display()))?;
            println!("{idx}");
        }
    }

    Ok(())
}


fn build_transcriptome(args: &TranscriptomeArgs) -> Result<()> {
    let index = if args
        .index
        .extension()
        .map(|e| e == "dat")
        .unwrap_or(false)
    {
        SpliceIndex::load(&args.index)
            .with_context(|| format!("reading splice index {}", args.index.display()))?
    } else {
        SpliceIndex::from_path(
            &args.index,
            args.bin_width,
            IdNameKeys::default(),
        )
        .with_context(|| format!("building index from {}", args.index.display()))?
    };

    let genome = Genome::from_fasta(&args.genome)
        .with_context(|| {
            format!(
                "failed to load genome FASTA {}",
                args.genome.display()
            )
        })?;

    let out = std::fs::File::create(&args.out)
        .with_context(|| {
            format!(
                "failed to create output FASTA {}",
                args.out.display()
            )
        })?;

    let writer: Box<dyn std::io::Write> = if args
        .out
        .extension()
        .map(|e| e == "gz")
        .unwrap_or(false)
    {
        Box::new(GzEncoder::new(
            out,
            Compression::default(),
        ))
    } else {
        Box::new(out)
    };

    let mut writer = std::io::BufWriter::new(writer);

    let mut n_written = 0usize;

    let region: Option<(String, u32, u32)> = args
        .region
        .as_deref()
        .map(parse_region)
        .transpose()?;


    for tx in &index.transcripts {
        
        if let Some((region_chr, region_start0, region_end0)) = &region {
            let region_chr_id = index
                .chr_id(region_chr)
                .with_context(|| {
                    format!("region chromosome not found in index: {region_chr}")
                })?;

            if tx.chr_id != region_chr_id {
                continue;
            }

            let Some((tx_start, tx_end)) = tx.span() else {
                continue;
            };

            if !(tx_start < *region_end0 && *region_start0 < tx_end) {
                continue;
            }
        }

        let chr_name = &index.chr_names[tx.chr_id];

        /*
        let Some(_tx_name) = tx.primary_name() else {
            continue;
        };

        let genome_chr_id = genome
            .chr_id(chr_name)
            .with_context(|| {
                format!("chromosome missing from genome FASTA: {chr_name}")
            })?;
        */

        let record = FastaRecord::from_transcript(
            &genome,
            tx,
            chr_name,
        )?;

        record.write(&mut writer, args.line_width)?;

        n_written += 1;
    }

    eprintln!(
        "Wrote {} transcript sequences to {}",
        n_written,
        args.out.display()
    );
    Ok(())
}


fn parse_region(region: &str) -> Result<(String, u32, u32)> {
    let (chr, rest) = region
        .split_once(':')
        .ok_or_else(|| anyhow!("region must look like chr:start-end"))?;

    let (start_s, end_s) = rest
        .split_once('-')
        .ok_or_else(|| anyhow!("region must look like chr:start-end"))?;

    let start1: u32 = start_s.replace(',', "").parse()
        .with_context(|| format!("invalid region start: {start_s}"))?;

    let end1: u32 = end_s.replace(',', "").parse()
        .with_context(|| format!("invalid region end: {end_s}"))?;

    if start1 == 0 || end1 == 0 {
        return Err(anyhow!("region coordinates must be > 0"));
    }

    if end1 < start1 {
        return Err(anyhow!("region end must be >= start"));
    }

    Ok((
        chr.to_string(),
        start1 - 1,
        end1,
    ))
}
