//! Convert a transcriptome-coordinate BAM/SAM/CRAM alignment back to genome coordinates.
//!
//! This tool is intended for alignments produced by mapping reads against a
//! transcriptome FASTA built from the same GTF/GFF + genome FASTA supplied here.
//!
//! The conversion is deliberately simple:
//!
//! 1. Build a `gtf_splice_index::SpliceIndex` from the annotation and genome FASTA.
//! 2. Replace transcriptome `@SQ` records in the BAM header with genome `@SQ` records.
//! 3. For each mapped read:
//!    - resolve the transcript name from the input BAM reference name,
//!    - find the matching transcript in the `SpliceIndex`,
//!    - translate transcriptome POS+CIGAR to genome POS+CIGAR using transcript exons,
//!    - reverse-complement sequence/quality for minus-strand transcripts,
//!    - update `tid`, `pos`, `cigar`, and write the record.
//! 4. Unmapped reads are passed through unchanged.
//!
//! Reads whose transcript name is not found in the index are skipped. This is safer
//! than writing incorrectly placed reads.
//!
//! Example:
//!
//! ```bash
//! bam-transcriptome-to-genome \
//!   --bam transcriptome.bam \
//!   --gtf genes.gtf \
//!   --genome genome.fa \
//!   --out genome.bam \
//!   --threads 8
//! ```
//!
//! The output BAM is not automatically indexed. Run:
//!
//! ```bash
//! samtools index genome.bam
//! ```

use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::bam;
use rust_htslib::bam::{Read, Reader, Writer};
use std::path::PathBuf;

use bam_tide::transcriptome_to_genome::BamTranscriptomeMapper;

/// Convert transcriptome-coordinate BAM records back to genome coordinates.
#[derive(Debug, Parser)]
#[command(author, version, about)]
struct Cli {
    /// Input BAM/SAM/CRAM aligned to transcriptome reference names.
    #[arg(short = 'b', long = "bam")]
    bam: PathBuf,

    /// GTF/GFF annotation used to build the transcriptome FASTA.
    ///
    /// Transcript IDs in this annotation must match the input BAM reference names.
    #[arg(short = 'g', long = "gtf")]
    gtf: PathBuf,

    /// Genome FASTA used to build the transcriptome FASTA.
    ///
    /// The chromosome names must match the annotation seqnames.
    #[arg(short = 'f', long = "genome")]
    genome: PathBuf,

    /// Output genome-coordinate BAM.
    #[arg(short = 'o', long = "out")]
    out: PathBuf,

    /// Number of BAM reader/writer threads.
    #[arg(short = 't', long = "threads", default_value_t = 4)]
    threads: usize,

    /// Write skipped/failed-read statistics to stderr.
    #[arg(long = "quiet", default_value_t = false)]
    quiet: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    run(cli)
}

fn run(cli: Cli) -> Result<()> {
    if cli.threads == 0 {
        anyhow::bail!("--threads must be at least 1");
    }

    let mut reader = Reader::from_path(&cli.bam)
        .with_context(|| format!("opening input BAM/SAM/CRAM {}", cli.bam.display()))?;

    reader
        .set_threads(cli.threads)
        .with_context(|| format!("setting {} reader threads", cli.threads))?;

    let mut mapper = BamTranscriptomeMapper::new(&cli.gtf, &cli.genome)
        .with_context(|| {
            format!(
                "building transcriptome-to-genome mapper from annotation {} and genome {}",
                cli.gtf.display(),
                cli.genome.display()
            )
        })?;

    let out_header = mapper
        .make_header(reader.header())
        .context("building genome-coordinate BAM header")?;

    let mut writer = Writer::from_path(&cli.out, &out_header, bam::Format::Bam)
        .with_context(|| format!("creating output BAM {}", cli.out.display()))?;

    writer
        .set_threads(cli.threads)
        .with_context(|| format!("setting {} writer threads", cli.threads))?;

    // Keep the old header alive and immutable while streaming records.
    // It is needed because input TIDs refer to transcriptome reference names,
    // while output TIDs refer to genome chromosome names.
    let old_header = reader.header().clone();

    for record in reader.records() {
        let mut record = record.with_context(|| {
            format!("reading record from input {}", cli.bam.display())
        })?;

        let write_record = mapper
            .map_record_in_place(&mut record, &old_header)
            .with_context(|| {
                format!(
                    "mapping read {:?} from transcriptome to genome coordinates",
                    String::from_utf8_lossy(record.qname())
                )
            })?;

        if write_record {
            writer
                .write(&record)
                .with_context(|| format!("writing record to {}", cli.out.display()))?;
        }
    }

    if !cli.quiet {
        eprintln!("{}", mapper.stats());
    }

    Ok(())
}
