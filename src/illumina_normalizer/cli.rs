use clap::{Parser, ValueEnum};
use std::path::PathBuf;
use fast_tag_mapper::FastMapperCli;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum PrimerRead {
    #[value(name = "r1")]
    R1,

    #[value(name = "r2")]
    R2,
}

impl PrimerRead {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::R1 => "r1",
            Self::R2 => "r2",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum InsertRead {
    /// Use the insert returned by the primer detector.
    ///
    /// This is useful when the biological sequence is on the same read that
    /// contains the primer/cell/UMI grammar.
    #[value(name = "detected")]
    Detected,

    /// Use the synchronized R1 record as the emitted biological read.
    #[value(name = "r1")]
    R1,

    /// Use the synchronized R2 record as the emitted biological read.
    ///
    /// This is the common Illumina single-cell mode:
    /// R1 is barcode/UMI, R2 is biological sequence.
    #[value(name = "r2")]
    R2,
}

impl InsertRead {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Detected => "detected",
            Self::R1 => "r1",
            Self::R2 => "r2",
        }
    }
}

#[derive(Debug, Clone, Parser)]
#[command(
    author,
    version,
    about = "Normalize Illumina FASTQ pairs into mapper FASTQ plus read-tag metadata",
    long_about = "\
Normalize Illumina paired FASTQ reads into one mapper-facing FASTQ record per
accepted molecule plus a read-tag TSV.

The primer grammar is supplied by sc_primer and can describe 10x, BD Rhapsody,
sample tags, feature tags, or custom structures.  Unlike the ONT normalizer,
Illumina normally has one molecule per FASTQ pair, so only the first valid
primer hit is used.

Typical 10x/BD-style use:

    bam-illumina-normalizer \\
      --r1 sample_R1.fastq.gz \\
      --r2 sample_R2.fastq.gz \\
      --out normalized_R2.fastq.gz \\
      --read-tags molecule_tags.tsv \\
      --primer-read r1 \\
      --insert-read r2 \\
      --chemistry bd-v2-384 \\
      --threads 8 \\
      --gzip-level 1

For chemistries where the biological insert is part of the same read as the
primer grammar, use:

    --insert-read detected
"
)]
pub struct Cli {
    #[arg(
        long,
        value_name = "FASTQ[.GZ]",
        help = "Input R1 FASTQ. Usually barcode/UMI for Illumina single-cell libraries."
    )]
    pub r1: PathBuf,

    #[arg(
        long,
        value_name = "FASTQ[.GZ]",
        help = "Input R2 FASTQ. Usually biological insert/read-to-map."
    )]
    pub r2: PathBuf,

    #[arg(
        long,
        short,
        value_name = "FASTQ[.GZ]",
        help = "Output normalized FASTQ for mapping. Usually synchronized R2."
    )]
    pub out: PathBuf,

    #[arg(
        long,
        short = 't',
        value_name = "TSV",
        help = "Output molecule metadata TSV. Contains read_id, original_read_id, orientation, raw_cb, quality_cb, raw_umi, quality_umi, and status."
    )]
    pub read_tags: PathBuf,

    #[arg(
        long,
        value_enum,
        default_value_t = PrimerRead::R1,
        help = "Which input read should be scanned with the sc_primer detector."
    )]
    pub primer_read: PrimerRead,

    #[arg(
        long,
        value_enum,
        default_value_t = InsertRead::R2,
        help = "Which read should become the emitted mapper FASTQ record. Use r2 for ordinary Illumina, detected for same-read primer+insert layouts."
    )]
    pub insert_read: InsertRead,

    #[command(flatten)]
    pub primer: sc_primer::PrimerCli,

    #[command(flatten)]
    pub feature_tags: FastMapperCli,

    #[arg(
        long,
        default_value_t = 20,
        value_name = "BP",
        help = "Minimum emitted insert/read length."
    )]
    pub min_insert_len: usize,

    #[arg(
        long,
        default_value_t = 4,
        value_name = "N",
        help = "Worker threads for paired FASTQ processing."
    )]
    pub threads: usize,

    #[arg(
        long,
        default_value_t = 1,
        value_name = "0-9",
        help = "Gzip compression level for FASTQ output."
    )]
    pub gzip_level: u32,

    #[arg(
        long,
        default_value_t = false,
        help = "Write plain FASTQ instead of gzip-compressed FASTQ."
    )]
    pub no_gzip: bool,
}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
