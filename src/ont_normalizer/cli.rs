use clap::Parser;
use std::path::PathBuf;
use sc_primer::PrimerCli;
use fast_tag_mapper::FastMapperCli;


#[derive(Debug, Clone, Parser)]
#[command(
    author,
    version,
    about = "Normalize messy ONT/Dorado BAM reads into one 10x-style molecule per FASTQ read",
    long_about = "\
Normalize messy ONT/Dorado BAM reads into clean FASTQ records.

The tool searches each BAM read in both orientations for 10x-style 3' barcode
cassettes:

    adapter + cell barcode + UMI + polyT + transcript

Each detected molecule is written as one FASTQ record:

    original_read_name/mol<N>

The output sequence is normalized to the expected molecule orientation. Reads
without a valid cassette are not emitted, but they are counted in the final
summary. A TSV sidecar file records the extracted CB/UMI, qualities, coordinates,
orientation, and status for every emitted molecule.

The normalizer does NOT perform whitelist correction. CB and UMI are raw slices
from the read. Barcode correction can be done later.

Typical use:

    bam-ont-normalizer \\
      --bam dorado.bam \\
      --out normalized.fastq.gz \\
      --tags molecule_tags.tsv \\
      --threads 8 \\
      --gzip-level 1
"
)]
pub struct Cli {
    #[arg(
        long,
        short,
        value_name = "BAM",
        help = "Input Dorado/ONT BAM file. The BAM may be unmapped; query sequence and qualities are used."
    )]
    pub bam: PathBuf,

    #[arg(
        long,
        short,
        value_name = "FASTQ[.GZ]",
        help = "Output normalized FASTQ file. By default this is gzip-compressed unless --no-gzip is set."
    )]
    pub out: PathBuf,

    #[arg(
        long,
        short,
        value_name = "TSV",
        help = "Output molecule metadata TSV. Contains one row per emitted molecule with CB/UMI, qualities, coordinates, orientation, and status."
    )]
    pub read_tags: PathBuf,

    #[command(flatten)]
    pub primer: PrimerCli,

    #[command(flatten)]
    pub feature_tags: FastMapperCli,

    #[arg(
        long,
        default_value_t = 4,
        value_name = "N",
        help = "Number of HTSlib threads used for BAM reading/decompression."
    )]
    pub threads: usize,

    #[arg(
        long,
        default_value_t = 1,
        value_name = "0-9",
        help = "Gzip compression level for FASTQ output. Level 1 is fast and recommended for large ONT files."
    )]
    pub gzip_level: u32,

    #[arg(
        long,
        default_value_t = false,
        help = "Write plain FASTQ instead of gzip-compressed FASTQ. Useful for speed benchmarking."
    )]
    pub no_gzip: bool,

    #[arg(
        long,
        default_value_t = 20,
        value_name = "BP",
        help = "Minimum transcript/insert length after primer extraction. Shorter molecules are discarded."
    )]
    pub min_transcript_len: usize,
}

impl Cli {
    pub fn parse_args() -> Self {
        Self::parse()
    }
}
