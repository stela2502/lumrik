//cli.rs
// src/quantification/cli.rs

use std::str::FromStr;
use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use read_tag_table::ReadTagTableCli;



#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum QuantMode {
    Gene,
    Transcript,
}


#[derive(Debug, Clone, Copy)]
pub struct BamAuxTag(pub [u8; 2]);

impl FromStr for BamAuxTag {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let bytes = value.as_bytes();

        if bytes.len() != 2 {
            return Err(format!(
                "BAM aux tag '{value}' must be exactly two ASCII characters"
            ));
        }

        if !bytes.iter().all(|b| b.is_ascii_alphanumeric()) {
            return Err(format!(
                "BAM aux tag '{value}' must contain only ASCII letters/numbers"
            ));
        }

        Ok(Self([bytes[0], bytes[1]]))
    }
}


#[derive(Parser, Debug, Clone)]
#[command(
    name = "bam-quant",
    about = "Quantify 10x BAM against splice index into scdata, optionally collecting SNP ref/alt matrices"
)]
pub struct QuantCli {
    /// Input BAM file(s) from a single cell mapping.
    #[arg(short = 'b', long = "bam", required = true, num_args = 1..)]
    pub bam: Vec<PathBuf>,

    /// Splice index path (built from GTF beforehand)
    #[arg(long, short)]
    pub index: PathBuf,

    /// Outpath for the 10x mtx-formatted outfiles
    #[arg(long, short)]
    pub outpath: PathBuf,

    /// Split Intronic from rest.
    ///
    /// This is currently not recommended as exon/intron detection seems to be too strict
    /// for normal sequencing data.
    #[arg(long, short, default_value_t = false)]
    pub split_intronic: bool,

    /// Minimum MAPQ
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Use only read1 (recommended for 10x; reduces duplicate mate-counting noise)
    #[arg(long, default_value_t = false)]
    pub read1_only: bool,

    /// Rayon thread count (0 = default)
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Collect Gene or Transcript names
    #[arg(long, value_enum, default_value_t = QuantMode::Gene)]
    pub quant_mode: QuantMode,

    /// Max reads to process (debug/dev)
    #[arg(long)]
    pub max_reads: Option<usize>,

    /// Min read counts per reported cell (debug/dev)
    #[arg(long, default_value_t = 400)]
    pub min_cell_counts: usize,

    /// Optional reference genome FASTA.
    ///
    /// If supplied, BAM-derived AlignedRead objects are refined against the genome.
    #[arg(long)]
    pub genome: Option<PathBuf>,

    #[command(flatten)]
    pub read_tags: ReadTagTableCli,

    /// Optional SNP VCF.
    ///
    /// If supplied, SNP ref/alt matrices are written in addition to the normal
    /// gene/transcript matrix. Requires --genome.
    #[arg(long)]
    pub vcf: Option<PathBuf>,

    /// Minimum SNP anchor/support passed to snp_index.match_read().
    #[arg(long, default_value_t = 20)]
    pub snp_min_anchor: u8,

    /// Disable genome-based AlignedRead refinement even if --genome is supplied.
    #[arg(long, default_value_t = false)]
    pub no_genome_refine: bool,

    // ------------------------------
    // MatchOptions exposed to user
    // ------------------------------
    /// If true, require read blocks to be on a compatible strand.
    #[arg(long, default_value_t = false)]
    pub require_strand: bool,

    /// If true, require the read to have the exact same splice junction chain as the transcript.
    #[arg(long, default_value_t = false)]
    pub require_exact_junction_chain: bool,

    /// Maximum allowed 5′ overhang (bp). If exceeded -> OverhangTooLarge.
    #[arg(long, default_value_t = 100)]
    pub max_5p_overhang_bp: u32,

    /// Maximum allowed 3′ overhang (bp). If exceeded -> OverhangTooLarge.
    #[arg(long, default_value_t = 100)]
    pub max_3p_overhang_bp: u32,

    /// Allowed sequencing error gap. If exceeded -> JunctionMismatch.
    #[arg(long, default_value_t = 5)]
    pub allowed_intronic_gap_size: u32,

    /// BAM aux tag containing the cell barcode.
    ///
    /// Examples:
    ///   CB (10x corrected)
    ///   CR (10x raw)
    ///   XC (custom)
    #[arg(long, default_value = "CB")]
    pub cell_tag: BamAuxTag,

    /// BAM aux tag containing the UMI.
    ///
    /// Examples:
    ///   UB (10x corrected)
    ///   UR (10x raw)
    ///   XM (custom)
    #[arg(long, default_value = "UB")]
    pub umi_tag: BamAuxTag,

}
