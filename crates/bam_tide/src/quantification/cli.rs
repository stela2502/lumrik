//cli.rs
// src/quantification/cli.rs

use std::path::PathBuf;
use std::str::FromStr;

use clap::{Args, ValueEnum};

use read_tag_table::ReadTagTableCli;

use crate::{
    cli::AnalysisType,
    quantification::bam_collector::config::{
        DEFAULT_ALLOWED_INTRONIC_GAP_SIZE, DEFAULT_MAX_3P_OVERHANG_BP, DEFAULT_MAX_5P_OVERHANG_BP,
        DEFAULT_SNP_MIN_ANCHOR,
    },
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum QuantMode {
    Gene,
    Transcript,
}

impl QuantMode {
    pub fn splice_match_mode(self) -> gtf_splice_index::SpliceMatchMode {
        match self {
            Self::Gene => gtf_splice_index::SpliceMatchMode::Gene,
            Self::Transcript => gtf_splice_index::SpliceMatchMode::Transcript,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum CellCallingMode {
    /// Historical fixed UMI threshold (`--min-cell-counts`).
    Fixed,
    /// sc-beacon two-Gaussian mixture over log10 per-cell exonic UMI totals.
    Beacon,
}

#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum GrammarSelection {
    Gex,
    Vdj,
    Other,
    All,
}

impl GrammarSelection {
    pub fn accepts(self, grammar_type: sc_primer::GrammarType) -> bool {
        match self {
            Self::All => true,
            Self::Gex => grammar_type == sc_primer::GrammarType::Gex,
            Self::Vdj => grammar_type == sc_primer::GrammarType::Vdj,
            Self::Other => grammar_type == sc_primer::GrammarType::Other,
        }
    }
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

#[derive(Args, Debug, Clone)]
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

    /// Cell-calling strategy. `fixed` preserves the historical hard UMI cutoff;
    /// `beacon` uses sc-beacon barcode-rank knee detection.
    #[arg(long, value_enum, default_value_t = CellCallingMode::Fixed)]
    pub cell_calling: CellCallingMode,

    /// Minimum UMI count for `--cell-calling fixed`.
    #[arg(long, default_value_t = 400)]
    pub min_umi_count: usize,

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
    #[arg(long, default_value_t = DEFAULT_SNP_MIN_ANCHOR)]
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
    #[arg(long, default_value_t = DEFAULT_MAX_5P_OVERHANG_BP)]
    pub max_5p_overhang_bp: u32,

    /// Maximum allowed 3′ overhang (bp). If exceeded -> OverhangTooLarge.
    #[arg(long, default_value_t = DEFAULT_MAX_3P_OVERHANG_BP)]
    pub max_3p_overhang_bp: u32,

    /// Allowed sequencing error gap. If exceeded -> JunctionMismatch.
    #[arg(long, default_value_t = DEFAULT_ALLOWED_INTRONIC_GAP_SIZE)]
    pub allowed_intronic_gap_size: u32,

    /// Primer-grammar provenance to quantify from Lumrik mapper QNAMEs.
    /// Legacy QNAME/CB+UB input without provenance is treated as GEX.
    #[arg(long, value_enum, default_value_t = GrammarSelection::Gex)]
    pub grammar_type: GrammarSelection,

    /// Optional read grammar for BAMs without cell/UMI metadata.
    ///
    /// `NONE` or another grammar without CELL/UMI enables sequence-based
    /// paired-fragment deduplication and therefore requires query-name sorted
    /// BAM input. If omitted, Lumrik-encoded QNAMEs, read-tag tables, or BAM
    /// CB/UB tags are used as before.
    #[arg(long)]
    pub primer_structure: Option<String>,

    /// Quantify a single-cell BAM using cell/UMI metadata, or a bulk BAM
    /// without cell/UMI tags.
    ///
    /// In bulk mode all reads from one BAM are assigned to one synthetic
    /// sample/cell and the read QNAME supplies a stable pseudo-UMI so paired
    /// mates count as one fragment. Multiple BAM inputs become separate
    /// synthetic samples.
    #[arg(long, value_enum, default_value_t = AnalysisType::SingleCell)]
    pub analysis_type: AnalysisType,

    /// BAM aux tag containing the cell barcode (single-cell mode only).
    ///
    /// Examples:
    ///   CB (10x corrected)
    ///   CR (10x raw)
    ///   XC (custom)
    #[arg(long, default_value = "CB")]
    pub cell_tag: BamAuxTag,

    /// BAM aux tag containing the UMI (single-cell mode only).
    ///
    /// Examples:
    ///   UB (10x corrected)
    ///   UR (10x raw)
    ///   XM (custom)
    #[arg(long, default_value = "UB")]
    pub umi_tag: BamAuxTag,
}
