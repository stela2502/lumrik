// crates/bam_tide/src/quantification/bam_collector/config.rs

use std::path::PathBuf;

use clap::Args;

use crate::quantification::cli::QuantMode;

#[derive(Debug, Clone, Args)]
pub struct BamCollectorConfig {
    #[arg(long, help = "Splice index used for gene/transcript assignment.")]
    pub index: PathBuf,

    #[arg(
        long,
        help = "Optional reference genome FASTA used for genome refinement."
    )]
    pub genome: Option<PathBuf>,

    #[arg(long, help = "Optional VCF used for SNP quantification.")]
    pub vcf: Option<PathBuf>,

    #[arg(
        long,
        value_enum,
        default_value_t = QuantMode::Gene,
        help = "Quantification mode."
    )]
    pub quant_mode: QuantMode,

    #[arg(
        long,
        default_value_t = 0,
        help = "Minimum mapping quality accepted for quantification."
    )]
    pub min_mapq: u8,

    #[arg(long, help = "Optional maximum number of accepted reads to process.")]
    pub max_reads: Option<usize>,

    #[arg(
        long,
        default_value_t = false,
        help = "Restrict quantification to read 1 of paired alignments."
    )]
    pub read1_only: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Disable genome refinement even when --genome is supplied."
    )]
    pub no_genome_refine: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Require strand agreement for feature assignment."
    )]
    pub require_strand: bool,

    #[arg(
        long,
        default_value_t = false,
        help = "Require the exact annotated splice-junction chain."
    )]
    pub require_exact_junction_chain: bool,

    #[arg(
        long,
        default_value_t = 0,
        help = "Maximum allowed 5-prime overhang in bp."
    )]
    pub max_5p_overhang_bp: u32,

    #[arg(
        long,
        default_value_t = 0,
        help = "Maximum allowed 3-prime overhang in bp."
    )]
    pub max_3p_overhang_bp: u32,

    #[arg(long, default_value_t = 0, help = "Maximum allowed intronic gap size.")]
    pub allowed_intronic_gap_size: u32,

    #[arg(
        long,
        default_value_t = 5,
        help = "Minimum aligned anchor length required for SNP support."
    )]
    pub snp_min_anchor: u8,

    #[arg(long, help = "Optional BAM output for a streamed mapper input.")]
    pub bam_out: Option<PathBuf>,
}
