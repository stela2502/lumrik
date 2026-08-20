//cli.rs

use std::path::PathBuf;
use std::str::FromStr;

use anyhow::{bail, Result};
use clap::Parser;

use bam_tide::quantification::cli::QuantMode;

use sc_beacon::cli::GuideModelCli;
use sc_mapper::StreamingMapperCli;
use sc_primer::PrimerCli;

#[derive(Debug, Clone)]
pub enum FastFeatureSource {
    BdSampleHuman,
    BdSampleMouse,
    Fasta(PathBuf),
}


impl FromStr for FastFeatureSource {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "bd_sample_human" => Ok(Self::BdSampleHuman),
            "bd_sample_mouse" => Ok(Self::BdSampleMouse),
            "" => bail!("empty fast-feature source"),
            path => Ok(Self::Fasta(PathBuf::from(path))),
        }
    }
}

#[derive(Parser, Debug)]
pub struct Cli {
    /// Barcode / UMI read.
    #[arg(long)]
    pub r1: PathBuf,

    /// Biological insert read.
    #[arg(long)]
    pub r2: PathBuf,

    // --------------------------------------------------------
    // Chemistry
    // --------------------------------------------------------

    #[command(flatten)]
    pub primer: PrimerCli,

    // --------------------------------------------------------
    // Genomic mapping
    // --------------------------------------------------------

    #[command(flatten)]
    pub mapper: StreamingMapperCli,

    /// Gene stats model
    #[command(flatten)]
    pub beacon: GuideModelCli,

    /// bam-tide splice index.
    #[arg(long, short)]
    pub index: PathBuf,


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
    pub fast_features: Vec<FastFeatureSource>,

    /// Minimum supporting 8-mer hits for fast feature mapping.
    #[arg(long, default_value_t = 4)]
    pub fast_feature_min_hits: u32,

    // --------------------------------------------------------
    // Quantification
    // --------------------------------------------------------

    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    #[arg(
        long,
        value_enum,
        default_value_t = QuantMode::Gene
    )]
    pub quant_mode: QuantMode,

    #[arg(long, default_value_t = 400)]
    pub min_cell_counts: usize,

    #[arg(long)]
    pub max_reads: Option<usize>,

    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    // --------------------------------------------------------
    // Optional genome / SNP support
    // --------------------------------------------------------

    #[arg(long)]
    pub genome: Option<PathBuf>,

    #[arg(long)]
    pub vcf: Option<PathBuf>,

    #[arg(long, default_value_t = 20)]
    pub snp_min_anchor: u8,

    #[arg(long, default_value_t = false)]
    pub no_genome_refine: bool,

    // --------------------------------------------------------
    // Splice matching
    // --------------------------------------------------------

    #[arg(long, default_value_t = false)]
    pub require_strand: bool,

    #[arg(long, default_value_t = false)]
    pub require_exact_junction_chain: bool,

    #[arg(long, default_value_t = 100)]
    pub max_5p_overhang_bp: u32,

    #[arg(long, default_value_t = 100)]
    pub max_3p_overhang_bp: u32,

    #[arg(long, default_value_t = 5)]
    pub allowed_intronic_gap_size: u32,

    // --------------------------------------------------------
    // Output
    // --------------------------------------------------------

    #[arg(long, short)]
    pub outpath: PathBuf,
}