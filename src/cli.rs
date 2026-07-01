//cli.rs
// src/cli.rs
use crate::bed_data::Normalize;
use clap::{Parser, ValueEnum};

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum AnalysisType {
    Bulk,
    #[value(name = "single-cell")]
    SingleCell,
}

/// Shared CLI options for coverage exporters (bedGraph / bigWig)
#[derive(Parser, Debug, Clone)]
pub struct CoverageCli {
    /// Input BAM file (no need to sort or index)
    #[arg(short = 'b', long)]
    pub bam: String,

    /// Output file (bedGraph or BigWig depending on binary)
    #[arg(short = 'o', long)]
    pub outfile: String,

    /*
    /// tag name for the CELL information (default CB for velocity default - change to CR for CellRanger)
    #[arg(short = 'c', long, default_value = "CB")]
    pub cell_tag: String,

    /// tag name for the UMI information (default UB for velocity default - change to UR for CellRanger)
    #[arg(short = 'u', long, default_value = "UB")]
    pub umi_tag: String,

    /// Collect single cell info or bulk
    #[arg(short = 'a', long, value_enum, default_value_t = AnalysisType::Bulk)]
    pub analysis_type: AnalysisType,
    */
    /// Normalize the data somehow
    #[arg(
        short = 'n',
        long = "normalize",
        default_value = "not",
        long_help = r#"Coverage normalization method.

Available modes:
  - not   : no normalization (raw read counts per bin)
  - rpkm  : Reads Per Kilobase per Million mapped reads
  - cpm   : Counts Per Million mapped reads
  - bpm   : Bins Per Million mapped reads (bin-level CPM)
  - rpgc  : Reads Per Genomic Content (1× genome coverage)

Default: not"#
    )]
    pub normalize: Normalize,

    /// Bin width for coverage calculation
    #[arg(short = 'w', long, default_value_t = 50)]
    pub width: u32,

    /// Collect only R1 areas
    #[arg(long)]
    pub only_r1: bool,

    /// Minimum mapping quality to include a read
    #[arg(long, default_value_t = 0)]
    pub min_mapping_quality: u8,

    /// Include secondary alignments
    #[arg(long, default_value_t = true)]
    pub include_secondary: bool,

    /// Include supplementary alignments
    #[arg(long, default_value_t = true)]
    pub include_supplementary: bool,

    /// Include duplicate-marked reads
    #[arg(long, default_value_t = true)]
    pub include_duplicates: bool,

    #[arg(
        long,
        long_help = r#"Exclude reads with ANY of these SAM flag bits set
(equivalent to deeptools --samFlagExclude).

The value is a bitmask of SAM flags. Any read with
(read_flag & mask) != 0 will be discarded.

Examples:
  - 256  -> exclude secondary alignments
  - 512  -> exclude QC-failed reads
  - 1024 -> exclude PCR/optical duplicates
  - 2048 -> exclude supplementary alignments
  - 2816 -> exclude secondary + QC-fail + supplementary

Default: None (no flag-based exclusion, matches bamCoverage defaults)"#
    )]
    pub sam_flag_exclude: Option<u16>,

    #[arg(
        long,
        long_help = r#"Include only reads that have ALL of these SAM flag bits set.
Applied after the exclusion test
(equivalent to deeptools --samFlagInclude).

The value is a bitmask of SAM flags. A read is kept only if:
  (read_flag & mask) == mask

Examples:
  - 64   -> include only read1
  - 128  -> include only read2
  - 2    -> include only properly paired reads

Default: None (no include constraint, matches bamCoverage defaults)"#
    )]
    pub sam_flag_include: Option<u16>,
}
