use std::path::PathBuf;

use clap::Parser;

use crate::encoder::model::EncoderConfig;

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
pub enum TruthFeatureMode {
    Gene,
    Transcript,
}

#[derive(Debug, Clone, Parser)]
#[command(
    name = "bam-quant-testdata",
    author,
    version,
    about = "Generate synthetic 10x-style SAM reads for bam-quant tests",
    long_about = "\
Generate artificial SAM reads from a reference FASTA, a serialized SpliceIndex,
and an optional VCF. The output is intended for technical tests of transcript,
strand, unspliced, antisense, sequencing-error, and SNP-handling logic."
)]
pub struct TestDataCli {
    /// Serialized SpliceIndex file built from the test GTF
    #[arg(long, value_name = "FILE")]
    pub gtf: PathBuf,

    /// Reference FASTA file used by the SpliceIndex
    #[arg(long, value_name = "FILE")]
    pub fasta: PathBuf,

    /// Output SAM file
    #[arg(short, long, value_name = "FILE")]
    pub outfile: PathBuf,

    /// Optional VCF with SNPs to hide in generated reads
    #[arg(long, value_name = "FILE")]
    pub vcf: Option<PathBuf>,

    /// Read length
    #[arg(long, default_value_t = 30)]
    pub seq_len: usize,

    /// Reads generated per cell per selected gene
    #[arg(long, default_value_t = 70)]
    pub reads_per_cell: usize,

    /// Fraction of reads emitted as antisense reads
    #[arg(long, default_value_t = 10.0 / 70.0)]
    pub antisense_fraction: f64,

    /// Fraction of reads emitted as unspliced/intronic reads
    #[arg(long, default_value_t = 10.0 / 70.0)]
    pub unspliced_fraction: f64,

    /// Probability of applying a covered VCF SNP to an eligible read
    #[arg(long, default_value_t = 0.20)]
    pub snp_fraction: f64,

    /// Sequencing error probability in the read body
    #[arg(long, default_value_t = 0.005)]
    pub body_error_rate: f64,

    /// Sequencing error probability at the outermost read bases
    #[arg(long, default_value_t = 0.15)]
    pub end_error_rate: f64,

    /// Number of bases at each read end with elevated error probability
    #[arg(long, default_value_t = 4)]
    pub bad_end_bases: usize,

    /// RNG seed for deterministic output
    #[arg(long, default_value_t = 1)]
    pub seed: u64,

    /// 10x cell barcodes to simulate.
    ///
    /// Accepts comma-separated values.
    /// Example:
    /// --cell-barcodes AAAC...-1,TTGC...-1
    #[arg(
        long,
        value_delimiter = ',',
        value_name = "CB",
        default_value = "AAACCGCTCTATCCTA-1,AAACGAATCACCAGAC-1,AAACGAATCCGCGAAT-1"
    )]
    pub cell_barcodes: Vec<String>,

    /// Number of genes randomly sampled per cell.
    #[arg(long, value_name = "N", default_value_t = 100)]
    pub genes_per_cell: usize,

    /// Optional output directory for generator-side truth matrices.
    ///
    /// If set, truth data is written into subfolders:
    /// `exonic/`, `intronic/`, `ref/`, and `alt/`.
    #[arg(long, value_name = "DIR")]
    pub truth_path: Option<std::path::PathBuf>,

    /// Minimum number of exonic molecules required for a cell to be exported
    /// in the truth matrices.
    #[arg(long, value_name = "N", default_value_t = 90)]
    pub min_cell_counts: usize,

    /// Feature namespace used for exonic/intronic truth matrices.
    #[arg(long, value_enum, default_value_t = TruthFeatureMode::Gene)]
    pub truth_feature_mode: TruthFeatureMode,
}

impl TestDataCli {
    pub fn config(&self) -> EncoderConfig {
        EncoderConfig {
            seq_len: self.seq_len,
            reads_per_cell: self.reads_per_cell,
            antisense_fraction: self.antisense_fraction,
            unspliced_fraction: self.unspliced_fraction,
            snp_fraction: self.snp_fraction,
            body_error_rate: self.body_error_rate,
            end_error_rate: self.end_error_rate,
            bad_end_bases: self.bad_end_bases,
            min_cell_counts: self.min_cell_counts,
            truth_path: self.truth_path.clone(),
            seed: self.seed,
            truth_feature_mode: self.truth_feature_mode,
            genes_per_cell: self.genes_per_cell,
            cell_barcodes: if self.cell_barcodes.is_empty() {
                vec!["AAACCGCTCTATCCTA-1".to_string()]
            } else {
                self.cell_barcodes.clone()
            },
        }
    }
}
