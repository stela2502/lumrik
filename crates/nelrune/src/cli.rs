use std::path::PathBuf;

use anyhow::{bail, Result};
use clap::Parser;

use bam_tide::quantification::bam_collector::BamCollectorConfig;
use bam_tide::AdditionalFeatureSource;
use sc_mapper::StreamingMapperCli;
use sc_primer::PrimerCli;


#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Stream ONT or Illumina single-cell reads through normalization, mapping and bam-tide quantification"
)]
pub struct Cli {
    /// Illumina barcode / UMI FASTQ(s).
    #[arg(long, num_args = 1..)]
    pub r1: Vec<PathBuf>,

    /// Illumina biological insert FASTQ(s).
    #[arg(long, num_args = 1..)]
    pub r2: Vec<PathBuf>,

    /// ONT/Dorado BAM input. Mutually exclusive with --r1/--r2.
    #[arg(long)]
    pub bam: Option<PathBuf>,

    #[command(flatten)]
    pub primer: PrimerCli,

    #[command(flatten)]
    pub mapper: StreamingMapperCli,

    /// All bam-tide quantification options live here; Nelrune does not duplicate them.
    #[command(flatten)]
    pub bam_collector: BamCollectorConfig,

    /// Additional short-feature references searched before genomic mapping.
    ///
    /// Built-ins:
    ///   bd_sample_human
    ///   bd_sample_mouse
    ///
    /// Any other value is interpreted as a FASTA path. For FASTA sources,
    /// the filename stem defines the output/feature type and FASTA headers
    /// define the individual feature names (for example hto.fa -> type hto).
    #[arg(long, num_args = 1..)]
    pub additional_features: Vec<AdditionalFeatureSource>,

    /// Minimum supporting 8-mer hits for additional-feature mapping.
    #[arg(long, default_value_t = 4)]
    pub additional_feature_min_hits: u32,

    /// Minimum Illumina insert length emitted by the normalizer.
    #[arg(long, default_value_t = 20)]
    pub min_insert_len: usize,

    /// Minimum ONT transcript length emitted by the normalizer.
    #[arg(long, default_value_t = 20)]
    pub min_transcript_len: usize,

    /// Minimum cell UMI count retained during sparse-matrix export.
    #[arg(long, default_value_t = 400)]
    pub min_cell_counts: usize,

    /// Worker threads used by the normalizers / export helpers.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Output directory.
    #[arg(long, short)]
    pub outpath: PathBuf,

    /// Port for the live health/progress server.
    #[arg(long, default_value_t = 8787)]
    pub health_port: u16,

    /// Hostname shown in the live-server URL. Useful on clusters when HOSTNAME is not externally resolvable.
    #[arg(long)]
    pub health_hostname: Option<String>,

    /// Disable the live health/progress server.
    #[arg(long, default_value_t = false)]
    pub no_health_server: bool,
}

impl Cli {
    pub fn validate(&self) -> Result<()> {
        let has_ont = self.bam.is_some();
        let has_illumina = !self.r1.is_empty() || !self.r2.is_empty();

        match (has_ont, has_illumina) {
            (true, true) => {
                bail!("choose exactly one input mode: --bam OR --r1/--r2")
            }
            (false, false) => {
                bail!("no input supplied: use --bam for ONT or --r1/--r2 for Illumina")
            }
            _ => {}
        }

        if has_illumina {
            if self.r1.is_empty() || self.r2.is_empty() {
                bail!("Illumina input requires both --r1 and --r2")
            }

            if self.r1.len() != self.r2.len() {
                bail!(
                    "--r1 and --r2 must be supplied the same number of times (got {} R1 files and {} R2 files)",
                    self.r1.len(),
                    self.r2.len(),
                );
            }
        }

        if self.mapper.mapper_paired {
            bail!(
                "Nelrune currently maps the normalized biological insert as single-end; --mapper-paired is not supported"
            );
        }

        Ok(())
    }
}
