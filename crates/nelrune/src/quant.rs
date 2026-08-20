use anyhow::{Context, Result};

use rust_htslib::bam::HeaderView;

use bam_tide::quantification::chunk_processor::ChunkProcessor;
use bam_tide::quantification::cli::QuantMode;
use bam_tide::quantification::job::Job;
use bam_tide::quantification::snp::SnpSideChannel;
use bam_tide::results::QuantData;

use snp_index::Genome;

use crate::cli::Cli;

pub const CHUNK: usize = 2_000_000;

pub fn configure_rayon(
    threads: usize,
) {
    if threads == 0 {
        return;
    }

    let _ =
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
}

pub fn load_genome(
    args: &Cli,
) -> Result<Option<Genome>> {
    match &args.genome {
        Some(path) => {
            Genome::from_fasta(path)
                .with_context(|| {
                    format!(
                        "reading genome FASTA {}",
                        path.display()
                    )
                })
                .map(Some)
        }

        None => Ok(None),
    }
}

pub fn load_snp_side_channel(
    args: &Cli,
    header: &HeaderView,
) -> Result<Option<SnpSideChannel>> {
    let Some(vcf) =
        &args.vcf
    else {
        return Ok(None);
    };

    println!("Loading SNP index");

    let chr_names:
        Vec<String> =
        (0..header.target_count())
            .map(|tid| {
                String::from_utf8_lossy(
                    header.tid2name(tid),
                )
                .to_string()
            })
            .collect();

    let chr_lengths:
        Vec<u32> =
        (0..header.target_count())
            .map(|tid| {
                header
                    .target_len(tid)
                    .unwrap_or(0)
                    as u32
            })
            .collect();

    SnpSideChannel::from_vcf_path(
        vcf,
        chr_names,
        chr_lengths,
        args.snp_min_anchor,
    )
    .map(Some)
}
