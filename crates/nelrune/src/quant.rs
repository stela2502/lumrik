use anyhow::{Context, Result};

use rust_htslib::bam::HeaderView;

use bam_tide::quantification::cli::QuantMode;
use bam_tide::quantification::snp::SnpSideChannel;
use bam_tide::results::QuantData;

use snp_index::Genome;

use crate::cli::Cli;

use bam_tide::index::{
    GeneFeatureIndex,
    TranscriptFeatureIndex,
};

use gtf_splice_index::SpliceIndex;



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


pub fn write_quantification(
    args: &Cli,
    idx: &SpliceIndex,
    snp: Option<&SnpSideChannel>,
    data: &mut QuantData,
) -> Result<()> {
    match args.quant_mode {
        QuantMode::Gene => {
            let feature_index =
                GeneFeatureIndex::new(idx);

            data.write_finalized(
                &args.outpath,
                &feature_index,
                snp.map(|s| &s.index),
            )
            .map_err(anyhow::Error::msg)
            .context(
                "writing gene quantification"
            )?;
        }

        QuantMode::Transcript => {
            let feature_index =
                TranscriptFeatureIndex::new(idx);

            data.write_finalized(
                &args.outpath,
                &feature_index,
                snp.map(|s| &s.index),
            )
            .map_err(anyhow::Error::msg)
            .context(
                "writing transcript quantification"
            )?;
        }
    }

    Ok(())
}