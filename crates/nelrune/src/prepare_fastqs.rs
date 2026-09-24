use std::fs;
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use bam_tide::AdditionalFeatureSource;
use bam_tide::illumina_normalizer::cli::{InsertRead, PrimerRead};
use bam_tide::illumina_normalizer::{IlluminaNormalizer, IlluminaNormalizerConfig};
use clap::Parser;
use lumrik_status::{public_hostname, spawn_status_server};
use nelrune::progress::RunProgress;
use sc_primer::PrimerCli;

#[derive(Debug, Parser)]
#[command(about = "Prepare mapper-ready single-cell FASTQ shards without running STAR")]
pub struct PrepareFastqsCli {
    #[arg(long, num_args = 1..)]
    r1: Vec<PathBuf>,
    #[arg(long, num_args = 1..)]
    r2: Vec<PathBuf>,
    #[command(flatten)]
    primer: PrimerCli,
    #[arg(long, num_args = 1..)]
    additional_features: Vec<AdditionalFeatureSource>,
    #[arg(long, default_value_t = 4)]
    additional_feature_min_hits: u32,
    #[arg(long, default_value_t = 20)]
    min_insert_len: usize,
    #[arg(long)]
    max_reads: Option<usize>,
    #[arg(long, default_value_t = 0)]
    threads: usize,
    #[arg(long, default_value_t = 8787)]
    health_port: u16,
    #[arg(long)]
    health_hostname: Option<String>,
    #[arg(long, default_value_t = false)]
    no_health_server: bool,
    #[arg(long, short)]
    outpath: PathBuf,
}

pub fn run() -> Result<()> {
    let args = PrepareFastqsCli::parse_from(std::env::args().skip(1));
    if args.r1.is_empty() || args.r2.is_empty() || args.r1.len() != args.r2.len() {
        bail!("prepare-fastqs requires matching --r1 and --r2 lists");
    }

    fs::create_dir_all(&args.outpath)
        .with_context(|| format!("creating {}", args.outpath.display()))?;
    super::configure_rayon(args.threads);

    let mut progress = RunProgress::new();
    progress.open_log(args.outpath.join("nelrune.log"))?;
    let _health_server = if args.no_health_server {
        progress.stage("health server disabled");
        None
    } else {
        let health_addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), args.health_port);
        let server = spawn_status_server(progress.state_handle(), health_addr)?;
        let hostname = public_hostname(args.health_hostname.as_deref());
        let url = format!("http://{}:{}", hostname, server.addr().port());
        progress.set_public_url(url.clone());
        eprintln!("[nelrune prepare-fastqs] health server: {url}");
        Some(server)
    };

    progress.stage("loading chemistry");
    let primer = args.primer.detector().map_err(anyhow::Error::msg)?;
    let cell_barcode_len = primer.cell_len();
    let config = IlluminaNormalizerConfig {
        out: args.outpath.join("unused.fastq"),
        read_tags: args.outpath.join("unused.read_tags.tsv"),
        primer_read: PrimerRead::R1,
        insert_read: InsertRead::R2,
        primer: primer.clone(),
        additional_features: args.additional_features.clone(),
        additional_feature_min_hits: args.additional_feature_min_hits,
        min_insert_len: args.min_insert_len,
        max_reads: args.max_reads,
        threads: args.threads,
        gzip_level: 1,
        gzip: true,
    };

    let mut normalizer = IlluminaNormalizer::new(config)?;
    progress.stage("preparing FASTQ shards");
    let inputs: Vec<_> = args.r1.into_iter().zip(args.r2).collect();
    let shard_dir = args.outpath.join("prepared_fastqs");
    let shards = normalizer.prepare_fastqs_sharded(&inputs, &shard_dir, |stats| {
        progress.update_from_mapping_info(stats);
    })?;
    progress.update_from_mapping_info(normalizer.stats());

    let feature_path = args.outpath.join("feature_observations.bin");
    normalizer
        .take_feature_tag_counts()
        .save_observations(&feature_path)?;

    fs::write(
        args.outpath.join("prepare-report.txt"),
        normalizer.stats().to_string(),
    )?;
    let mut manifest = format!(
        "format\tnelrune-prepare-v1\ncell_barcode_len\t{cell_barcode_len}\nfeature_observations\t{}\n",
        feature_path.display()
    );
    for grammar in primer.grammars() {
        manifest.push_str(&format!(
            "grammar\t{}\t{}\n",
            grammar.name,
            grammar.grammar_type.code()
        ));
    }
    for shard in &shards {
        manifest.push_str(&format!("fastq\t{}\n", shard.display()));
    }
    fs::write(args.outpath.join("prepare-manifest.tsv"), manifest)?;

    eprintln!(
        "[nelrune prepare-fastqs] wrote {} FASTQ shards to {}",
        shards.len(),
        shard_dir.display()
    );
    progress.finish();
    progress.write_final_status(&args.outpath)?;
    Ok(())
}
