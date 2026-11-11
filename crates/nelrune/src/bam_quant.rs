//! Re-quantify a Nelrune mapper BAM, preserving encoded QNAME cell/UMI identity.
//! No normalization, mapper, completion-marker or FASTQ inputs are involved.
use std::fs;
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use anyhow::{Context, Result, bail};
use clap::Parser;
use bam_tide::FeatureTagCounts;
use bam_tide::quantification::bam_collector::BamCollectorConfig;
use lumrik_status::{public_hostname, spawn_status_server};
use nelrune::progress::RunProgress;

#[derive(Debug, Parser)]
#[command(version = env!("CARGO_PKG_VERSION"), about = "Re-quantify an existing Nelrune mapper BAM without remapping")]
struct BamQuantCli {
    /// Existing Nelrune mapper BAM (encoded cell/UMI QNAMEs or CB/UB tags).
    #[arg(long)]
    mapped_bam: lumrik_clap::ReadableFile,
    #[command(flatten)]
    collector: BamCollectorConfig,
    #[arg(long)]
    min_umi_count: Option<usize>,
    /// Cell barcode length for matrix export (10x 3' v3/v4: 16).
    #[arg(long, default_value_t = 16)]
    cell_barcode_len: usize,
    #[arg(long, default_value_t = 0)]
    threads: usize,
    #[arg(long, short)]
    outpath: lumrik_clap::OutputDir,
    #[arg(long, default_value_t = 8787)]
    health_port: u16,
    #[arg(long)]
    health_hostname: Option<String>,
    #[arg(long, default_value_t = false)]
    no_health_server: bool,
}

pub fn run() -> Result<()> {
    let args = BamQuantCli::parse_from(std::env::args().skip(1));
    if args.cell_barcode_len == 0 { bail!("--cell-barcode-len must be positive"); }
    if args.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global();
    }
    fs::create_dir_all(&args.outpath)?;
    super::check_bam(&args.mapped_bam)?;
    let mut progress = RunProgress::new();
    progress.start_timer("nelrune/startup");
    progress.open_log(args.outpath.join("nelrune.log"))?;
    let _server = if args.no_health_server { None } else {
        let addr = SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), args.health_port);
        let server = spawn_status_server(progress.state_handle(), addr)
            .context("starting health server")?;
        eprintln!("[nelrune bam-quant] health server: http://{}:{}",
            public_hostname(args.health_hostname.as_deref()), server.addr().port());
        Some(server)
    };
    progress.stage("reusing existing mapped BAM; no normalization or mapping");
    super::quantify_mapped_bam(
        &args.mapped_bam, &args.collector, &args.outpath,
        args.min_umi_count, None, args.cell_barcode_len,
        FeatureTagCounts::default(), &mut progress,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn bam_quant_has_no_mapper_or_normalizer_requirements() {
        let command = BamQuantCli::command();
        command.clone().debug_assert();
        let names: Vec<_> = command.get_arguments().map(|a| a.get_id().as_str().to_owned()).collect();
        assert!(names.iter().any(|n| n == "mapped_bam"));
        assert!(names.iter().any(|n| n == "index"));
        assert!(!names.iter().any(|n| n == "mapper" || n == "mapper_index" || n == "bam"));
    }
}
