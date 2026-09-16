use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use anyhow::{Context, Result};
use clap::Parser;
use fast_tag_mapper::{FastMapperCli, MapStatus};
use flate2::read::MultiGzDecoder;
use mapping_info::MappingInfo;

#[derive(Debug, Parser)]
#[command(author, version, about = "Benchmark FastTagMapper on FASTQ reads")]
struct Cli {
    #[command(flatten)]
    mapper: FastMapperCli,

    /// FASTQ/FASTQ.GZ whose sequence lines are mapped.
    #[arg(long)]
    fastq: PathBuf,

    /// Maximum number of reads to keep in memory and benchmark.
    #[arg(long, default_value_t = 1_000_000)]
    max_reads: usize,

    /// Number of timed passes over the in-memory reads.
    #[arg(long, default_value_t = 5)]
    iterations: usize,

    /// Print the first N mapped reads after the timed benchmark.
    #[arg(long, default_value_t = 0)]
    show_mapped: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    anyhow::ensure!(cli.iterations > 0, "--iterations must be greater than zero");
    anyhow::ensure!(cli.max_reads > 0, "--max-reads must be greater than zero");

    eprintln!("loading {} ...", cli.fastq.display());
    let load_started = Instant::now();
    let reads = load_fastq_sequences(&cli.fastq, cli.max_reads)?;
    let load_elapsed = load_started.elapsed();

    anyhow::ensure!(!reads.is_empty(), "no FASTQ reads loaded");

    eprintln!(
        "loaded {} reads in {:.3}s ({:.0} reads/s; I/O + decompression excluded from benchmark)",
        reads.len(),
        load_elapsed.as_secs_f64(),
        reads.len() as f64 / load_elapsed.as_secs_f64()
    );

    let mapper = cli.mapper.mapper()?;
    eprintln!("features: {}", mapper.feature_count());
    eprintln!("min_hits: {}", mapper.min_hits());
    eprintln!("iterations: {}", cli.iterations);

    let mut timings = Vec::with_capacity(cli.iterations);
    let mut expected_calls = None;

    for iteration in 0..cli.iterations {
        // MappingInfo is intentionally one object per pass, matching normal use.
        // Its ticker/reporting cost is therefore included in map_feature_id().
        let mut info = MappingInfo::new(None, 0.0, 0);
        let started = Instant::now();
        let mut calls = 0usize;
        let mut checksum = 0u64;

        for seq in &reads {
            if let Some(feature_id) = mapper.map_feature_id(seq, &mut info) {
                calls += 1;
                // Keep the returned value observable so the optimizer cannot
                // discard the mapping work.
                checksum = checksum.wrapping_add(feature_id);
            }
        }

        let elapsed = started.elapsed();
        timings.push(elapsed);

        if let Some(expected) = expected_calls {
            anyhow::ensure!(
                calls == expected,
                "mapping call count changed between iterations: {expected} -> {calls}"
            );
        } else {
            expected_calls = Some(calls);
        }

        eprintln!(
            "  pass {}: {:.3}s = {:.0} reads/s = {:.3} us/read (checksum={checksum})",
            iteration + 1,
            elapsed.as_secs_f64(),
            reads.len() as f64 / elapsed.as_secs_f64(),
            elapsed.as_secs_f64() * 1_000_000.0 / reads.len() as f64,
        );
    }

    report_summary(&timings, reads.len(), expected_calls.unwrap_or(0));

    if cli.show_mapped > 0 {
        report_mapped_reads(&mapper, &reads, cli.show_mapped);
    }

    Ok(())
}

fn report_mapped_reads(mapper: &fast_tag_mapper::FastTagMapper, reads: &[Vec<u8>], limit: usize) {
    eprintln!("diagnostic: first {limit} mapped reads (untimed)");

    let mut info = MappingInfo::new(None, 0.0, 0);
    let mut shown = 0usize;

    for (read_index, seq) in reads.iter().enumerate() {
        let MapStatus::Hit {
            feature_id,
            feature_index,
            hits,
        } = mapper.map_status(seq, &mut info)
        else {
            continue;
        };

        let feature = mapper.feature(feature_index);
        let feature_name = feature
            .map(|feature| feature.name.as_str())
            .unwrap_or("<unknown>");
        let feature_type = feature
            .map(|feature| feature.feature_type.as_str())
            .unwrap_or("<unknown>");

        eprintln!(
            "mapped[{}] read_index={} feature_id={} feature_index={} feature_name={} feature_type={} hits={}\n  read={}",
            shown + 1,
            read_index,
            feature_id,
            feature_index,
            feature_name,
            feature_type,
            hits,
            String::from_utf8_lossy(seq),
        );

        shown += 1;
        if shown >= limit {
            break;
        }
    }

    if shown == 0 {
        eprintln!("  no mapped reads found");
    } else if shown < limit {
        eprintln!("  only {shown} mapped reads found in the loaded fixture");
    }
}

fn report_summary(timings: &[Duration], reads: usize, calls: usize) {
    let best = timings.iter().copied().min().unwrap();
    let mean_secs = timings.iter().map(Duration::as_secs_f64).sum::<f64>() / timings.len() as f64;

    eprintln!("FastTagMapper::map_feature_id");
    eprintln!(
        "  calls: {calls}/{reads} ({:.2}%)",
        calls as f64 * 100.0 / reads as f64
    );
    eprintln!(
        "  best: {:.3}s = {:.0} reads/s = {:.3} us/read",
        best.as_secs_f64(),
        reads as f64 / best.as_secs_f64(),
        best.as_secs_f64() * 1_000_000.0 / reads as f64,
    );
    eprintln!(
        "  mean: {:.3}s = {:.0} reads/s = {:.3} us/read",
        mean_secs,
        reads as f64 / mean_secs,
        mean_secs * 1_000_000.0 / reads as f64,
    );
}

fn load_fastq_sequences(path: &Path, max_reads: usize) -> Result<Vec<Vec<u8>>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut lines = reader.lines();
    let mut reads = Vec::with_capacity(max_reads.min(1_000_000));

    while reads.len() < max_reads {
        let Some(header) = lines.next() else { break };
        let header = header?;
        anyhow::ensure!(header.starts_with('@'), "invalid FASTQ header: {header}");

        let seq = lines
            .next()
            .context("truncated FASTQ: missing sequence")??
            .into_bytes();
        let plus = lines.next().context("truncated FASTQ: missing + line")??;
        anyhow::ensure!(plus.starts_with('+'), "invalid FASTQ + line: {plus}");
        let qual = lines.next().context("truncated FASTQ: missing quality")??;
        anyhow::ensure!(
            seq.len() == qual.len(),
            "FASTQ sequence/quality length mismatch"
        );

        reads.push(seq);
    }

    Ok(reads)
}
