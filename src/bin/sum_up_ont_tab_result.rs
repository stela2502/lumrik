use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
};

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Input TSV, optionally .gz
    input: PathBuf,

    /// Minimum number of observations for one raw_cb + raw_umi pair
    #[arg(long, default_value_t = 2)]
    min_pair_count: u64,

    /// Minimum number of surviving UMIs per raw_cb after pair filtering
    #[arg(long, default_value_t = 3)]
    min_cb_umis: u64,
}

#[derive(Debug)]
struct Summary {
    n: usize,
    mean: f64,
    median: f64,
    min: u64,
    max: u64,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let reader = open_maybe_gz(&cli.input)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(true)
        .from_reader(reader);

    let headers = rdr.headers()?.clone();

    let raw_cb_ix = headers
        .iter()
        .position(|h| h == "raw_cb")
        .context("Could not find column 'raw_cb'")?;

    let raw_umi_ix = headers
        .iter()
        .position(|h| h == "raw_umi")
        .context("Could not find column 'raw_umi'")?;

    let mut cb_umi_counts: HashMap<(String, String), u64> = HashMap::new();

    let mut total_rows = 0u64;
    let mut used_rows = 0u64;
    let mut skipped_missing = 0u64;
    let mut skipped_empty = 0u64;

    for rec in rdr.records() {
        total_rows += 1;

        let rec = rec?;

        let raw_cb = match rec.get(raw_cb_ix).map(str::trim) {
            Some(v) => v,
            None => {
                skipped_missing += 1;
                continue;
            }
        };

        let raw_umi = match rec.get(raw_umi_ix).map(str::trim) {
            Some(v) => v,
            None => {
                skipped_missing += 1;
                continue;
            }
        };

        if raw_cb.is_empty() || raw_umi.is_empty() {
            skipped_empty += 1;
            continue;
        }

        used_rows += 1;

        *cb_umi_counts
            .entry((raw_cb.to_owned(), raw_umi.to_owned()))
            .or_insert(0) += 1;
    }

    let pre = build_stats(&cb_umi_counts, 1, 1);

    let post = build_stats(
        &cb_umi_counts,
        cli.min_pair_count,
        cli.min_cb_umis,
    );

    println!("Input: {}", cli.input.display());
    println!();
    println!("Rows");
    println!("  total           : {total_rows}");
    println!("  used            : {used_rows}");
    println!("  skipped missing : {skipped_missing}");
    println!("  skipped empty   : {skipped_empty}");
    println!();

    println!("Filtering");
    println!("  min_pair_count  : {}", cli.min_pair_count);
    println!("  min_cb_umis     : {}", cli.min_cb_umis);
    println!();

    print_stats("Pre-filtered", &pre);
    println!();
    print_stats("Post-filtered", &post);

    Ok(())
}

#[derive(Debug)]
struct Stats {
    raw_cb_entries: usize,
    unique_cb_umi_combos: usize,
    total_pair_observations: u64,
    umis_per_cb: Summary,
    detections_per_cb_umi: Summary,
}

fn build_stats(
    cb_umi_counts: &HashMap<(String, String), u64>,
    min_pair_count: u64,
    min_cb_umis: u64,
) -> Stats {
    let mut cb_to_umis: HashMap<&str, HashSet<&str>> = HashMap::new();
    let mut pair_counts = Vec::new();

    for ((cb, umi), count) in cb_umi_counts {
        if *count < min_pair_count {
            continue;
        }

        cb_to_umis
            .entry(cb.as_str())
            .or_default()
            .insert(umi.as_str());

        pair_counts.push(*count);
    }

    let valid_cbs: HashSet<&str> = cb_to_umis
        .iter()
        .filter_map(|(cb, umis)| {
            if umis.len() as u64 >= min_cb_umis {
                Some(*cb)
            } else {
                None
            }
        })
        .collect();

    let mut final_cb_to_umi_count: HashMap<&str, u64> = HashMap::new();
    let mut final_pair_counts = Vec::new();
    let mut total_pair_observations = 0u64;

    for ((cb, _umi), count) in cb_umi_counts {
        if *count < min_pair_count {
            continue;
        }

        if !valid_cbs.contains(cb.as_str()) {
            continue;
        }

        *final_cb_to_umi_count.entry(cb.as_str()).or_insert(0) += 1;
        final_pair_counts.push(*count);
        total_pair_observations += *count;
    }

    let umis_per_cb: Vec<u64> = final_cb_to_umi_count.values().copied().collect();

    Stats {
        raw_cb_entries: final_cb_to_umi_count.len(),
        unique_cb_umi_combos: final_pair_counts.len(),
        total_pair_observations,
        umis_per_cb: summarize(umis_per_cb),
        detections_per_cb_umi: summarize(final_pair_counts),
    }
}

fn print_stats(name: &str, stats: &Stats) {
    println!("{name}");
    println!("  Raw CB entries              : {}", stats.raw_cb_entries);
    println!(
        "  Unique raw_cb + raw_umi     : {}",
        stats.unique_cb_umi_combos
    );
    println!(
        "  Total pair observations     : {}",
        stats.total_pair_observations
    );
    println!();

    print_summary("  UMIs per raw_cb", &stats.umis_per_cb);
    println!();
    print_summary(
        "  Detections per raw_cb + raw_umi",
        &stats.detections_per_cb_umi,
    );
}

fn print_summary(name: &str, s: &Summary) {
    println!("{name}");
    println!("    n      : {}", s.n);
    println!("    mean   : {:.3}", s.mean);
    println!("    median : {:.3}", s.median);
    println!("    min    : {}", s.min);
    println!("    max    : {}", s.max);
}

fn summarize(mut values: Vec<u64>) -> Summary {
    if values.is_empty() {
        return Summary {
            n: 0,
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
        };
    }

    values.sort_unstable();

    let n = values.len();
    let min = values[0];
    let max = values[n - 1];

    let sum: u128 = values.iter().map(|&x| x as u128).sum();
    let mean = sum as f64 / n as f64;

    let median = if n % 2 == 0 {
        let a = values[n / 2 - 1];
        let b = values[n / 2];
        (a as f64 + b as f64) / 2.0
    } else {
        values[n / 2] as f64
    };

    Summary {
        n,
        mean,
        median,
        min,
        max,
    }
}

fn open_maybe_gz(path: &PathBuf) -> Result<Box<dyn Read>> {
    let file = File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;

    let reader = BufReader::new(file);

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(reader)))
    } else {
        Ok(Box::new(reader))
    }
}
