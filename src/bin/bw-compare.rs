use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::Parser;

use bam_tide::compare_report::CompareReport;

// bigtools
use bigtools::utils::reopen::ReopenableFile;
use bigtools::{BigWigRead, ChromInfo};

use std::fs::File;
use std::io::{BufWriter, Write};

type BwReader = BigWigRead<ReopenableFile>;

/// Compare two BigWig files (typically python bamCoverage vs. rust bam-coverage)
/// and report per-chromosome and global differences.
///
/// The tool bins both BigWigs with the same bin width and compares the values
/// position by position. It reports several statistics describing how different
/// the signals are.
///
/// Output columns:
///   n_over_eps       Number of bins where |python - rust| > eps
///   frac_n_over_eps  Fraction of bins over eps
///   mean_abs         Mean absolute difference
///   var_abs          Variance of absolute differences
///   rmse             Root mean squared error
///   max_abs          Maximum absolute difference
///   pearson_rho      Pearson correlation between tracks
///
/// A final TOTAL line summarizes all chromosomes.
///
/// If --outfile is not given, a report file will be created automatically:
///   bw_compare_<rust_basename>_w<bin_width>.txt
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Compare two BigWig files and report coverage differences"
)]
struct Args {
    /// Python-generated BigWig (reference)
    #[arg(long, short = 'a', value_name = "FILE")]
    python_bw: PathBuf,

    /// Rust-generated BigWig (to be validated)
    #[arg(long, short = 'b', value_name = "FILE")]
    rust_bw: PathBuf,

    /// Bin width used during coverage generation (must match both files)
    #[arg(long, default_value_t = 50, value_name = "INT")]
    bin_width: u32,

    /// Epsilon threshold for counting a bin as different
    /// (|python - rust| > eps)
    #[arg(long, default_value_t = 1e-5, value_name = "FLOAT")]
    eps: f64,

    /// Optional output file for the comparison report.
    ///
    /// If not provided, the report is written to STDOUT
    #[arg(long, short, value_name = "FILE")]
    outfile: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let mut py = BigWigRead::open_file(&args.python_bw)
        .with_context(|| format!("open python BW {}", args.python_bw.display()))?;

    let mut rs = BigWigRead::open_file(&args.rust_bw)
        .with_context(|| format!("open rust BW {}", args.rust_bw.display()))?;

    let mut out: Box<dyn Write> = if let Some(p) = args.outfile.clone() {
        Box::new(BufWriter::new(File::create(p)?))
    } else {
        // fallback: write nowhere
        Box::new(std::io::stdout())
    };

    // Read chrom lists
    let py_chroms: Vec<ChromInfo> = py.chroms().to_vec();
    let rs_chroms: Vec<ChromInfo> = rs.chroms().to_vec();

    // Build rust chrom lookup (borrowed, no clones)
    let mut rs_map: HashMap<String, u32> = HashMap::new();
    for c in &rs_chroms {
        rs_map.insert(c.name.clone(), c.length);
    }

    // Validate chrom sets + lengths
    for c in &py_chroms {
        match rs_map.get(c.name.as_str()) {
            None => bail!("chrom {} missing in rust BigWig", c.name),
            Some(&len2) if len2 != c.length => bail!(
                "chrom {} length mismatch: python={} rust={}",
                c.name,
                c.length,
                len2
            ),
            _ => {}
        }
    }

    let mut total = CompareReport::default();

    // Header
    let (n_over, frac, mean, var, rmse, max) = CompareReport::finish_names();
    writeln!(
        out,
        "CHR\t{}\t{}\t{}\t{}\t{}\t{}\tpearson_rho",
        n_over, frac, mean, var, rmse, max
    )?;

    for c in &py_chroms {
        let chr = c.name.as_str();
        let chr_len = c.length;

        let py_bins = bins_from_bigwig(&mut py, chr, chr_len, args.bin_width)?;
        let rs_bins = bins_from_bigwig(&mut rs, chr, chr_len, args.bin_width)?;

        let mut chr_rep = CompareReport::default();
        chr_rep.update_from_bins(
            chr,
            chr_len as u64,
            args.bin_width as u64,
            &py_bins,
            &rs_bins,
            args.eps,
        );

        // NEW finish()
        let (n_over_eps, frac_n_over_eps, mean_abs, var_abs, rmse, max_abs) = chr_rep.finish();

        let _ = writeln!(
            out,
            "{}\t{}\t{:.6}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.4}",
            chr,
            n_over_eps,
            frac_n_over_eps,
            mean_abs,
            var_abs,
            rmse,
            max_abs,
            chr_rep.pearson()
        );

        total.merge(&chr_rep);
    }
    let (n_over_eps, frac_n_over_eps, mean_abs, var_abs, rmse, max_abs) = total.finish();
    let _ = writeln!(
        out,
        "TOTAL\t{}\t{:.6}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.4}",
        n_over_eps,
        frac_n_over_eps,
        mean_abs,
        var_abs,
        rmse,
        max_abs,
        total.pearson()
    );

    println!("Most divergent bin:\n{:?}", total.max_bin);
    Ok(())
}

/// Create a binned representation of a BigWig chromosome
fn bins_from_bigwig(bw: &mut BwReader, chr: &str, chr_len: u32, bin_w: u32) -> Result<Vec<f64>> {
    let nbins = (chr_len as u64).div_ceil(bin_w as u64) as usize;
    let mut bins = vec![0.0_f64; nbins];

    let it = bw
        .get_interval(chr, 0, chr_len)
        .with_context(|| format!("intervals {chr}:0-{chr_len}"))?;

    for iv in it {
        let iv = iv.unwrap();
        let mut s = iv.start;
        let e = iv.end.min(chr_len);
        let v = iv.value as f64;

        if e <= s {
            continue;
        }

        while s < e {
            let bin_id = (s / bin_w) as usize;
            let bin_start = (bin_id as u32) * bin_w;
            let bin_end = (bin_start + bin_w).min(chr_len);

            let seg_end = e.min(bin_end);
            let overlap = (seg_end - s) as f64;

            bins[bin_id] += v * overlap;
            s = seg_end;
        }
    }

    // normalize to mean per base in bin
    for (i, x) in bins.iter_mut().enumerate() {
        let bin_start = (i as u32) * bin_w;
        let bin_end = (bin_start + bin_w).min(chr_len);
        let denom = (bin_end - bin_start) as f64;
        if denom > 0.0 {
            *x /= denom;
        }
    }

    Ok(bins)
}
