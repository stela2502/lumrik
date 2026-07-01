use std::fs::File;
use std::path::PathBuf;
use std::process::Command;

//use tempfile::tempdir;

// bigtools read API
use bigtools::BigWigRead;

fn exe_path(name: &str) -> PathBuf {
    let base = if cfg!(debug_assertions) {
        PathBuf::from("./target/debug")
    } else {
        PathBuf::from("./target/release")
    };
    base.join(name)
}

fn run_ok(cmd: &mut Command) -> Result<(), String> {
    let out = cmd
        .output()
        .map_err(|e| format!("failed to spawn: {e:?}"))?;
    if !out.status.success() {
        return Err(format!(
            "command failed (exit={:?})\ncmd: {:?}\nstdout:\n{}\nstderr:\n{}",
            out.status.code(),
            cmd,
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    Ok(())
}

#[allow(dead_code)]
/// Compute weighted mean value in [start,end) from bigWig.
/// If no intervals overlap, returns 0.0 (bin is empty).
fn bw_weighted_mean(
    bw: &mut BigWigRead<&File>,
    chr: &str,
    start: u32,
    end: u32,
) -> Result<f32, String> {
    let mut it = bw
        .get_interval(chr, start, end)
        .map_err(|e| format!("bigwig get_interval failed for {chr}:{start}-{end}: {e:?}"))?;

    let mut sum = 0.0_f64;
    let mut covered = 0u64;

    while let Some(item) = it.next() {
        let v = item.map_err(|e| format!("bigwig interval item error: {e:?}"))?;
        // overlap with requested window
        let s = start.max(v.start);
        let e = end.min(v.end);
        if e > s {
            let w = (e - s) as u64;
            covered += w;
            sum += (v.value as f64) * (w as f64);
        }
    }

    if covered == 0 {
        Ok(0.0)
    } else {
        Ok((sum / covered as f64) as f32)
    }
}

#[allow(dead_code)]
/// Pearson correlation over paired vectors (ignores all-zero case gracefully).
fn pearson_corr(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }
    let mean_x = x.iter().sum::<f64>() / n as f64;
    let mean_y = y.iter().sum::<f64>() / n as f64;

    let mut num = 0.0;
    let mut den_x = 0.0;
    let mut den_y = 0.0;
    for i in 0..n {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }
    if den_x == 0.0 || den_y == 0.0 {
        return f64::NAN;
    }
    num / (den_x.sqrt() * den_y.sqrt())
}

#[test]
fn rust_bigwig_matches_golden_deeptools_bwcompare() -> Result<(), String> {
    use std::fs;
    use std::path::Path;
    use std::process::Command;
    use tempfile::tempdir;

    if cfg!(debug_assertions) {
        eprintln!("Skipping slow bigWig comparison test in debug mode");
        return Ok(());
    }

    // --- config ---
    let bam = Path::new("legacy/testData/subset.bam");
    let golden = Path::new("legacy/testData/subset.deepTools.bw");
    let bin_width: u32 = 50;

    // thresholds
    let eps_abs: f64 = 0.0;
    let max_frac_bad: f64 = 0.0;
    let min_corr: f64 = 1.0;

    if !bam.exists() {
        return Err(format!("missing test BAM: {}", bam.display()));
    }
    if !golden.exists() {
        return Err(format!("missing golden bigWig: {}", golden.display()));
    }

    let tmp = tempdir().map_err(|e| format!("tempdir: {e:?}"))?;
    let out_bw = tmp.path().join("rust.bw");
    let out_tsv = Path::new("legacy/testData/cmp.tsv");

    // --- produce rust bigWig ---
    let exe_cov = exe_path("bam-coverage");
    if !exe_cov.exists() {
        return Err(format!("missing binary: {}", exe_cov.display()));
    }

    run_ok(Command::new(&exe_cov).args([
        "-b",
        bam.to_str().unwrap(),
        "-o",
        out_bw.to_str().unwrap(),
        "-w",
        &bin_width.to_string(),
        "-n",
        "not",
        "--min-mapping-quality",
        "0",
    ]))?;

    if !out_bw.exists()
        || fs::metadata(&out_bw)
            .map_err(|e| format!("stat rust.bw: {e:?}"))?
            .len()
            == 0
    {
        return Err("rust bigWig not created or empty".into());
    }

    // --- compare with bw-compare ---
    let exe_cmp = exe_path("bw-compare");
    if !exe_cmp.exists() {
        return Err(format!("missing binary: {}", exe_cmp.display()));
    }

    // IMPORTANT: adapt args to your bw-compare CLI.
    // Goal: emit a TSV with a TOTAL row containing fields:
    //   n_over_eps, frac_n_over_eps, mean_abs, rmse, max_abs, pearson_rho
    run_ok(Command::new(&exe_cmp).args([
        "-a",
        golden.to_str().unwrap(),
        "-b",
        out_bw.to_str().unwrap(),
        "--eps",
        &eps_abs.to_string(),
        "-o",
        out_tsv.to_str().unwrap(),
    ]))?;

    let tsv = fs::read_to_string(&out_tsv).map_err(|e| format!("read cmp.tsv: {e:?}"))?;

    // --- parse TOTAL row ---
    // Expected header like:
    // flag  CHR  n_over_eps  frac_n_over_eps  mean_abs  var_abs  rmse  max_abs  pearson_rho ...
    let mut total_line: Option<&str> = None;
    for line in tsv.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        // Accept either "TOTAL" in CHR col, or line starting with "TOTAL"
        if line.contains("TOTAL\t") || line.starts_with("TOTAL\t") {
            total_line = Some(line);
            break;
        }
    }

    let cmp_cmdline = format!(
        "{} -a {} -b {}",
        exe_cmp.display(),
        golden.display(),
        out_bw.display()
    );

    let line = total_line
        .ok_or_else(|| format!("no TOTAL row found in bw-compare output:\n{cmp_cmdline}\n{tsv}"))?;

    let cols: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();

    // Try to locate columns by header if available; otherwise assume fixed indices.
    // Here we assume:
    // 0 flag, 1 CHR, 2 n_over_eps, 3 frac_n_over_eps, 4 mean_abs, 6 rmse, 7 max_abs, 8 pearson_rho
    let frac_bad: f64 = cols
        .get(2)
        .ok_or("missing frac_n_over_eps")?
        .parse()
        .map_err(|_| format!("bad frac_n_over_eps in: {line}"))?;
    let mean_abs: f64 = cols
        .get(3)
        .ok_or("missing mean_abs")?
        .parse()
        .map_err(|_| format!("bad mean_abs in: {line}"))?;
    let max_abs: f64 = cols
        .get(6)
        .ok_or("missing max_abs")?
        .parse()
        .map_err(|_| format!("bad max_abs in: {line}"))?;
    let corr: f64 = cols
        .get(7)
        .ok_or_else(|| format!("missing pearson_rho; line={line}"))?
        .parse()
        .map_err(|_| format!("bad pearson_rho in: {line}"))?;

    eprintln!(
        "bw-compare TOTAL: frac_bad={:.3e} mean_abs={:.6} max_abs={:.6} corr={:.6}",
        frac_bad, mean_abs, max_abs, corr
    );

    if frac_bad > max_frac_bad {
        return Err(format!(
            "frac_bad too high: {frac_bad:.3e} > {max_frac_bad:.3e}"
        ));
    }
    // optional: mean_abs/max_abs checks, if you want
    if corr.is_nan() || corr < min_corr {
        return Err(format!("corr too low: {corr} < {min_corr}"));
    }

    Ok(())
}
