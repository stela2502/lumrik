use std::path::{Path, PathBuf};
use std::process::Command;

#[test]
fn test_bam_coverage_matches_deeptools_bigwig() {
    let bam = "legacy/testData/subset.bam";

    if cfg!(debug_assertions) {
        eprintln!("Skipping slow bigWig comparison test in debug mode");
        return;
    }

    let out_dir = PathBuf::from("legacy/testData");
    let py_bw = out_dir.join("_ref_deeptools.bw");
    let rs_bw = out_dir.join("_test_bam_coverage.bw");

    let _eps_abs: f64 = 0.0;
    let max_frac_bad: f64 = 0.0;
    let min_corr: f64 = 1.0;

    // ------------------------------------------------------------------
    // 1) Create reference BigWig with deeptools bamCoverage
    // ------------------------------------------------------------------
    if !Path::new(&py_bw).exists() {
        let status = Command::new("bamCoverage")
            .args([
                "-b",
                bam,
                "-o",
                py_bw.to_str().unwrap(),
                "--binSize",
                "50",
                "--normalizeUsing",
                "None",
                "--minMappingQuality",
                "0",
                //"--ignoreDuplicates",
                "--extendReads",
                "0",
                "--numberOfProcessors",
                "4",
            ])
            .status()
            .expect("failed to run deeptools bamCoverage");

        assert!(status.success(), "deeptools bamCoverage failed");
    }

    // ------------------------------------------------------------------
    // 2) Create BigWig with our bam-coverage
    // ------------------------------------------------------------------
    let exe = env!("CARGO_BIN_EXE_bam-coverage");

    let status = Command::new(exe)
        .args([
            "-b",
            bam,
            "-o",
            rs_bw.to_str().unwrap(),
            "-w",
            "50",
            "--min-mapping-quality",
            "0",
        ])
        .status()
        .expect("failed to run bam-coverage");

    assert!(status.success(), "bam-coverage failed");

    // ------------------------------------------------------------------
    // 3) Compare using bw-compare
    // ------------------------------------------------------------------
    let cmp = env!("CARGO_BIN_EXE_bw-compare");

    let output = Command::new(cmp)
        .args([
            "--python-bw",
            py_bw.to_str().unwrap(),
            "--rust-bw",
            rs_bw.to_str().unwrap(),
            "--bin-width",
            "50",
            "--eps",
            "1e-4",
        ])
        .output()
        .expect("failed to run bw-compare");

    if !output.status.success() {
        panic!(
            "bw-compare failed\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }
    let stdout = String::from_utf8_lossy(&output.stdout);
    let cmp_cmdline = format!(
        "{} -a {} -b {} --bin-width 50 --eps 1e-4",
        cmp,
        py_bw.display(),
        rs_bw.display()
    );
    match parse_total_line(&stdout) {
        Ok((frac_bad, _mean_abs, _max_abs, corr)) => {
            if frac_bad > max_frac_bad {
                panic!(
                    "{}",
                    format!("frac_bad too high: {frac_bad:.3e} > {max_frac_bad:.3e}")
                );
            }
            // optional: mean_abs/max_abs checks, if you want
            if corr.is_nan() || corr < min_corr {
                panic!("{}", format!("corr too low: {corr} < {min_corr}"));
            }
        }
        Err(e) => panic!(
            "{}",
            format!("{cmp_cmdline} create the wrong output {e:?} ")
        ),
    };
}

fn parse_total_line(stdout: &str) -> Result<(f64, f64, f64, f64), String> {
    for line in stdout.lines() {
        if line.contains("TOTAL\t") {
            let cols: Vec<&str> = line.split('\t').collect();
            assert!(cols.len() >= 5, "Malformed TOTAL line: {line}");

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

            return Ok((frac_bad, mean_abs, max_abs, corr));
        }
    }
    Err(format!(
        "No TOTAL line found in bw-compare output:\n{stdout}"
    ))
}
