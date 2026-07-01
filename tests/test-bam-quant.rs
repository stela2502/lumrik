use assert_cmd::Command;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command as StdCommand;

const PBMC_URL: &str = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam";

// GENCODE Human Release 49, “Basic gene annotation CHR” GTF on EMBL-EBI FTP. :contentReference[oaicite:1]{index=1}
const GENCODE_BASIC_CHR_GTF_GZ: &str = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz";

const GRCH38_FA_URL: &str =
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz";

fn run(cmd: &mut StdCommand) {
    let status = cmd.status().expect("failed to run command");
    assert!(status.success(), "command failed: {:?}", cmd);
}

fn download_if_needed(url: &str, path: &Path) {
    if path.exists() {
        eprintln!("✓ cached: {}", path.display());
        return;
    }
    eprintln!("⬇ downloading: {}", url);
    run(StdCommand::new("wget").args(["-O", path.to_str().unwrap(), url]));
}

fn require_exe(name: &str) {
    let exists = std::process::Command::new("which")
        .arg(name)
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    assert!(exists, "required executable not found in PATH: {}", name);
}

fn require_tools() {
    require_exe("gtf_splice_index");
    require_exe("samtools");
    require_exe("wget");
    require_exe("bash");
    require_exe("gzip");
    require_exe("awk");
    require_exe("bcftools");
}

/// Create chr1-only GTF (keeps all header lines, filters feature lines by seqname == "chr1").
/// Works for typical GENCODE GTFs where seqname is column 1.
fn subset_gtf_chr1_if_needed(gtf_gz: &Path, out_gtf: &Path) {
    if out_gtf.exists() {
        eprintln!("✓ cached: {}", out_gtf.display());
        return;
    }

    eprintln!("⚙ subsetting chr1 GTF: {}", out_gtf.display());

    // gzip -dc <in> | awk '(^#) || ($1=="chr1")' > out
    let script = format!(
        r#"
set -euo pipefail
gzip -dc "{in_gz}" \
  | awk 'BEGIN{{FS="\t"}} /^#/ {{print; next}} ($1=="chr1") {{print}}' \
  > "{out_gtf}"
"#,
        in_gz = gtf_gz.display(),
        out_gtf = out_gtf.display(),
    );

    run(StdCommand::new("bash").arg("-c").arg(script));
}

/// Build your splice index if missing.
/// Assumes you have a workspace binary named `gtf_splice_index` with a `build` subcommand:
///   gtf_splice_index build -a <annotation.gtf> -i <index.dat>
fn build_index_if_needed(gtf_chr1: &Path, index_path: &Path) {
    if index_path.exists() {
        eprintln!("✓ cached: {}", index_path.display());
        return;
    }

    eprintln!("⚙ building index: {}", index_path.display());

    let mut cmd = assert_cmd::Command::new("gtf_splice_index");

    cmd.args([
        "build",
        "-a",
        gtf_chr1.to_str().unwrap(),
        "-i",
        index_path.to_str().unwrap(),
    ]);

    cmd.assert().success();
}

fn ensure_clean_dir(outdir: &Path) -> std::io::Result<()> {
    if outdir.exists() {
        fs::remove_dir_all(outdir)?; // removes dir + all contents
    }
    fs::create_dir_all(outdir)?;
    Ok(())
}

fn add_suffix(path: &Path, suffix: &str) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let name = path
        .file_name()
        .expect("path has no filename")
        .to_string_lossy();

    parent.join(format!("{name}{suffix}"))
}

fn get_matrix_nnz(matrix_gz: &Path) -> usize {
    let output = StdCommand::new("bash")
        .arg("-c")
        .arg(format!(
            r#"gzip -dc "{}" | awk 'BEGIN{{n=0}} !/^%/ {{n++; if (n==2) {{print $3; exit}}}}'"#,
            matrix_gz.display()
        ))
        .output()
        .expect("failed to inspect matrix");

    assert!(
        output.status.success(),
        "failed reading matrix {}: {}",
        matrix_gz.display(),
        String::from_utf8_lossy(&output.stderr)
    );

    String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .unwrap_or_else(|_| {
            panic!(
                "could not parse matrix nnz from {}. stdout was: {:?}",
                matrix_gz.display(),
                String::from_utf8_lossy(&output.stdout)
            )
        })
}

#[test]
#[ignore = "Downloads data + requires wget/gzip/awk/samtools. Run manually (~10min)."]
fn bam_quant_manual_real_data_end_to_end() {
    // ---- knobs ----
    let n_reads: u64 = std::env::var("BAM_QUANT_N")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(2_000_000);

    require_tools(); // check for all external tools

    // ---- persistent cache ----
    let cache_dir = PathBuf::from("/tmp/bam_quant_cache");
    std::fs::create_dir_all(&cache_dir).unwrap();

    // Inputs
    let pbmc_bam = cache_dir.join("pbmc_1k.bam");
    let gtf_gz = cache_dir.join("gencode.basic.chr.gtf.gz");
    let gtf_chr1 = cache_dir.join("gencode.basic.chr1.gtf");

    // SNP analysis
    let genome_fa_gz = cache_dir.join("genome.fa.gz");
    let vcf_chr1 = PathBuf::from("tests/data/chr1_observed_snps.vcf.gz");

    // Derived artifacts
    let _subset_bam = cache_dir.join(format!("pbmc_subset_{n_reads}.bam"));
    let index_path = cache_dir.join("splice_index_chr1.dat");

    // Output
    let outdir = cache_dir.join(format!("out_{n_reads}_chr1"));
    std::fs::create_dir_all(&outdir).unwrap();

    // ---- step 1: download inputs ----
    download_if_needed(PBMC_URL, &pbmc_bam);
    download_if_needed(GENCODE_BASIC_CHR_GTF_GZ, &gtf_gz);

    let enable_snps = cfg!(feature = "manual-snp-tests");

    if enable_snps {
        download_if_needed(GRCH38_FA_URL, &genome_fa_gz);
        //download_if_needed(DBSNP_VCF_URL, &vcf_gz);
    }
    // ---- step 2: subset GTF to chr1 ----
    subset_gtf_chr1_if_needed(&gtf_gz, &gtf_chr1);

    // ---- step 3: build index (via gtf_splice_index) ----
    build_index_if_needed(&gtf_chr1, &index_path);

    // ---- step 5: run bam-quant - clean outdirs if necessary ----
    let _ = ensure_clean_dir(&outdir);
    let outdir2 = cache_dir.join(format!("out_{n_reads}_chr1_intronic"));
    let _ = ensure_clean_dir(&outdir2);

    eprintln!("🚀 running bam-quant…");

    let mut cmd = Command::cargo_bin("bam-quant").unwrap();
    let mut args = vec![
        "--bam".to_string(),
        pbmc_bam.to_string_lossy().to_string(),
        "--index".to_string(),
        index_path.to_string_lossy().to_string(),
        "--outpath".to_string(),
        outdir.to_string_lossy().to_string(),
        "--max-reads".to_string(),
        n_reads.to_string(),
        "--min-cell-counts".to_string(),
        "10".to_string(),
    ];

    if enable_snps {
        args.extend([
            "--genome".to_string(),
            genome_fa_gz.to_string_lossy().to_string(),
            "--vcf".to_string(),
            vcf_chr1.to_string_lossy().to_string(),
        ]);
    }

    // Apply to command
    cmd.args(&args);
    // Print runnable command
    eprintln!(
        "\n🔥 Run manually:\n{}\n",
        format!("bam-quant {}", shell_escape_args(&args))
    );

    cmd.assert().success();

    assert_10x_output_nonempty(&outdir);
    assert_matrix_has_entries(&outdir.join("matrix.mtx.gz"));

    if enable_snps {
        let ref_out = add_suffix(&outdir, "_ref");
        let alt_out = add_suffix(&outdir, "_alt");

        eprintln!("🔬 checking SNP outputs:");
        eprintln!("  ref: {}", ref_out.display());
        eprintln!("  alt: {}", alt_out.display());

        // ---- existence ----
        assert!(ref_out.exists(), "missing SNP ref output dir");
        assert!(alt_out.exists(), "missing SNP alt output dir");

        // ---- basic 10x structure ----
        assert_10x_output_nonempty(&ref_out);
        assert_10x_output_nonempty(&alt_out);

        // ---- matrix sanity ----
        let ref_matrix = ref_out.join("matrix.mtx.gz");
        let alt_matrix = alt_out.join("matrix.mtx.gz");

        // ensure files exist + non-empty
        assert_file_nonempty(&ref_matrix);
        assert_file_nonempty(&alt_matrix);

        // ---- inspect nnz (non-zero entries) ----
        let nnz_ref = get_matrix_nnz(&ref_matrix);
        let nnz_alt = get_matrix_nnz(&alt_matrix);

        eprintln!("  SNP nnz: ref={nnz_ref}, alt={nnz_alt}");

        // At least one of them should have signal
        assert!(
            nnz_ref > 0 || nnz_alt > 0,
            "both SNP matrices are empty (no signal detected)"
        );

        // ---- optional stronger check ----
        // If both are non-zero, ensure they are not identical (sanity)
        if nnz_ref > 0 && nnz_alt > 0 {
            assert!(
                nnz_ref != nnz_alt,
                "ref and alt matrices have identical nnz — suspicious"
            );
        }
    }
}

fn assert_file_nonempty(path: &Path) {
    let meta =
        std::fs::metadata(path).unwrap_or_else(|_| panic!("missing file: {}", path.display()));
    assert!(meta.len() > 0, "empty file: {}", path.display());
}

fn assert_10x_output_nonempty(outdir: &Path) {
    assert_file_nonempty(&outdir.join("matrix.mtx.gz"));
    assert_file_nonempty(&outdir.join("barcodes.tsv.gz"));
    assert_file_nonempty(&outdir.join("features.tsv.gz"));
}

fn assert_matrix_has_entries(matrix_gz: &Path) {
    let output = StdCommand::new("bash")
        .arg("-c")
        .arg(format!(
            r#"gzip -dc "{}" | awk 'BEGIN{{n=0}} !/^%/ {{n++; if (n==2) {{print $3; exit}}}}'"#,
            matrix_gz.display()
        ))
        .output()
        .expect("failed to inspect matrix");

    assert!(output.status.success());

    let nnz: usize = String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .expect("could not parse matrix nnz");

    assert!(nnz > 0, "matrix has zero non-zero entries");
}

fn shell_escape_args(args: &[String]) -> String {
    args.iter()
        .map(|a| {
            if a.contains(' ') {
                format!("'{}'", a.replace('\'', "'\\''"))
            } else {
                a.clone()
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}
