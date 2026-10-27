// Verify that every executable in this crate reports its Cargo package version.
// These are integration tests: Cargo builds the binaries and supplies their paths.
use std::process::Command;

fn check_version(binary: &str, path: &str) {
    let output = Command::new(path)
        .arg("--version")
        .output()
        .unwrap_or_else(|e| panic!("could not execute {binary}: {e}"));
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "{binary} --version failed ({}): stdout={stdout:?}, stderr={stderr:?}",
        output.status
    );
    let version = env!("CARGO_PKG_VERSION");
    assert!(
        stdout.contains(version),
        "{binary} --version should contain Cargo package version {version:?}; got {stdout:?}"
    );
}

#[test]
fn version_bam_coverage_0() {
    check_version("bam-coverage", env!("CARGO_BIN_EXE_bam-coverage"));
}

#[test]
fn version_bam_quant_1() {
    check_version("bam-quant", env!("CARGO_BIN_EXE_bam-quant"));
}

#[test]
fn version_bam_read_tag_stats_2() {
    check_version("bam-read-tag-stats", env!("CARGO_BIN_EXE_bam-read-tag-stats"));
}

#[test]
fn version_bam_subset_tag_3() {
    check_version("bam-subset-tag", env!("CARGO_BIN_EXE_bam-subset-tag"));
}

#[test]
fn version_bam_transcript_prob_4() {
    check_version("bam-transcript-prob", env!("CARGO_BIN_EXE_bam-transcript-prob"));
}

#[test]
fn version_bam_transcriptome_to_genome_5() {
    check_version("bam-transcriptome-to-genome", env!("CARGO_BIN_EXE_bam-transcriptome-to-genome"));
}

#[test]
fn version_bw_compare_6() {
    check_version("bw-compare", env!("CARGO_BIN_EXE_bw-compare"));
}

#[test]
fn version_fastq_primer_restore_7() {
    check_version("fastq-primer-restore", env!("CARGO_BIN_EXE_fastq-primer-restore"));
}

#[test]
fn version_find_sequences_8() {
    check_version("find-sequences", env!("CARGO_BIN_EXE_find-sequences"));
}

#[test]
fn version_gtf_splice_index_9() {
    check_version("gtf-splice-index", env!("CARGO_BIN_EXE_gtf-splice-index"));
}

#[test]
fn version_sum_up_ont_tab_result_10() {
    check_version("sum_up_ont_tab_result", env!("CARGO_BIN_EXE_sum_up_ont_tab_result"));
}

