use assert_cmd::Command;
use predicates::prelude::*;

const SEQ: &str = "GTTGCCATTATAGTGAGTTGAATTCGACAATCACGCTTATTAAACGTGGAGTCGTGATTA";

#[test]
fn bd_cell_followed_by_fixed_adapter_still_matches() {
    let mut cmd = Command::cargo_bin("identify_primers").unwrap();

    cmd.args([
        "--primer-structure",
        "SEARCH:0..5+BD_CELL:v2.384+INSERT:CGTGGAGTCGTGATTA:mm2",
        "--seq",
        "AGTGGTTAGTGTGATTCTAATCGGACATGGTTCACTTTCGGACCGTGGAGTCGTGATTA",
    ]);

    cmd.assert().success().stdout(predicate::str::contains(
        "summary: 1 complete primer match(es)",
    ));
}

#[test]
fn bd_cell_followed_by_fixed_adapter_fails_if_sequences_fails() {
    let mut cmd = Command::cargo_bin("identify_primers").unwrap();

    cmd.args([
        "--primer-structure",
        "SEARCH:0..5+BD_CELL:v2.384+INSERT:CGTGGAGGGTGGATTA:mm2",
        "--seq",
        SEQ,
    ]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains(
            "summary: no complete primer match",
        ))
        .stdout(predicate::str::contains("status: OK").not())
        .stdout(predicate::str::contains("summary: 1 complete primer match(es)").not());
}
#[test]
fn identify_primers_reports_bd_whitelist_failure_reason() {
    let mut cmd = Command::cargo_bin("identify_primers").unwrap();

    cmd.args([
        "--chemistry",
        "bd-v2-384",
        "--seq",
        "GTTATTTTCGTGAATCTCAGAAGACATGTACAACGACCAGCCATTTTTTT",
    ]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("summary: no complete primer match"))
        .stdout(predicate::str::contains(
            "reason: BD_CELL: C1 is not exact and one-mismatch whitelist correction failed or was ambiguous",
        ));
}

#[test]
fn identify_primers_accepts_sequence_file_and_summarizes_errors() {
    let path = std::env::temp_dir().join(format!(
        "sc_primer_identify_{}_{}.txt",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));

    std::fs::write(
        &path,
        concat!(
            "TGCTGGCACGTGAATCTCAGAAGACATGTACAACGACCAGCCATTTTTTT\n",
            "GTTATTTTCGTGAATCTCAGAAGACATGTACAACGACCAGCCATTTTTTT\n",
        ),
    )
    .unwrap();

    let mut cmd = Command::cargo_bin("identify_primers").unwrap();
    cmd.args(["--chemistry", "bd-v2-384", "--seq", path.to_str().unwrap()]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("sequences: 2"))
        .stdout(predicate::str::contains("valid: 1"))
        .stdout(predicate::str::contains("invalid: 1"))
        .stdout(predicate::str::contains(
            "1\tBD_CELL: C1 is not exact and one-mismatch whitelist correction failed or was ambiguous",
        ));

    let _ = std::fs::remove_file(path);
}
