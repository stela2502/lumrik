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
fn identify_primers_reports_rescued_bd_cell() {
    let mut cmd = Command::cargo_bin("identify_primers").unwrap();

    cmd.args([
        "--chemistry",
        "bd-v2-384",
        "--seq",
        "GTTATTTTCGTGAATCTCAGAAGACATGTACAACGACCAGCCATTTTTTT",
    ]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("summary: 1 complete primer match(es)"))
        .stdout(predicate::str::contains("cell_seq: GTTAATTCCATCTCAGAATGTACAACG"));
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
            "ANAGGAAACTCATGGTGCGTGGATCTGGCAATGAGCCTGCCGCCACTATCAGTCGTGGCATATGTGAGTCGTGATTATAGAGAGAGAGACCAAAATTCAAAGAGAAAATGGATTTTCAGGTGCAGATTTTCAGCTTCCTGCTAATCAGTGC\n",
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
            "1\tBD_CELL: combined 8 bp linker has more than two mismatches",
        ));

    let _ = std::fs::remove_file(path);
}

#[test]
fn identify_primers_rejects_mouse_igk_transcript() {
    const MOUSE_IGK: &str = "ANAGGAAACTCATGGTGCGTGGATCTGGCAATGAGCCTGCCGCCACTATCAGTCGTGGCATATGTGAGTCGTGATTATAGAGAGAGAGACCAAAATTCAAAGAGAAAATGGATTTTCAGGTGCAGATTTTCAGCTTCCTGCTAATCAGTGC";

    let mut cmd = Command::cargo_bin("identify_primers").unwrap();
    cmd.args(["--chemistry", "bd-v2-384", "--seq", MOUSE_IGK]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("summary: no complete primer match"))
        .stdout(predicate::str::contains(
            "reason: BD_CELL: combined 8 bp linker has more than two mismatches",
        ))
        .stdout(predicate::str::contains("orientation: ReverseComplement").not());
}
