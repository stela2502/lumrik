use assert_cmd::Command;
use predicates::prelude::*;

const SEQ: &str =
    "GTTGCCATTATAGTGAGTTGAATTCGACAATCACGCTTATTAAACGTGGAGTCGTGATTA";

#[test]
fn bd_cell_followed_by_fixed_adapter_still_matches() {
    let mut cmd = Command::cargo_bin("identify_primers").unwrap();

    cmd.args([
        "--primer-structure",
        "SEARCH:0..5+BD_CELL:v2.384+INSERT:CGTGGAGTCGTGATTA:mm2",
        "--seq",
        "AGTGGTTAGTGTGATTCTAATCGGACATGGTTCACTTTCGGACCGTGGAGTCGTGATTA",
    ]);

    cmd.assert()
        .success()
        .stdout(predicate::str::contains(
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
        .stdout(predicate::str::contains("summary: no complete primer match"))
        .stdout(predicate::str::contains("status: OK").not())
        .stdout(predicate::str::contains("summary: 1 complete primer match(es)").not());
}