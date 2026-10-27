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
fn version_bd_linker_survey_0() {
    check_version("bd-linker-survey", env!("CARGO_BIN_EXE_bd-linker-survey"));
}

#[test]
fn version_identify_primers_1() {
    check_version("identify_primers", env!("CARGO_BIN_EXE_identify_primers"));
}

#[test]
fn version_primer_diagnostics_2() {
    check_version("primer-diagnostics", env!("CARGO_BIN_EXE_primer-diagnostics"));
}

#[test]
fn version_primer_stress_3() {
    check_version("primer-stress", env!("CARGO_BIN_EXE_primer-stress"));
}

#[test]
fn version_sc_cell_id_convert_4() {
    check_version("sc-cell-id-convert", env!("CARGO_BIN_EXE_sc-cell-id-convert"));
}

