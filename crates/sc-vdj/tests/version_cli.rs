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
fn version_compare_to_airr_bd_0() {
    check_version("compare-to-airr-bd", env!("CARGO_BIN_EXE_compare-to-airr-bd"));
}

#[test]
fn version_create_test_vdj_index_1() {
    check_version("create-test-vdj-index", env!("CARGO_BIN_EXE_create-test-vdj-index"));
}

#[test]
fn version_vdj_decode_2() {
    check_version("vdj-decode", env!("CARGO_BIN_EXE_vdj-decode"));
}

#[test]
fn version_vdj_index_3() {
    check_version("vdj-index", env!("CARGO_BIN_EXE_vdj-index"));
}

#[test]
fn version_vdj_rich_cell_4() {
    check_version("vdj-rich-cell", env!("CARGO_BIN_EXE_vdj-rich-cell"));
}

#[test]
fn version_vdj_summary_5() {
    check_version("vdj-summary", env!("CARGO_BIN_EXE_vdj-summary"));
}

