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
fn version_clonomap_0() {
    check_version("clonomap", env!("CARGO_BIN_EXE_clonomap"));
}

#[test]
fn version_clonomap_batch_1() {
    check_version("clonomap_batch", env!("CARGO_BIN_EXE_clonomap_batch"));
}

#[test]
fn version_tree_viewer_2() {
    check_version("tree_viewer", env!("CARGO_BIN_EXE_tree_viewer"));
}

