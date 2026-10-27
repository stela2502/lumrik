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
fn version_bd_fast_map_0() {
    check_version("bd-fast-map", env!("CARGO_BIN_EXE_bd-fast-map"));
}

#[test]
fn version_fast_tag_stress_1() {
    check_version("fast-tag-stress", env!("CARGO_BIN_EXE_fast-tag-stress"));
}

