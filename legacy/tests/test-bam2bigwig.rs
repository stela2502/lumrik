
use std::process::Command;
use std::path::Path;
use std::fs;
#[allow(unused_imports)]
use std::process::exit;
use std::io::Read;

#[test]
fn test_bam2bigwig_creates_output() {
    let input_bam = "testData/test.bam";
    let output_bw = "testData/output/test.bw";

    // Remove output file if it exists (cleanup from prior runs)
    if Path::new(output_bw).exists() {
        fs::remove_file(output_bw).expect("Failed to remove old test output file");
    }

    let is_release_mode = !cfg!(debug_assertions);

    if ! is_release_mode {
        eprintln!("Test should be re-run in release mode (speed!)");
        return;
    }

    // Cargo sets this for integration tests when the binary exists in Cargo.toml [[bin]]
    let exe = env!("CARGO_BIN_EXE_bam2bigwig");

    // Run the compiled binary (assumes it's in target/debug or target/release)
    let status = Command::new(exe)
        .args([
            "-b", input_bam,
            "-o", output_bw,
            "--umi-tag", "No",
        ])
        .status()
        .expect("Failed to execute bam2bigwig");

    assert!(status.success(), "bam2bigwig did not exit successfully");

    // Check if the output file was created
    assert!(
        Path::new(output_bw).exists(),
        "Output file was not created: {}",
        output_bw
    );

    // BigWig is a UCSC binary format with a known magic in the first 4 bytes.
    // This is a fast sanity check that you're not writing bedGraph/wig text by mistake.
    let mut f = fs::File::open(output_bw).unwrap();
    let mut magic = [0u8; 4];
    f.read_exact(&mut magic).unwrap();

    // BigWig magic is typically 0x888FFC26 (endianness can differ when reading raw bytes)
    // Accept either byte order so the test isn't platform-endian brittle.
    let magic_be = u32::from_be_bytes(magic);
    let magic_le = u32::from_le_bytes(magic);
    assert!(
        magic_be == 0x888F_FC26 || magic_le == 0x888F_FC26,
        "Output does not look like BigWig (bad magic): bytes={:?} be=0x{:08X} le=0x{:08X}",
        magic, magic_be, magic_le
    );

    // Clean up the output file after test
    fs::remove_file(output_bw).expect("Failed to remove test output file after check");
}
