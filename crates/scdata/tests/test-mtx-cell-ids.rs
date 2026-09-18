use std::fs::{self, File};
use std::io::Write;

use flate2::Compression;
use flate2::write::GzEncoder;
use int_to_str::IntToStr;
use scdata::read_mtx_cell_ids;
use tempfile::tempdir;

#[test]
fn reads_cell_ids_from_plain_mex_barcodes() {
    let tmp = tempdir().unwrap();
    let mex = tmp.path().join("mex");
    fs::create_dir_all(&mex).unwrap();
    fs::write(mex.join("barcodes.tsv"), "ACGT\nTGCA\n").unwrap();

    let cells = read_mtx_cell_ids(&mex).unwrap();
    assert_eq!(cells.len(), 2);
    assert!(cells.contains(&IntToStr::new(b"ACGT").into_u64()));
    assert!(cells.contains(&IntToStr::new(b"TGCA").into_u64()));
}

#[test]
fn reads_cell_ids_from_nelrune_analysis_exonic_directory() {
    let tmp = tempdir().unwrap();
    let exonic = tmp.path().join("exonic");
    fs::create_dir_all(&exonic).unwrap();

    let file = File::create(exonic.join("barcodes.tsv.gz")).unwrap();
    let mut gz = GzEncoder::new(file, Compression::default());
    writeln!(gz, "AACCGGTT").unwrap();
    writeln!(gz, "TTGGCCAA").unwrap();
    gz.finish().unwrap();

    let cells = read_mtx_cell_ids(tmp.path()).unwrap();
    assert_eq!(cells.len(), 2);
    assert!(cells.contains(&IntToStr::new(b"AACCGGTT").into_u64()));
    assert!(cells.contains(&IntToStr::new(b"TTGGCCAA").into_u64()));
}

#[test]
fn missing_mex_barcodes_is_an_error() {
    let tmp = tempdir().unwrap();
    let err = read_mtx_cell_ids(tmp.path()).unwrap_err().to_string();
    assert!(err.contains("neither a MEX directory"));
}
