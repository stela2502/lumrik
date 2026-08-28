use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::MultiGzDecoder;
use gtf_splice_index::{AnnotationBuilder, SpliceIndex};

struct Cleanup(PathBuf);

impl Drop for Cleanup {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

const TEST_DATA: &str = "tests/data";

const BD_MOUSE_SAMPLE_07_READ: &str =
    "GTTGTCAAGATGCTACCGTTCAGAGACCGGAGGCGTGTGTACGTGCGTTTCGAATTCCTGTAAGCCCACC";

const HTO_SEQUENCE: &str =
    "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATA";

const GENOMIC_READ: &str =
    "ATCGATGCTAGCTACGATCGTACGCTAGCATGCTACGATCGTAGCTACGATGCTAGCATC";

fn test_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(TEST_DATA)
        .join(name)
}

fn star_available() -> bool {
    Command::new("STAR")
        .arg("--version")
        .output()
        .map(|out| out.status.success())
        .unwrap_or(false)
}

fn ensure_star_index(fasta: &Path, gtf: &Path, star_index: &Path) {
    if star_index.join("Genome").is_file() {
        return;
    }

    assert!(
        star_available(),
        "STAR is required for this integration test but was not found in PATH"
    );

    fs::create_dir_all(star_index)
        .expect("failed to create STAR index directory");

    let status = Command::new("STAR")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            star_index.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            "--genomeSAindexNbases",
            "2",
            "--sjdbOverhang",
            "59",
        ])
        .status()
        .expect("failed to start STAR genomeGenerate");

    assert!(status.success(), "STAR genomeGenerate failed");
}

fn ensure_splice_index(gtf: &Path, index_path: &Path) {
    if index_path.is_file() {
        return;
    }

    let index = AnnotationBuilder::new(1_000)
        .build_from_path(gtf)
        .expect("failed to build tiny splice index");

    index
        .save(index_path)
        .expect("failed to write tiny splice index");

    SpliceIndex::load(index_path)
        .expect("created splice index cannot be loaded");
}

fn write_fastq(path: &Path, reads: &[(&str, &str)]) {
    let mut text = String::new();

    for (name, seq) in reads {
        text.push('@');
        text.push_str(name);
        text.push('\n');
        text.push_str(seq);
        text.push_str("\n+\n");
        text.push_str(&"I".repeat(seq.len()));
        text.push('\n');
    }

    fs::write(path, text)
        .unwrap_or_else(|err| panic!("failed to write {}: {err}", path.display()));
}

fn read_gzip_text(path: &Path) -> String {
    let file = fs::File::open(path)
        .unwrap_or_else(|err| panic!("failed to open {}: {err}", path.display()));

    let mut decoder = MultiGzDecoder::new(file);
    let mut text = String::new();
    decoder
        .read_to_string(&mut text)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));

    text
}

fn assert_feature_output(out: &Path, feature_type: &str, feature_name: &str) {
    let dir = out.join(feature_type);
    let features = dir.join("features.tsv.gz");

    assert!(
        dir.is_dir(),
        "Nelrune did not create expected additional-feature output folder {}",
        dir.display()
    );

    for file in ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"] {
        assert!(
            dir.join(file).is_file(),
            "missing expected additional-feature matrix output: {}",
            dir.join(file).display()
        );
    }

    let text = read_gzip_text(&features);
    assert!(
        text.lines().any(|line| {
            let mut fields = line.split('\t');
            fields.next() == Some(feature_name)
                && fields.next() == Some(feature_name)
                && fields.next() == Some(feature_type)
        }),
        "expected feature {feature_name:?} with type {feature_type:?} in {}:\n{text}",
        features.display()
    );
}

#[test]
fn integration_additional_features_preserves_bd_rhapsody_and_hto_names() {
    let r1 = test_path("additional_features_R1.fastq");
    let r2 = test_path("additional_features_R2.fastq");
    let hto = test_path("hto.fasta");

    let fasta = test_path("tiny.fa");
    let gtf = test_path("tiny.gtf");
    let splice_index = test_path("tiny.gtf.dat");
    let star_index = test_path("star_index");

    for path in [&fasta, &gtf] {
        assert!(
            path.is_file(),
            "required test input is missing: {}",
            path.display()
        );
    }

    ensure_star_index(&fasta, &gtf, &star_index);
    ensure_splice_index(&gtf, &splice_index);

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("test-integration-additional-features");

    // Remove artifacts from the PREVIOUS test run.
    // Deliberately do not clean up at the end so failed runs can be inspected.
    if out.exists() {
        fs::remove_dir_all(&out)
            .expect("failed to remove previous integration-test directory");
    }

    fs::create_dir_all(&out)
        .expect("failed to create integration-test directory");

    let output = Command::new(env!("CARGO_BIN_EXE_nelrune"))
        .args([
            "--r1",
            r1.to_str().unwrap(),
            "--r2",
            r2.to_str().unwrap(),
            "--primer-structure",
            "CELL:4+UMI:4",
            "--additional-features",
            "bd_sample_mouse",
            hto.to_str().unwrap(),
            "--additional-feature-min-hits",
            "4",
            "--mapper",
            "star",
            "--mapper-index",
            star_index.to_str().unwrap(),
            "--mapper-threads",
            "2",
            "--index",
            splice_index.to_str().unwrap(),
            "--require-strand",
            "--min-mapq",
            "0",
            "--min-cell-counts",
            "1",
            "--min-insert-len",
            "20",
            "--threads",
            "2",
            "--outpath",
            out.to_str().unwrap(),
            "--no-health-server",
        ])
        .output()
        .expect("failed to start nelrune");

    if !output.status.success() {
        panic!(
            "\nNelrune additional-feature integration test failed.\n\
             status: {}\n\
             stdout:\n{}\n\
             stderr:\n{}\n",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }

    // Built-in BD Rhapsody sample tags must retain both their real sample name
    // and the built-in data type used for the output folder.
    assert_feature_output(
        &out,
        "bd_sample_mouse",
        "SampleTag07_mm",
    );

    // For user FASTA input, the filename defines the data type/output folder
    // while the FASTA header defines the individual feature name.
    assert_feature_output(
        &out,
        "hto",
        "Donor_A_HTO",
    );
}
