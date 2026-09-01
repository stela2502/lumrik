use std::fs;
use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, Result, bail};
use tempfile::TempDir;

use bam_tide::fastq::FastqRecord;

use sc_mapper::core::MapperLaunch;
use sc_mapper::traits::ExternalMapper;

use sc_mapper::Bwa;

fn bwa_available() -> bool {
    match Command::new("bwa").arg("2>&1").output() {
        Ok(_) => true,

        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("skipping BWA integration test: bwa not found");
            false
        }

        Err(err) => {
            eprintln!("skipping BWA integration test: {err}");
            false
        }
    }
}

fn prepare_bwa_reference() -> Result<(TempDir, PathBuf)> {
    let tmpdir = tempfile::tempdir().context("failed to create temporary BWA test directory")?;

    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/tiny.fa");

    let reference = tmpdir.path().join("tiny.fa");

    fs::copy(&source, &reference).with_context(|| {
        format!(
            "failed to copy test reference {} -> {}",
            source.display(),
            reference.display()
        )
    })?;

    let status = Command::new("bwa")
        .arg("index")
        .arg(&reference)
        .status()
        .context("failed to run `bwa index`")?;

    if !status.success() {
        bail!("`bwa index` failed with status {status}");
    }

    Ok((tmpdir, reference))
}

#[test]
fn bwa_maps_a_real_fastq_record() -> Result<()> {
    if !bwa_available() {
        eprintln!("bwa is not available - abort");
        return Ok(());
    }

    let (_tmpdir, reference) = prepare_bwa_reference()?;

    let launch = MapperLaunch {
        mapper_bin: "bwa".into(),
        index: reference,
        options: vec!["mem".into()],
        threads: 2,
        paired: false,
    };

    let mapper = Bwa::from_launch(launch);

    let mut mapper = mapper.spawn()?;

    let read = fq("read_001", "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA");

    let mut calls = Vec::new();

    if let Some(call) = mapper.process(&read, None)? {
        calls.push(call);
    }

    calls.extend(mapper.finish()?);

    assert_eq!(
        calls.len(),
        1,
        "expected exactly one MappingCall for one submitted BWA read"
    );

    let call = &calls[0];

    assert_eq!(call.read_id, "read_001");

    assert!(
        !call.records.records.is_empty(),
        "BWA returned no SAM records"
    );

    assert!(
        call.records
            .records
            .iter()
            .any(|rec| !rec.record.is_unmapped()),
        "read_001 should map to the tiny test reference"
    );

    Ok(())
}

fn fq(id: &str, seq: &str) -> FastqRecord {
    FastqRecord {
        id: id.to_string(),
        seq: seq.as_bytes().to_vec(),
        qual: vec![b'I'; seq.len()],
    }
}
