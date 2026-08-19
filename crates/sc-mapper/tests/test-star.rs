use std::fs;
use std::path::PathBuf;
use std::process::Command;

use anyhow::{bail, Context, Result};
use bam_tide::fastq::FastqRecord;
use tempfile::TempDir;

use sc_mapper::core::MapperLaunch;
use sc_mapper::process::Star;
use sc_mapper::traits::ExternalMapper;


fn fq(id: &str, seq: &str) -> FastqRecord {
    FastqRecord {
        id: id.to_string(),
        seq: seq.as_bytes().to_vec(),
        qual: vec![b'I'; seq.len()],
    }
}


fn star_available() -> bool {
    match Command::new("STAR")
        .arg("--version")
        .output()
    {
        Ok(_) => true,

        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("skipping STAR integration test: STAR not found");
            false
        }

        Err(err) => {
            eprintln!("skipping STAR integration test: {err}");
            false
        }
    }
}


fn prepare_star_reference() -> Result<(TempDir, PathBuf)> {
    let tmpdir = tempfile::tempdir()
        .context("failed to create temporary STAR test directory")?;

    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/tiny.fa");

    let fasta = tmpdir.path().join("tiny.fa");

    fs::copy(&source, &fasta).with_context(|| {
        format!(
            "failed to copy STAR test reference {} -> {}",
            source.display(),
            fasta.display()
        )
    })?;

    let genome_dir = tmpdir.path().join("star-index");

    fs::create_dir(&genome_dir)
        .context("failed to create temporary STAR genome directory")?;

    let status = Command::new("STAR")
        .arg("--runMode")
        .arg("genomeGenerate")
        .arg("--runThreadN")
        .arg("1")
        .arg("--genomeDir")
        .arg(&genome_dir)
        .arg("--genomeFastaFiles")
        .arg(&fasta)

        // Tiny genomes need a tiny suffix-array index parameter.
        .arg("--genomeSAindexNbases")
        .arg("1")

        .status()
        .context("failed to run STAR genomeGenerate")?;

    if !status.success() {
        bail!(
            "STAR genomeGenerate failed with status {status}"
        );
    }

    Ok((tmpdir, genome_dir))
}


#[test]
fn star_maps_a_real_fastq_record() -> Result<()> {
    if !star_available() {
        return Ok(());
    }

    let (_tmpdir, genome_dir) =
        prepare_star_reference()?;

    let launch = MapperLaunch {
        mapper_bin: "STAR".into(),
        index: genome_dir,
        options: vec![
            "--runThreadN".into(),
            "1".into(),

            "--outSAMtype".into(),
            "SAM".into(),

            "--outStd".into(),
            "SAM".into(),
        ],
        threads: 2,
        paired: false,
    };

    let mapper = Star::from_launch(launch);

    mapper.check()?;

    let mut mapper = mapper.spawn()?;

    let read = fq(
        "read_001",
        "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA",
    );

    let mut calls = Vec::new();

    if let Some(call) =
        mapper.process(&read, None)?
    {
        calls.push(call);
    }

    calls.extend(mapper.finish()?);

    assert_eq!(
        calls.len(),
        1,
        "expected exactly one MappingCall for one STAR input read"
    );

    let call = &calls[0];

    assert_eq!(
        call.read_id,
        "read_001"
    );

    assert!(
        !call.records.records.is_empty(),
        "STAR returned no SAM records"
    );

    assert!(
        call.records
            .records
            .iter()
            .any(|rec| !rec.record.is_unmapped()),
        "read_001 should map to the tiny STAR reference"
    );

    Ok(())
}