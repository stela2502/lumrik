use std::fs;
use std::path::PathBuf;
use std::process::Command;

use anyhow::{Context, Result, bail};
use bam_tide::fastq::FastqRecord;
use tempfile::TempDir;

use sc_mapper::core::MapperLaunch;
use sc_mapper::process::Star;
use sc_mapper::traits::ExternalMapper;
use sc_mapper::{MapperKind, StreamingMapperCli};

fn fq(id: &str, seq: &str) -> FastqRecord {
    FastqRecord {
        id: id.to_string(),
        seq: seq.as_bytes().to_vec(),
        qual: vec![b'I'; seq.len()],
    }
}

fn star_available() -> bool {
    match Command::new("STAR").arg("--version").output() {
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

#[test]
fn star_from_cli_requires_input_completion_before_header() -> Result<()> {
    if !star_available() {
        return Ok(());
    }

    let (_tmpdir, genome_dir) = prepare_star_reference()?;

    let cli = StreamingMapperCli {
        mapper: MapperKind::Star,
        mapper_bin: None,
        mapper_index: genome_dir,
        mapper_options: None,
        mapper_threads: 1,
        mapper_paired: false,
        mapper_keep_multimappers: true,
    };

    let mut mapper = cli.from_cli().context("starting STAR through from_cli")?;

    let read = fq("read_001", "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA");

    mapper.submit(&read, None)?;

    std::thread::sleep(std::time::Duration::from_secs(1));

    eprintln!("header before finish: {}", mapper.header_loaded());

    /*let calls =
        mapper.finish()
            .context("finishing STAR")?;

    eprintln!(
        "calls after finish: {}",
        calls.len()
    );
    */
    assert!(
        mapper.header_loaded(),
        "STAR still produced no header after input was closed"
    );

    let header = mapper.header()?.clone();

    let header_text = String::from_utf8(header.to_bytes())?;

    assert!(
        header_text.contains("@SQ"),
        "STAR header contains no @SQ records:\n{header_text}"
    );

    Ok(())
}

fn prepare_star_reference() -> Result<(TempDir, PathBuf)> {
    let tmpdir = tempfile::tempdir().context("failed to create temporary STAR test directory")?;

    let source = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/tiny.fa");

    let fasta = tmpdir.path().join("tiny.fa");

    fs::copy(&source, &fasta).with_context(|| {
        format!(
            "failed to copy STAR test reference {} -> {}",
            source.display(),
            fasta.display()
        )
    })?;

    let genome_dir = tmpdir.path().join("star-index");

    fs::create_dir(&genome_dir).context("failed to create temporary STAR genome directory")?;

    let output = Command::new("STAR")
        .current_dir(tmpdir.path())
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
        .output()
        .context("failed to run STAR genomeGenerate")?;

    if !output.status.success() {
        bail!(
            "STAR genomeGenerate failed with status {}\n\nstdout:\n{}\n\nstderr:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }

    Ok((tmpdir, genome_dir))
}

#[test]
fn star_maps_a_real_fastq_record() -> Result<()> {
    if !star_available() {
        return Ok(());
    }

    let (_tmpdir, genome_dir) = prepare_star_reference()?;

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

    let mapper = Star::from_launch(launch)?;

    std::thread::sleep(std::time::Duration::from_millis(500));

    mapper.check()?;

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
        "expected exactly one MappingCall for one STAR input read"
    );

    let call = &calls[0];

    assert_eq!(call.read_id, "read_001");

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

#[test]
fn star_from_cli_maps_after_input_is_closed() -> Result<()> {
    if !star_available() {
        return Ok(());
    }

    let (_tmpdir, genome_dir) = prepare_star_reference()?;

    let cli = StreamingMapperCli {
        mapper: MapperKind::Star,
        mapper_bin: None,
        mapper_index: genome_dir,
        mapper_options: None,
        mapper_threads: 1,
        mapper_paired: false,
        mapper_keep_multimappers: true,
    };

    let mut mapper = cli
        .from_cli()
        .context("starting STAR through StreamingMapperCli::from_cli")?;

    assert!(
        mapper.is_running()?,
        "STAR exited immediately after startup"
    );

    let read = fq("read_001", "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA");

    mapper.submit(&read, None)?;

    let calls = mapper.finish().context("finishing STAR")?;

    assert_eq!(calls.len(), 1, "expected exactly one STAR MappingCall");

    assert_eq!(calls[0].read_id, "read_001");

    assert!(
        !calls[0].records.records.is_empty(),
        "STAR returned no SAM records"
    );

    Ok(())
}
