use std::path::Path;

use anyhow::Result;
use bam_tide::fastq::FastqRecord;

use sc_mapper::core::MapperLaunch;
use sc_mapper::process::{MapperProcess, Minimap2};
use sc_mapper::traits::ExternalMapper;
use sc_mapper::core::MapperProcessLike;

fn fq(id: &str, seq: &str) -> FastqRecord {
    FastqRecord {
        id: id.to_string(),
        seq: seq.as_bytes().to_vec(),
        qual: vec![b'I'; seq.len()],
    }
}

#[test]
fn single_stdin_process_writes_fastq_and_drains_all_sam_records() -> Result<()> {
    // Fake mapper:
    // 1. consumes FASTQ from stdin until EOF
    // 2. emits a minimal SAM header
    // 3. emits two alignments for the same read
    //
    // Because output is only written after stdin reaches EOF, this also checks
    // that finish() closes the mapper input before draining stdout.
    let args = vec![
        "-c".to_string(),
        r#"
            cat >/dev/null
            printf '@HD\tVN:1.6\tSO:unsorted\n'
            printf '@SQ\tSN:chr1\tLN:1000\n'
            printf '@SQ\tSN:chr2\tLN:1000\n'
            printf 'read1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
            printf 'read1\t256\tchr2\t1\t20\t4M\t*\t0\t0\tACGT\tIIII\n'
        "#
        .to_string(),
    ];

    let mut process =
        MapperProcess::spawn_single_stdin(Path::new("sh"), &args, None)?;

    let r1 = fq("read1", "ACGT");

    process.write_fastq(&r1, None)?;

    // The fake mapper deliberately produces nothing before EOF.
    assert!(
        process.next_cluster()?.is_none(),
        "fake mapper should not emit a cluster before input is closed"
    );

    let clusters = Box::new(process).finish()?;

    assert_eq!(clusters.len(), 1);

    let cluster = &clusters[0];

    assert_eq!(cluster.read_id, "read1");
    assert_eq!(cluster.records.len(), 2);
    assert_eq!(cluster.records[0].qname(), b"read1");
    assert_eq!(cluster.records[1].qname(), b"read1");

    Ok(())
}


fn minimap2() -> Minimap2 {
    Minimap2::from_launch(MapperLaunch {
        mapper_bin: "minimap2".into(),
        index: "tests/data/tiny.fa".into(),
        options: vec!["-ax".into(), "sr".into()],
        threads: 2,
        paired: false,
    })
}

fn minimap2_available(mapper: &Minimap2) -> bool {
    match mapper.check() {
        Ok(()) => true,
        Err(err) => {
            eprintln!("skipping minimap2 integration test: {err}");
            false
        }
    }
}

#[test]
fn minimap2_check_is_useful_when_binary_is_present() -> Result<()> {
    let mapper = minimap2();

    if !minimap2_available(&mapper) {
        return Ok(());
    }

    mapper.check()?;

    Ok(())
}

#[test]
fn minimap2_maps_a_real_fastq_record() -> Result<()> {
    let mapper = minimap2();

    if !minimap2_available(&mapper) {
        return Ok(());
    }

    let mut mapper = mapper.spawn()?;

    // This sequence must occur uniquely in tests/data/tiny.fa.
    //
    // Make tiny.fa deterministic for this test, for example:
    //
    // >chrTest
    // TTTTTTTTTTACGTTGCAACGTTGCAACGTTGCAACGTTGCAAAAAAAAAAA
    //
    // The read below then has one obvious alignment.
    let read = fq(
        "read_001",
        "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA",
    );

    let mut calls = Vec::new();

    // process() submits the FASTQ record and may opportunistically return one
    // completed result. For a one-read test this will often be None, which is
    // perfectly valid.
    if let Some(call) = mapper.process(&read, None)? {
        calls.push(call);
    }

    // Closing mapper input is essential: minimap2 may buffer output until EOF.
    // finish() must return every MappingCall that was still outstanding.
    calls.extend(mapper.finish()?);

    assert_eq!(
        calls.len(),
        1,
        "expected exactly one completed MappingCall for one submitted read"
    );

    let call = &calls[0];

    assert_eq!(call.read_id, "read_001");
    assert_eq!(call.records.records.len(), 1);

    let record = &call.records.records[0].record;

    assert!(
        !record.is_unmapped(),
        "read_001 should map to tests/data/tiny.fa"
    );

    Ok(())
}

#[test]
fn minimap2_streaming_returns_every_submitted_read_exactly_once() -> Result<()> {
    let mapper = minimap2();

    if !minimap2_available(&mapper) {
        return Ok(());
    }

    let mut mapper = mapper.spawn()?;

    let reads = [
        fq(
            "read_001",
            "ACGTTGCAACGTTGCAACGTTGCAACGTTGCA",
        ),
        fq(
            "read_002",
            "CGTTGCAACGTTGCAACGTTGCAACGTTGCAA",
        ),
        fq(
            "read_003",
            "GTTGCAACGTTGCAACGTTGCAACGTTGCAAC",
        ),
    ];

    let mut calls = Vec::new();

    for read in &reads {
        if let Some(call) = mapper.process(read, None)? {
            calls.push(call);
        }

        // Drain anything else that is already ready without blocking.
        while let Some(call) = mapper.try_next()? {
            calls.push(call);
        }
    }

    calls.extend(mapper.finish()?);

    assert_eq!(
        calls.len(),
        reads.len(),
        "every submitted FASTQ record must yield exactly one MappingCall"
    );

    let mut received_ids: Vec<_> =
        calls.iter().map(|call| call.read_id.as_str()).collect();

    received_ids.sort_unstable();

    let mut expected_ids: Vec<_> =
        reads.iter().map(|read| read.clean_id()).collect();

    expected_ids.sort_unstable();

    assert_eq!(received_ids, expected_ids);

    for call in &calls {
        assert!(
            !call.records.records.is_empty(),
            "{} returned no mapper records",
            call.read_id
        );
    }

    Ok(())
}