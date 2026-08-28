use anyhow::Result;

use bam_tide::fastq::FastqRecord;
use bam_tide::illumina_normalizer::{
    IlluminaNormalizerConfig,
    IlluminaPartial,
};
use bam_tide::illumina_normalizer::cli::{
    InsertRead,
    PrimerRead,
};
use sc_primer::{Grammar, PrimerDetector};

use std::path::PathBuf;


fn test_primer_detector() -> PrimerDetector {
    let grammar = Grammar::parse(
        "illumina-test",
        "CELL:8+UMI:6",
    )
    .expect("failed to create test primer grammar");

    PrimerDetector::from_grammar(grammar)
        .expect("failed to create test primer detector")
}


fn test_config() -> IlluminaNormalizerConfig {
    IlluminaNormalizerConfig {
        out: PathBuf::from("unused-out.fastq"),
        read_tags: PathBuf::from("unused-read-tags.tsv"),

        primer_read: PrimerRead::R1,
        insert_read: InsertRead::R2,

        primer: test_primer_detector(),
        additional_features: Vec::new(),
        additional_feature_min_hits: 4,

        max_reads: Some(10),

        min_insert_len: 30,
        threads: 1,
        gzip_level: 1,
        gzip: false,
    }
}


fn fastq(id: &str, seq: &[u8]) -> FastqRecord {
    FastqRecord::new(
        id,
        seq,
        &vec![b'I'; seq.len()],
    )
}


#[test]
fn illumina_normalize_pair_splits_cell_umi_and_insert() -> Result<()> {
    let config = test_config();

    /*
        CELL      UMI      INSERT
        AAAAAAAA  CCCCCC   GATTACAGATTACAGATTACAGATTACAGATTACA
    */
    let r1 = fastq(
        "read1",
        b"AAAAAAAACCCCCCGATTACAGATTACAGATTACAGATTACAGATTACA",
    );

    let r2_seq =
        b"TGCATGCATGCATGCATGCATGCATGCATGCA";

    let r2 = fastq(
        "read1",
        r2_seq,
    );

    let mut partial = IlluminaPartial::new();

    let feature_mapper = fast_tag_mapper::FastTagMapper::new();
    partial.normalize_pair(
        &r1,
        &r2,
        &config,
        &feature_mapper,
    )?;

    assert_eq!(
        partial.candidates.len(),
        1,
        "expected exactly one Illumina candidate"
    );

    let candidate = &partial.candidates[0];

    /*
     * R2 must survive unchanged apart from its normalized ID.
     */
    assert_eq!(
        candidate.fastq_record.seq,
        r2_seq,
        "R2 sequence was changed"
    );

    //eprintln!("reads_log   = {:?}", partial.stats.reads_log);
    //eprintln!("error_counts = {:?}", partial.stats.error_counts);
    /*
     * The remainder of R1 after CELL+UMI must become
     * the paired R1 insert.
     */
    let paired_r1 = candidate
        .paired_r1_record
        .as_ref()
        .expect("expected usable paired R1 insert");

    assert_eq!(
        paired_r1.seq,
        b"GATTACAGATTACAGATTACAGATTACAGATTACA",
        "wrong R1 insert extracted"
    );

    /*
     * The read tag must contain the CELL and UMI extracted
     * from R1.
     */
    assert_eq!(
        candidate.read_tag.cell_seq,
        b"AAAAAAAA"
    );

    assert_eq!(
        candidate.read_tag.umi_seq,
        b"CCCCCC"
    );


    assert_eq!(
        partial
            .stats
            .reads_log
            .get("total_pairs")
            .copied()
            .unwrap_or(0),
        1
    );

    assert_eq!(
        partial
            .stats
            .reads_log
            .get("candidate_pairs")
            .copied()
            .unwrap_or(0),
        1
    );

    Ok(())
}


#[test]
fn illumina_normalize_pair_keeps_r2_when_r1_insert_is_too_short() -> Result<()> {
    let config = test_config();

    /*
        CELL      UMI      INSERT
        AAAAAAAA  CCCCCC   GATTACA

        INSERT is far shorter than the current 30 bp
        usable_insert() requirement.
    */
    let r1 = fastq(
        "read1",
        b"AAAAAAAACCCCCCGATTACA",
    );

    let r2_seq =
        b"TGCATGCATGCATGCATGCATGCATGCATGCA";

    let r2 = fastq(
        "read1",
        r2_seq,
    );

    let mut partial = IlluminaPartial::new();

    let feature_mapper = fast_tag_mapper::FastTagMapper::new();
    partial.normalize_pair(
        &r1,
        &r2,
        &config,
        &feature_mapper,
    )?;

    assert_eq!(
        partial.candidates.len(),
        1,
        "R2 must still become a candidate"
    );

    let candidate = &partial.candidates[0];

    assert!(
        candidate.paired_r1_record.is_none(),
        "short R1 insert must not become a paired read"
    );

    assert_eq!(
        candidate.fastq_record.seq,
        r2_seq,
        "R2 must survive even when the R1 insert is unusable"
    );

    assert_eq!(
        candidate.read_tag.cell_seq,
        b"AAAAAAAA"
    );

    assert_eq!(
        candidate.read_tag.umi_seq,
        b"CCCCCC"
    );

    assert_eq!(
        partial
            .stats
            .reads_log
            .get("candidate_pairs")
            .copied()
            .unwrap_or(0),
        1
    );

    assert_eq!(
        partial
            .stats
            .reads_log
            .get("no_usable_paired_r1_insert")
            .copied()
            .unwrap_or(0),
        1
    );

    Ok(())
}


#[test]
fn illumina_normalize_pair_rejects_read_without_complete_primer() {
    let config = test_config();

    /*
     * Only four bases: impossible to extract CELL:8 + UMI:6.
     */
    let r1 = fastq(
        "read1",
        b"AAAA",
    );

    let r2 = fastq(
        "read1",
        b"TGCATGCATGCATGCATGCATGCATGCATGCA",
    );

    let mut partial = IlluminaPartial::new();

    let feature_mapper = fast_tag_mapper::FastTagMapper::new();
    let result = partial.normalize_pair(
        &r1,
        &r2,
        &config,
        &feature_mapper,
    );

    assert!(
        result.is_err(),
        "incomplete primer must be rejected"
    );

    assert!(
        partial.candidates.is_empty(),
        "rejected read must not create a candidate"
    );

    assert_eq!(
        partial
            .stats
            .reads_log
            .get("no_primer_match")
            .copied()
            .unwrap_or(0),
        1
    );
}


#[test]
fn illumina_normalize_pair_preserves_cell_and_umi_qualities() -> Result<()> {
    let config = test_config();

    let seq =
        b"AAAAAAAACCCCCCGATTACAGATTACAGATTACAGATTACAGATTACA";

    let mut qual = vec![b'I'; seq.len()];

    /*
     * Give CELL and UMI unmistakably different qualities so
     * we can prove the slices come from the correct positions.
     */
    qual[0..8].fill(b'A');
    qual[8..14].fill(b'B');

    let r1 = FastqRecord::new(
        "read1",
        seq,
        &qual,
    );

    let r2 = fastq(
        "read1",
        b"TGCATGCATGCATGCATGCATGCATGCATGCA",
    );

    let mut partial = IlluminaPartial::new();

    let feature_mapper = fast_tag_mapper::FastTagMapper::new();
    partial.normalize_pair(
        &r1,
        &r2,
        &config,
        &feature_mapper,
    )?;

    assert_eq!(partial.candidates.len(), 1);

    let candidate = &partial.candidates[0];

    assert_eq!(
        candidate.read_tag.cell_qual,
        vec![b'A'; 8]
    );

    assert_eq!(
        candidate.read_tag.umi_qual,
        vec![b'B'; 6]
    );

    Ok(())
}