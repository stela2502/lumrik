//! Regression for the five real ONT reads containing 15 independently barcoded
//! 10x cassettes. No mapper or quantifier is involved.
//! Place the five-read fixture in tests/data/ont_cassette_candidates.bam.
use rust_htslib::bam::{Read, Reader};
use sc_primer::{Grammar, PrimerDetector};

#[test]
fn five_ont_reads_have_fifteen_10x_cassettes() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/ont_cassette_candidates.bam");
    let grammar = Grammar::parse(
        "ont-10x-v4-multiplex-regression",
        "TYPE:GEX+FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:3p-v4+UMI:12",
    ).unwrap();
    let detector = PrimerDetector::from_grammar(grammar).unwrap()
        .with_read_wide_primer_search(true);
    let mut bam = Reader::from_path(path).unwrap();
    let mut reads = 0;
    let mut cassettes = 0;
    let mut multi = 0;
    for record in bam.records() {
        let record = record.unwrap();
        let seq = record.seq().as_bytes();
        let hits = detector.detect_all(&seq, &record.qual()).unwrap();
        for hit in &hits {
            assert!(hit.get_cell(&seq, &record.qual()).is_ok());
            assert!(hit.get_umi(&seq, &record.qual()).is_ok());
        }
        reads += 1;
        cassettes += hits.len();
        multi += usize::from(hits.len() > 1);
    }
    assert_eq!(reads, 5, "regression fixture must contain exactly five reads");
    assert_eq!(multi, 5, "all five reads must contain multiple cassettes");
    assert_eq!(cassettes, 15, "regression: ONT must find 15, not 1");
}
