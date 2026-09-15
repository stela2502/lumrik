use pretty_assertions::assert_eq;
use read_tag_table::ReadTagRecord;

use sc_primer::{
    BdCellVersion, Chemistry, Grammar, Orientation, PrimerDetector, RhapsodyWhitelist,
};

struct TestData;

impl TestData {
    fn qual(len: usize) -> Vec<u8> {
        vec![b'I'; len]
    }

    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        let mut out = Vec::new();
        for part in parts {
            out.extend_from_slice(part);
        }
        out
    }

    fn rc(seq: &[u8]) -> Vec<u8> {
        PrimerDetector::reverse_complement(seq)
    }

    fn tenx_3p_v3_grammar(name: &str) -> Grammar {
        Grammar::parse(
            name,
            "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10",
        )
        .unwrap()
    }

    fn tenx_3p_v3_no_polyt_grammar(name: &str) -> Grammar {
        Grammar::parse(name, "FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12").unwrap()
    }

    fn tenx_cell() -> &'static [u8] {
        b"AACCGGTTAACCGGTT"
    }

    fn tenx_umi() -> &'static [u8] {
        b"TTAACCGGTTAA"
    }

    fn tenx_insert() -> &'static [u8] {
        b"GACCTGACTGACTGACCTGA"
    }

    fn synthesized_tenx_read(grammar: &Grammar) -> Vec<u8> {
        let mut read = grammar
            .synthesize(Self::tenx_cell(), Self::tenx_umi())
            .unwrap();
        read.extend_from_slice(Self::tenx_insert());
        read
    }

    fn bd_v2_384_grammar(name: &str) -> Grammar {
        Grammar::parse(
            name,
            "FIXED:ATAGGAAACTCATGGT:mm=2+BD_CELL:v2.384+POLYT:min=0",
        )
        .unwrap()
    }

    fn bd_v2_384_cell(index: u64) -> Vec<u8> {
        let wl = RhapsodyWhitelist::bd_v2_384();
        let (c1, c2, c3) = wl
            .cell_id_to_parts_ids(index + 1)
            .expect("test requested invalid BD cell id");
        wl.create_cell_cassette(c1, c2, c3)
    }
}

struct TenxOntMultimerTest;

struct BuiltRead {
    seq: Vec<u8>,
    qual: Vec<u8>,
    cells: Vec<Vec<u8>>,
    umis: Vec<Vec<u8>>,
    inserts: Vec<Vec<u8>>,
}

impl TenxOntMultimerTest {
    fn dna_tag(index: usize, len: usize) -> Vec<u8> {
        let alphabet = b"ACGT";
        let mut out = vec![b'A'; len];
        let mut x = index;

        for pos in (0..len).rev() {
            out[pos] = alphabet[x & 3];
            x >>= 2;
        }

        out
    }

    fn cell(index: usize) -> Vec<u8> {
        let mut v = b"ACGTACGTACGT".to_vec();
        v.extend(Self::dna_tag(index, 4));
        v
    }

    fn umi(index: usize) -> Vec<u8> {
        let mut v = b"TTGGAACC".to_vec();
        v.extend(Self::dna_tag(index, 4));
        v
    }

    fn insert(index: usize) -> Vec<u8> {
        let mut v = b"GATCGATCGATCGATCGATCGATC".to_vec();
        v.extend(Self::dna_tag(index, 4));
        v
    }

    fn build_read(grammar: &Grammar) -> BuiltRead {
        let mut out = BuiltRead {
            seq: Vec::new(),
            qual: Vec::new(),
            cells: Vec::new(),
            umis: Vec::new(),
            inserts: Vec::new(),
        };

        for index in 0..10 {
            let cell = Self::cell(index);
            let umi = Self::umi(index);
            let insert = Self::insert(index);
            let monomer_quality = b'!'.saturating_add(index as u8);

            let mut primer = grammar.synthesize(&cell, &umi).unwrap();
            primer.extend_from_slice(&insert);

            out.qual
                .extend(std::iter::repeat_n(monomer_quality, primer.len()));
            out.seq.extend_from_slice(&primer);

            out.cells.push(cell);
            out.umis.push(umi);
            out.inserts.push(insert);
        }

        out
    }

    fn build_chemistry_monomers(detector: &PrimerDetector, fuzzy: bool) -> Vec<(Vec<u8>, Vec<u8>)> {
        let mut reads = Vec::with_capacity(10);

        for index in 0..10 {
            let mut cell = detector
                .cell_seq_for_index(index as u64)
                .expect("missing chemistry whitelist cell");

            if fuzzy {
                cell[0] = b'N';
            }

            let umi = detector.grammar().umi_from_u64(index as u64);
            let insert = Self::insert(index);
            let monomer_quality = b'!'.saturating_add(index as u8);

            let mut seq = detector.grammar().synthesize(&cell, &umi).unwrap();
            seq.extend_from_slice(&insert);
            let qual = vec![monomer_quality; seq.len()];

            reads.push((seq, qual));
        }

        reads
    }

    fn run() {
        let grammar = TestData::tenx_3p_v3_no_polyt_grammar("tenx-ont-stress");
        let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();
        let built = Self::build_read(&grammar);

        let matches = detector.detect_all(&built.seq, &built.qual).unwrap();

        assert_eq!(matches.len(), 10);

        for (index, primer_match) in matches.iter().enumerate() {
            let cell = primer_match.get_cell(&built.seq, &built.qual).unwrap();
            let umi = primer_match.get_umi(&built.seq, &built.qual).unwrap();
            let insert = primer_match.get_insert(&built.seq, &built.qual).unwrap();

            let expected_quality = vec![b'!'.saturating_add(index as u8); cell.seq.len()];

            assert_eq!(
                primer_match.orientation,
                Orientation::Forward,
                "{index}: orientation mismatch"
            );

            assert_eq!(
                cell.seq,
                built.cells[index],
                "{index}: cell seq got {:?} expected {:?}",
                std::str::from_utf8(&cell.seq),
                std::str::from_utf8(&built.cells[index]),
            );

            assert_eq!(
                umi.seq,
                built.umis[index],
                "{index}: umi seq got {:?} expected {:?}",
                std::str::from_utf8(&umi.seq),
                std::str::from_utf8(&built.umis[index]),
            );

            assert_eq!(
                insert.seq,
                built.inserts[index],
                "{index}: insert seq got {:?} expected {:?}",
                std::str::from_utf8(&insert.seq),
                std::str::from_utf8(&built.inserts[index]),
            );

            assert_eq!(cell.qual, expected_quality, "{index}: cell qual mismatch");
        }
    }
}

pub fn bench<F: FnMut()>(name: &str, iterations: usize, hits_per_op: usize, mut f: F) {
    use std::time::Instant;

    let start = Instant::now();

    for _ in 0..iterations {
        std::hint::black_box(f());
    }

    let elapsed = start.elapsed();

    let us_per_op = elapsed.as_secs_f64() * 1_000_000.0 / iterations as f64;

    let us_per_hit = us_per_op / hits_per_op as f64;

    eprintln!(
        "{:<40} {:>10.3} µs/op {:>10.3} µs/hit",
        name, us_per_op, us_per_hit,
    );
}

#[test]
fn test_parse_complicated_custom_grammar() {
    let grammar = Grammar::parse(
        "custom",
        "SEARCH:0..4+FIXED:CTACACGACGCTCTTCCGATCT:mm=2+CELL:16+UMI:12+POLYT:min=10",
    )
    .unwrap();

    assert_eq!(grammar.ops.len(), 5);
    assert_eq!(grammar.cell_len(), 16);
    assert_eq!(grammar.umi_len(), 12);
}

#[test]
fn test_tenx_roundtrip_synthesize_detect_polyt_insert() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-roundtrip");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let seq = TestData::synthesized_tenx_read(&grammar);
    let qual = TestData::qual(seq.len());

    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.orientation, Orientation::Forward);
    assert_eq!(
        hit.get_cell(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_cell()
    );
    assert_eq!(
        hit.get_umi(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_umi()
    );
    assert_eq!(
        hit.get_insert(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_insert()
    );
}

#[test]
fn test_custom_10x_adapter_allows_one_fixed_error_after_synthesis() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-fixed-error");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let mut seq = TestData::cat(&[b"NN", &TestData::synthesized_tenx_read(&grammar)]);
    seq[2 + 5] = b'A';

    let qual = TestData::qual(seq.len());
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.primer_start, 2);
    assert_eq!(
        hit.get_cell(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_cell()
    );
    assert_eq!(
        hit.get_umi(&seq, &qual).unwrap().seq.as_slice(),
        TestData::tenx_umi()
    );
}

#[test]
fn test_custom_10x_adapter_rejects_too_many_fixed_errors() {
    let grammar = Grammar::parse(
        "tenx-fixed-error-reject",
        "FIXED:CTACACGACGCTCTTCCGATCT:mm=1+CELL:16+UMI:12+POLYT:min=10",
    )
    .unwrap();
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let mut seq = TestData::synthesized_tenx_read(&grammar);
    seq[5] = b'A';
    seq[6] = b'A';

    let qual = TestData::qual(seq.len());

    assert!(detector.detect_first(&seq, &qual).unwrap().is_none());
}

#[test]
fn bd_v2_384_index_c1_accepts_one_base_fuzzy_rescue() {
    let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

    let exact = wl.index_c1(b"CGGAGAGAT").expect("exact C1 should exist");

    let rescued = wl
        .index_c1(b"CGGTGAGAT")
        .expect("C1 with one mismatch should be rescued");

    assert_eq!(rescued, exact);
}

#[test]
fn test_bd_v2_384_roundtrip_synthesize_detect() {
    let grammar = TestData::bd_v2_384_grammar("bd-roundtrip");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let cell = TestData::bd_v2_384_cell(0);
    let umi = b"AACCGG";
    let insert = b"GATCGATCGATC";

    let mut seq = grammar.synthesize(&cell, umi).unwrap();
    seq.extend_from_slice(insert);

    let qual = TestData::qual(seq.len());
    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.bd_cell_id, Some(1));
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().seq.as_slice(), umi);
    assert_eq!(hit.get_insert(&seq, &qual).unwrap().seq.as_slice(), insert);
}

#[test]
fn test_bd_v2_384_detection_returns_match_diagnostics() {
    let grammar = TestData::bd_v2_384_grammar("bd-diagnostics");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let cell = TestData::bd_v2_384_cell(0);
    let umi = b"AACCGG";
    let mut seq = grammar.synthesize(&cell, umi).unwrap();
    seq.extend_from_slice(b"GATCGATC");
    let qual = TestData::qual(seq.len());

    let (hit, diagnostics) = detector.detect_first_with_diagnostics(&seq, &qual).unwrap();
    let hit = hit.expect("BD primer should match");
    let bd = diagnostics
        .and_then(|diagnostics| diagnostics.bd)
        .expect("BD diagnostics should accompany the match");

    assert_eq!(hit.bd_cell_id, Some(1));
    assert_eq!(bd.shift, 0);
    assert_eq!(&bd.linker_signature, b"GTGAGACA");
    assert_eq!(bd.linker_mismatches, 0);
    assert_eq!(
        (bd.c1_mismatches, bd.c2_mismatches, bd.c3_mismatches),
        (0, 0, 0)
    );
}

#[test]
fn test_bd_v2_384_invalid_cell_is_translated_to_valid_cell() {
    let grammar = TestData::bd_v2_384_grammar("bd-repair");
    let mut detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let mut source_cell = TestData::bd_v2_384_cell(0);

    // Damage C1 enough that this exact cassette is not valid anymore.
    source_cell[0] = b'T';
    source_cell[1] = b'T';
    source_cell[2] = b'T';

    let record = ReadTagRecord {
        read_id: "read1".to_string(),
        original_read_id: None,
        cell_seq: source_cell.clone(),
        cell_qual: vec![b'I'; source_cell.len()],
        umi_seq: b"AACCGG".to_vec(),
        umi_qual: vec![b'I'; 6],
    };

    let (target_cell, mut primer) = detector.generate(&record).unwrap();
    primer.extend_from_slice(b"GATCGATCGATC");

    let qual = TestData::qual(primer.len());

    let hit = detector.detect_first(&primer, &qual).unwrap().unwrap();

    assert_ne!(target_cell, source_cell);
    assert_eq!(target_cell, TestData::bd_v2_384_cell(0));
    assert_eq!(hit.bd_cell_id, Some(1));

    let translation = detector.primer_translation();

    assert_eq!(translation.len(), 1);
}

#[test]
fn test_reverse_complement_detection_for_tenx_like_read() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-rc");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let forward = TestData::synthesized_tenx_read(&grammar);
    let reverse = TestData::rc(&forward);
    let qual = TestData::qual(reverse.len());

    let hit = detector.detect_first(&reverse, &qual).unwrap().unwrap();

    assert_eq!(hit.orientation, Orientation::ReverseComplement);
    assert!(hit.get_cell(&reverse, &qual).is_ok());
    assert!(hit.get_umi(&reverse, &qual).is_ok());
}

#[test]
fn test_detect_all_three_tenx_monomers_with_damaged_junctions() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-three");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let one = TestData::synthesized_tenx_read(&grammar);
    let two = TestData::synthesized_tenx_read(&grammar);
    let three = TestData::synthesized_tenx_read(&grammar);
    let seq = TestData::cat(&[&one, b"NNNCCT", &two, b"GGGGNN", &three]);
    let qual = TestData::qual(seq.len());

    let hits = detector.detect_all(&seq, &qual).unwrap();

    assert_eq!(hits.len(), 3);
    for hit in hits {
        assert_eq!(
            hit.get_cell(&seq, &qual).unwrap().seq.as_slice(),
            TestData::tenx_cell()
        );
    }
}

#[test]
fn test_detect_all_tenx_ont_read_with_ten_primer_insert_monomers() {
    TenxOntMultimerTest::run();
}

#[test]
fn test_quality_is_sliced_with_cell_and_umi() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-quality");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let seq = TestData::synthesized_tenx_read(&grammar);
    let mut qual = TestData::qual(seq.len());

    let fixed_len = b"CTACACGACGCTCTTCCGATCT".len();
    let cell_start = fixed_len;
    let umi_start = cell_start + 16;

    for q in qual.iter_mut().skip(cell_start).take(16) {
        *q = b'!';
    }
    for q in qual.iter_mut().skip(umi_start).take(12) {
        *q = b'#';
    }

    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();

    assert_eq!(hit.get_cell(&seq, &qual).unwrap().qual, vec![b'!'; 16]);
    assert_eq!(hit.get_umi(&seq, &qual).unwrap().qual, vec![b'#'; 12]);
}

#[test]
fn test_sequence_quality_length_mismatch_is_error() {
    let grammar = TestData::tenx_3p_v3_grammar("tenx-length-error");
    let detector = PrimerDetector::from_grammar(grammar.clone()).unwrap();

    let seq = TestData::synthesized_tenx_read(&grammar);
    let qual = TestData::qual(seq.len() - 1);

    assert!(detector.detect_all(&seq, &qual).is_err());
}

#[test]
fn test_bd_cell_version_parser() {
    assert_eq!(BdCellVersion::parse("v1").unwrap(), BdCellVersion::V1);
    assert_eq!(BdCellVersion::parse("v2.96").unwrap(), BdCellVersion::V2_96);
    assert_eq!(
        BdCellVersion::parse("v2.384").unwrap(),
        BdCellVersion::V2_384
    );
    assert!(BdCellVersion::parse("v3").is_err());
}

/*
#[test]
//#[ignore = "timing/stress test; run manually"]
fn stress_detect_all_tenx_v3_1000x() {
    use std::time::Instant;

    let detector = PrimerDetector::from_chemistry(Chemistry::TenxThreePrimeV3).unwrap();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    for i in 0..10 {
        let cell = detector
            .single_cell_system
            .as_ref()
            .unwrap()
            .cell_seq_for_index(i)
            .unwrap();

        let umi = detector.grammar().umi_from_u64(i as u64);

        let mut primer = detector
            .grammar()
            .synthesize(&cell, &umi)
            .unwrap();

        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");

        qual.extend(std::iter::repeat_n(b'I', primer.len()));
        seq.extend_from_slice(&primer);
    }

    let first = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(first.len(), 10);

    let start = Instant::now();

    for _ in 0..1000 {
        let hits = detector.detect_all(&seq, &qual).unwrap();
        assert_eq!(hits.len(), 10);
        std::hint::black_box(hits);
    }

    let elapsed = start.elapsed();

    eprintln!(
        "detect_all 10x-v3 stress: 1000 reads / 10000 primer matches in {:?}; {:.3} µs/read; {:.3} µs/match",
        elapsed,
        elapsed.as_secs_f64() * 1_000_000.0 / 1000.0,
        elapsed.as_secs_f64() * 1_000_000.0 / 10_000.0,
    );
}

#[test]
//#[ignore = "timing/stress test; run manually"]
fn stress_detect_all_bd_v2_384_1000x() {
    use std::time::Instant;

    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    let wl = RhapsodyWhitelist::bd_v2_384();

    for i in 0..10 {
        let cell_id = (i + 1) as u64;

        let (c1, c2, c3) = wl
            .cell_id_to_parts_ids(cell_id)
            .expect("invalid BD cell id");

        let cell = wl.create_cell_cassette(c1, c2, c3);
        let umi = detector.grammar().umi_from_u64(i as u64);

        let mut primer = detector
            .grammar()
            .synthesize(&cell, &umi)
            .unwrap();

        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");

        qual.extend(std::iter::repeat_n(b'I', primer.len()));
        seq.extend_from_slice(&primer);
    }

    let first = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(first.len(), 10);

    let start = Instant::now();

    for _ in 0..1000 {
        let hits = detector.detect_all(&seq, &qual).unwrap();
        assert_eq!(hits.len(), 10);
        std::hint::black_box(hits);
    }

    let elapsed = start.elapsed();

    eprintln!(
        "detect_all bd-v2-384 stress: 1000 reads / 10000 primer matches in {:?}; {:.3} µs/read; {:.3} µs/match",
        elapsed,
        elapsed.as_secs_f64() * 1_000_000.0 / 1000.0,
        elapsed.as_secs_f64() * 1_000_000.0 / 10_000.0,
    );
}*/

#[test]
fn benchmark_detect_all_tenx_positional() {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxThreePrimeV3).unwrap();
    let reads = TenxOntMultimerTest::build_chemistry_monomers(&detector, false);

    for (seq, qual) in &reads {
        let hits = detector.detect_all(seq, qual).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].primer_start, 0);
    }

    bench("10x detect_all positional", 10_000, reads.len(), || {
        for (seq, qual) in &reads {
            std::hint::black_box(detector.detect_all(seq, qual).unwrap());
        }
    });
}

#[test]
#[ignore = "benchmark"]
fn benchmark_detect_all_tenx_positional_long() {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxThreePrimeV3).unwrap();
    let reads = TenxOntMultimerTest::build_chemistry_monomers(&detector, false);

    for (seq, qual) in &reads {
        let hits = detector.detect_all(seq, qual).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].primer_start, 0);
    }

    bench(
        "10x detect_all positional [long]",
        100_000,
        reads.len(),
        || {
            for (seq, qual) in &reads {
                std::hint::black_box(detector.detect_all(seq, qual).unwrap());
            }
        },
    );
}

#[test]
fn benchmark_detect_all_bd_v2_384_multimer() {
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    let wl = RhapsodyWhitelist::bd_v2_384();

    for i in 0..100 {
        let cell_id = (i + 1) as u64;

        let (c1, c2, c3) = wl.cell_id_to_parts_ids(cell_id).expect("invalid cell id");

        let cell = wl.create_cell_cassette(c1, c2, c3);

        let umi = detector.grammar().umi_from_u64(i as u64);

        let mut primer = detector.grammar().synthesize(&cell, &umi).unwrap();

        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");

        seq.extend_from_slice(&primer);
        qual.extend(std::iter::repeat_n(b'I', primer.len()));
    }

    let hits = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].bd_cell_id, Some(1));

    bench("BD detect_all positional multimer", 10_000, 1, || {
        std::hint::black_box(detector.detect_all(&seq, &qual).unwrap());
    });
}

#[test]
#[ignore = "benchmark"]
fn benchmark_detect_all_bd_v2_384_multimer_long() {
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    let wl = RhapsodyWhitelist::bd_v2_384();

    for i in 0..100 {
        let cell_id = (i + 1) as u64;
        let (c1, c2, c3) = wl.cell_id_to_parts_ids(cell_id).expect("invalid cell id");
        let cell = wl.create_cell_cassette(c1, c2, c3);
        let umi = detector.grammar().umi_from_u64(i as u64);
        let mut primer = detector.grammar().synthesize(&cell, &umi).unwrap();

        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");
        seq.extend_from_slice(&primer);
        qual.extend(std::iter::repeat_n(b'I', primer.len()));
    }

    let hits = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].bd_cell_id, Some(1));

    bench("BD detect_all positional multimer [long]", 100_000, 1, || {
        std::hint::black_box(detector.detect_all(&seq, &qual).unwrap());
    });
}

#[test]
fn benchmark_detect_all_tenx_positional_fuzzy() {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxThreePrimeV3).unwrap();
    let reads = TenxOntMultimerTest::build_chemistry_monomers(&detector, true);

    for (seq, qual) in &reads {
        let hits = detector.detect_all(seq, qual).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].primer_start, 0);
    }

    bench(
        "10x fuzzy detect_all positional",
        10_000,
        reads.len(),
        || {
            for (seq, qual) in &reads {
                std::hint::black_box(detector.detect_all(seq, qual).unwrap());
            }
        },
    );
}

#[test]
#[ignore = "benchmark"]
fn benchmark_detect_all_tenx_positional_fuzzy_long() {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxThreePrimeV3).unwrap();
    let reads = TenxOntMultimerTest::build_chemistry_monomers(&detector, true);

    for (seq, qual) in &reads {
        let hits = detector.detect_all(seq, qual).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].primer_start, 0);
    }

    bench(
        "10x fuzzy detect_all positional [long]",
        100_000,
        reads.len(),
        || {
            for (seq, qual) in &reads {
                std::hint::black_box(detector.detect_all(seq, qual).unwrap());
            }
        },
    );
}

#[test]
fn benchmark_detect_all_bd_v2_384_multimer_fuzzy() {
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    let wl = RhapsodyWhitelist::bd_v2_384();
    let hits_per_read = 100;

    for i in 0..hits_per_read {
        let cell_id = (i + 1) as u64;

        let (c1, c2, c3) = wl
            .cell_id_to_parts_ids(cell_id)
            .expect("invalid BD cell id");

        let mut cell = wl.create_cell_cassette(c1, c2, c3);

        // One damaged base in C1. This should force fuzzy rescue.
        cell[0] = b'N';

        let umi = detector.grammar().umi_from_u64(i as u64);

        let mut primer = detector.grammar().synthesize(&cell, &umi).unwrap();

        // Grammar::synthesize() reserves four bases for SEARCH before the
        // BD cassette. Damage one base in each linker at the actual cassette
        // coordinates, in addition to the damaged barcode base above.
        // Layout after SEARCH(4): C1(9) + L1(4) + C2(9) + L2(4) + C3(9).
        const SEARCH_PREFIX: usize = 4;
        primer[SEARCH_PREFIX + 9 + 2] = b'A';
        primer[SEARCH_PREFIX + 9 + 4 + 9 + 2] = b'T';

        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");

        qual.extend(std::iter::repeat_n(b'I', primer.len()));
        seq.extend_from_slice(&primer);
    }

    let hits = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].bd_cell_id, Some(1));

    bench(
        "BD fuzzy detect_all positional multimer",
        10_000,
        1,
        || {
            std::hint::black_box(detector.detect_all(&seq, &qual).unwrap());
        },
    );
}
#[test]
#[ignore = "benchmark"]
fn benchmark_detect_all_bd_v2_384_multimer_fuzzy_long() {
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    let wl = RhapsodyWhitelist::bd_v2_384();
    let hits_per_read = 100;

    for i in 0..hits_per_read {
        let cell_id = (i + 1) as u64;
        let (c1, c2, c3) = wl
            .cell_id_to_parts_ids(cell_id)
            .expect("invalid BD cell id");
        let mut cell = wl.create_cell_cassette(c1, c2, c3);

        cell[0] = b'N';

        let umi = detector.grammar().umi_from_u64(i as u64);
        let mut primer = detector.grammar().synthesize(&cell, &umi).unwrap();

        // Grammar::synthesize() reserves four bases for SEARCH before the
        // BD cassette. Damage one base in each linker at the actual cassette
        // coordinates, in addition to the damaged barcode base above.
        // Layout after SEARCH(4): C1(9) + L1(4) + C2(9) + L2(4) + C3(9).
        const SEARCH_PREFIX: usize = 4;
        primer[SEARCH_PREFIX + 9 + 2] = b'A';
        primer[SEARCH_PREFIX + 9 + 4 + 9 + 2] = b'T';

        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");
        qual.extend(std::iter::repeat_n(b'I', primer.len()));
        seq.extend_from_slice(&primer);
    }

    let hits = detector.detect_all(&seq, &qual).unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].bd_cell_id, Some(1));

    bench(
        "BD fuzzy detect_all positional multimer [long]",
        100_000,
        1,
        || {
            std::hint::black_box(detector.detect_all(&seq, &qual).unwrap());
        },
    );
}

#[test]
fn tenx_builtin_is_positional_and_does_not_scan_for_a_later_barcode() {
    let detector = PrimerDetector::from_chemistry(Chemistry::TenxThreePrimeV3).unwrap();
    let mut reads = TenxOntMultimerTest::build_chemistry_monomers(&detector, false);
    let (valid_seq, valid_qual) = reads.remove(0);

    let mut shifted_seq = b"GATCGATC".to_vec();
    shifted_seq.extend_from_slice(&valid_seq);
    let mut shifted_qual = vec![b'I'; 8];
    shifted_qual.extend_from_slice(&valid_qual);

    assert!(detector
        .detect_all(&shifted_seq, &shifted_qual)
        .unwrap()
        .is_empty());
}

#[test]
fn bd_v2_384_detect_all_does_not_repeat_search_window_hits() {
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384).unwrap();
    let wl = RhapsodyWhitelist::bd_v2_384();

    let mut seq = Vec::new();
    let mut qual = Vec::new();

    for index in 0..2u64 {
        let (c1, c2, c3) = wl
            .cell_id_to_parts_ids(index + 1)
            .expect("invalid BD cell id");
        let cell = wl.create_cell_cassette(c1, c2, c3);
        let umi = detector.grammar().umi_from_u64(index);

        let mut primer = detector.grammar().synthesize(&cell, &umi).unwrap();
        primer.extend_from_slice(b"GATCGATCGATCGATCGATCGATCGATCG");

        qual.extend(std::iter::repeat_n(b'I', primer.len()));
        seq.extend_from_slice(&primer);
    }

    let hits = detector.detect_all(&seq, &qual).unwrap();

    assert_eq!(
        hits.len(),
        1,
        "BD detection must stay inside the leading SEARCH window"
    );
    assert_eq!(hits[0].bd_cell_id, Some(1));
    assert!(hits[0].primer_end > hits[0].primer_start);
}

#[test]
fn molecule_identity_preserves_barcoded_hard_umi_rule() {
    use int_to_str::IntToStr;

    let grammar = Grammar::parse("identity", "CELL:4+UMI:4").unwrap();
    let cell = b"ACGT";
    let umi = b"TGCA";
    let r2 = b"AACCGGTTAACCGGTTAACCGGTTAACCGGTT";

    let identity = grammar
        .molecule_identity(Some(cell), Some(umi), b"IGNORED_R1", r2)
        .unwrap();

    let old_cell = IntToStr::new(cell).into_u64();
    let mut old_hard = umi.to_vec();
    old_hard.extend_from_slice(&r2[..28]);
    let old_hard = IntToStr::new(&old_hard).into_u64();

    assert_eq!(identity.cell_id, old_cell);
    assert_eq!(identity.molecule_id, old_hard);
}

#[test]
fn none_and_insert_only_grammars_use_r1_and_r2_for_molecule_identity() {
    let none = Grammar::parse("none", "NONE").unwrap();
    let insert_only = Grammar::parse("insert", "INSERT:ACGT").unwrap();

    assert!(none.is_unbarcoded());
    assert!(insert_only.is_unbarcoded());

    let r1 = b"AAAACCCCGGGGTTTTAAAA";
    let r2 = b"TTTTGGGGCCCCAAAATTTT";

    let none_id = none.molecule_identity(None, None, r1, r2).unwrap();
    let insert_id = insert_only.molecule_identity(None, None, r1, r2).unwrap();
    assert_eq!(none_id, insert_id);

    let changed_r1 = b"CAAACCCCGGGGTTTTAAAA";
    assert_ne!(
        none_id.molecule_id,
        none.molecule_identity(None, None, changed_r1, r2)
            .unwrap()
            .molecule_id
    );

    let changed_r2 = b"CTTTGGGGCCCCAAAATTTT";
    assert_ne!(
        none_id.molecule_id,
        none.molecule_identity(None, None, r1, changed_r2)
            .unwrap()
            .molecule_id
    );
}

#[test]
fn none_grammar_treats_the_whole_read_as_insert() {
    let grammar = Grammar::parse("none", "NONE").unwrap();
    let detector = PrimerDetector::from_grammar(grammar).unwrap();
    let seq = b"ACGTACGTACGTACGTACGT";
    let qual = vec![b'I'; seq.len()];

    let hit = detector.detect_first(seq, &qual).unwrap().unwrap();
    let insert = hit.get_insert(seq, &qual).unwrap();

    assert_eq!(insert.seq, seq);
    assert_eq!(hit.primer_start, 0);
    assert_eq!(hit.primer_end, 0);
    assert_eq!(hit.insert_start, 0);
    assert_eq!(hit.insert_end, seq.len());
}

#[test]
fn bd_vdj_uses_its_own_cc_handle_and_linkers() {
    let detector = PrimerDetector::from_chemistry(Chemistry::BdV2_384Vdj).unwrap();
    let wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384Vdj);
    let cell = wl.cell_id_to_cassette(1).unwrap();
    let mut seq = detector.grammar().synthesize(&cell, b"ACGTAC").unwrap();
    seq.extend_from_slice(b"TATGCGTAGTAGGTATGACGTACGTACGT");
    let qual = vec![40u8; seq.len()];

    let hit = detector.detect_first(&seq, &qual).unwrap().unwrap();
    assert_eq!(hit.chemistry_name, "bd-v2-384-vdj");
    assert_eq!(hit.bd_cell_id, Some(1));
    assert_eq!(hit.cell_seq.as_deref(), wl.cell_id_to_seq(1).as_deref());
}

#[test]
fn multiple_chemistries_prefer_fast_fixed_anchor_but_detect_both_bd_structures() {
    let detector = PrimerDetector::from_chemistries([
        Chemistry::BdV2_384,
        Chemistry::BdV2_384Vdj,
    ])
    .unwrap();

    let names = detector.grammars().map(|g| g.name.as_str()).collect::<Vec<_>>();
    assert_eq!(names, vec!["bd-v2-384-vdj", "bd-v2-384"]);

    let vdj_wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384Vdj);
    let vdj_cell = vdj_wl.cell_id_to_cassette(1).unwrap();
    let vdj_grammar = detector
        .grammars()
        .find(|g| g.name == "bd-v2-384-vdj")
        .unwrap();
    let mut vdj = vdj_grammar.synthesize(&vdj_cell, b"ACGTAC").unwrap();
    vdj.extend_from_slice(b"TATGCGTAGTAGGTATGACGTACGT");
    let vdj_hit = detector.detect_first(&vdj, &vec![40u8; vdj.len()]).unwrap().unwrap();
    assert_eq!(vdj_hit.chemistry_name, "bd-v2-384-vdj");

    let dt_wl = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);
    let dt_cell = dt_wl.cell_id_to_cassette(2).unwrap();
    let dt_grammar = detector
        .grammars()
        .find(|g| g.name == "bd-v2-384")
        .unwrap();
    let mut dt = dt_grammar.synthesize(&dt_cell, b"TGCATG").unwrap();
    dt.extend_from_slice(b"ACGTACGTACGTACGT");
    let dt_hit = detector.detect_first(&dt, &vec![40u8; dt.len()]).unwrap().unwrap();
    assert_eq!(dt_hit.chemistry_name, "bd-v2-384");
}
