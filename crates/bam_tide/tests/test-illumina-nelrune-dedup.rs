// crates/bam_tide/tests/test-illumina-nelrune-dedup.rs

use anyhow::Result;
use tempfile::tempdir;

use bam_tide::illumina_normalizer::{IlluminaNormalizer, IlluminaNormalizerConfig};

use bam_tide::illumina_normalizer::cli::{InsertRead, PrimerRead};

use sc_primer::{Grammar, PrimerDetector};

#[test]
fn nelrune_run_does_not_emit_pcr_duplicates() -> Result<()> {
    let dir = tempdir()?;

    let r1_path = dir.path().join("R1.fastq");

    let r2_path = dir.path().join("R2.fastq");

    // --------------------------------------------------------
    // Two reads:
    //
    // same cell
    // same UMI
    // same biological insert
    //
    // Only the FASTQ read names differ.
    //
    // Therefore they must generate the same DedupKey and only
    // ONE mapper-facing molecule may leave nelrune_run().
    // --------------------------------------------------------

    std::fs::write(
        &r1_path,
        "\
@read1
ACGTAAAA
+
IIIIIIII
@read2
ACGTAAAA
+
IIIIIIII
",
    )?;

    std::fs::write(
        &r2_path,
        "\
@read1
TACGCTAGCATGCTACGATCGTAGCTACGAATGCTACGTAGCTACGATCGTAGCTAGCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
TACGCTAGCATGCTACGATCGTAGCTACGAATGCTACGTAGCTACGATCGTAGCTAGCAT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
",
    )?;

    // --------------------------------------------------------
    // Minimal artificial primer grammar
    //
    // R1:
    //
    // ACGT AAAA
    // ---- ----
    // cell UMI
    // --------------------------------------------------------

    let grammar = Grammar::parse("dedup-test", "CELL:4+UMI:4").map_err(anyhow::Error::msg)?;

    let primer = PrimerDetector::from_grammar(grammar).map_err(anyhow::Error::msg)?;

    let config = IlluminaNormalizerConfig {
        /*
         * nelrune_run() does not write these files.
         */
        out: dir.path().join("unused.fastq"),

        read_tags: dir.path().join("unused.tags"),

        primer_read: PrimerRead::R1,

        insert_read: InsertRead::R2,

        primer,

        additional_features: Vec::new(),

        additional_feature_min_hits: 4,

        min_insert_len: 20,

        threads: 1,

        gzip_level: 1,

        max_reads: Some(10),

        gzip: false,
    };

    let mut normalizer = IlluminaNormalizer::new(config)?;

    // --------------------------------------------------------
    // This Vec represents exactly what Nelrune would submit to
    // STAR.
    // --------------------------------------------------------

    let mut emitted = Vec::new();

    normalizer.nelrune_run(
        &r1_path,
        &r2_path,
        |batch| {
            emitted.extend(
                batch
                    .iter()
                    .map(|(r1, r2)| (r1.as_ref().map(|read| read.id.clone()), r2.id.clone())),
            );

            Ok(true)
        },
        |_| {},
    )?;

    // --------------------------------------------------------
    // The critical contract:
    //
    // only one PCR-equivalent molecule is allowed through the
    // Nelrune emission boundary.
    // --------------------------------------------------------

    assert_eq!(
        emitted.len(),
        1,
        "PCR-equivalent reads must only be emitted once to Nelrune/STAR"
    );

    Ok(())
}

#[test]
fn nelrune_run_none_grammar_dedups_by_r1_and_r2_sequence() -> Result<()> {
    let dir = tempdir()?;
    let r1_path = dir.path().join("R1.none.fastq");
    let r2_path = dir.path().join("R2.none.fastq");

    // pair 1 + pair 2: identical sequences, different QNAME => PCR duplicate
    // pair 3: same R1 but changed R2 prefix => distinct molecule
    std::fs::write(
        &r1_path,
        "\
@read1\nAAAACCCCGGGGTTTTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nAAAACCCCGGGGTTTTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read3\nAAAACCCCGGGGTTTTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )?;

    std::fs::write(
        &r2_path,
        "\
@read1\nTTTTGGGGCCCCAAAAACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read2\nTTTTGGGGCCCCAAAAACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@read3\nCTTTGGGGCCCCAAAAACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )?;

    let grammar = Grammar::parse("none", "NONE").map_err(anyhow::Error::msg)?;
    let primer = PrimerDetector::from_grammar(grammar).map_err(anyhow::Error::msg)?;

    let config = IlluminaNormalizerConfig {
        out: dir.path().join("unused.fastq"),
        read_tags: dir.path().join("unused.tags"),
        primer_read: PrimerRead::R1,
        insert_read: InsertRead::R2,
        primer,
        additional_features: Vec::new(),
        additional_feature_min_hits: 4,
        min_insert_len: 20,
        threads: 1,
        gzip_level: 1,
        max_reads: Some(10),
        gzip: false,
    };

    let mut normalizer = IlluminaNormalizer::new(config)?;
    let mut emitted = Vec::new();

    normalizer.nelrune_run(
        &r1_path,
        &r2_path,
        |batch| {
            emitted.extend(batch.iter().map(|(_, r2)| r2.id.clone()));
            Ok(true)
        },
        |_| {},
    )?;

    assert_eq!(
        emitted.len(),
        2,
        "NONE grammar must collapse identical R1+R2 pairs but retain a pair whose R2 sequence differs"
    );

    Ok(())
}

#[test]
fn bd_cell_umi_primer_qname_round_trip_preserves_identity() -> Result<()> {
    use std::collections::HashMap;

    use sc_primer::{BdCellVersion, Chemistry, ReadTagRecord, RhapsodyWhitelist};

    let detector =
        PrimerDetector::from_chemistry(Chemistry::BdV2_384).map_err(anyhow::Error::msg)?;
    let grammar = detector.grammar();
    let whitelist = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);
    let r2 = b"ACGTGTCAGTACGATCGTAGCTAGCATCGATGCTAGCTACGATCGTAGCTAGCATCGATGCTAG";

    let umis: [&[u8]; 10] = [
        b"AAAAAA", b"AAAAAC", b"AAAAAG", b"AAAAAT", b"AAAACA", b"AAAACC", b"AAAACG", b"AAAACT",
        b"AAAAGA", b"AAAAGC",
    ];

    let mut recovered: HashMap<(Vec<u8>, Vec<u8>), usize> = HashMap::new();

    for cell_index in 0..10_u64 {
        let cell_id = cell_index + 1;
        let cell_cassette = detector
            .cell_seq_for_index(cell_index)
            .map_err(anyhow::Error::msg)?;
        let canonical_cell = whitelist
            .cell_id_to_seq(cell_id)
            .expect("synthetic BD cell id must have a canonical 27-base barcode");

        for (umi_index, umi) in umis.iter().enumerate() {
            for copy in 0..2 {
                let primer = grammar
                    .synthesize(&cell_cassette, umi)
                    .map_err(anyhow::Error::msg)?;
                let qual = vec![b'I'; primer.len()];
                let hit = detector
                    .detect_first(&primer, &qual)
                    .map_err(anyhow::Error::msg)?
                    .unwrap_or_else(|| {
                        panic!(
                            "primer detection failed for cell_index={cell_index} umi_index={umi_index} copy={copy} cell={} umi={}",
                            String::from_utf8_lossy(&canonical_cell),
                            String::from_utf8_lossy(umi),
                        )
                    });

                let detected_cell = hit.cell_seq.clone().unwrap_or_else(|| {
                    hit.get_cell(&primer, &qual)
                        .expect("synthetic primer CELL slice must be valid")
                        .seq
                        .to_vec()
                });
                let detected_umi = hit
                    .get_umi(&primer, &qual)
                    .expect("synthetic primer UMI slice must be valid")
                    .seq
                    .to_vec();

                assert_eq!(
                    detected_cell, canonical_cell,
                    "primer CELL round trip changed identity at cell_index={cell_index} cell_id={cell_id} umi_index={umi_index} copy={copy}"
                );
                assert_eq!(
                    detected_umi, *umi,
                    "primer UMI round trip changed identity at cell_index={cell_index} umi_index={umi_index} copy={copy}"
                );

                let before_identity = grammar
                    .molecule_identity_if_exact(
                        Some(&detected_cell),
                        Some(&detected_umi),
                        &primer,
                        r2,
                    )
                    .map_err(anyhow::Error::msg)?
                    .expect("synthetic A/C/G/T molecule identity must be exact");

                let tag = ReadTagRecord::new(
                    format!("synthetic-{cell_index}-{umi_index}-{copy}"),
                    None,
                    &detected_cell,
                    vec![b'I'; detected_cell.len()],
                    &detected_umi,
                    vec![b'I'; detected_umi.len()],
                );
                let qname = tag.extend_qname(&tag.read_id);
                let decoded = ReadTagRecord::from_qname(&qname)?;

                assert_eq!(
                    decoded.cell_seq, canonical_cell,
                    "QNAME hex CELL round trip changed identity at cell_index={cell_index} cell_id={cell_id} umi_index={umi_index} copy={copy}; qname={qname}"
                );
                assert_eq!(
                    decoded.umi_seq, *umi,
                    "QNAME hex UMI round trip changed identity at cell_index={cell_index} umi_index={umi_index} copy={copy}; qname={qname}"
                );

                let after_identity = grammar
                    .molecule_identity_if_exact(
                        Some(&decoded.cell_seq),
                        Some(&decoded.umi_seq),
                        &primer,
                        r2,
                    )
                    .map_err(anyhow::Error::msg)?
                    .expect("QNAME-decoded molecule identity must remain exact");

                assert_eq!(
                    after_identity, before_identity,
                    "numeric molecule identity changed across primer/QNAME round trip at cell_index={cell_index} umi_index={umi_index} copy={copy}"
                );

                *recovered
                    .entry((decoded.cell_seq, decoded.umi_seq))
                    .or_insert(0) += 1;
            }
        }
    }

    assert_eq!(
        recovered.len(),
        100,
        "10 cells x 10 UMIs must yield 100 identities"
    );
    assert!(
        recovered.values().all(|&count| count == 2),
        "every synthetic cell+UMI identity must survive exactly twice; bad counts: {:?}",
        recovered
            .iter()
            .filter(|&(_, &count)| count != 2)
            .collect::<Vec<_>>()
    );

    Ok(())
}
