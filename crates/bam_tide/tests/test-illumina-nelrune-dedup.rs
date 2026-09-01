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
