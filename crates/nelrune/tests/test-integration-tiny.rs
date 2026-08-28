use std::path::{Path, PathBuf};
use std::process::Command;

use bam_tide::index::GeneFeatureIndex;
use gtf_splice_index::{
    AnnotationBuilder,
    SpliceIndex,
};
use scdata::{
    feature_index::FeatureIndex,
    Scdata,
};

const TEST_DATA: &str = "tests/data";

fn test_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join(TEST_DATA)
        .join(name)
}

fn require_file(path: &Path) {
    assert!(
        path.is_file(),
        "required test input is missing: {}",
        path.display()
    );
}

fn star_available() -> bool {
    Command::new("STAR")
        .arg("--version")
        .output()
        .map(|out| out.status.success())
        .unwrap_or(false)
}

fn ensure_star_index(
    fasta: &Path,
    gtf: &Path,
    star_index: &Path,
) {
    if star_index.join("Genome").is_file() {
        return;
    }

    assert!(
        star_available(),
        "STAR is required for this integration test but was not found in PATH"
    );

    std::fs::create_dir_all(star_index)
        .expect("failed to create STAR index directory");

    let status = Command::new("STAR")
        .args([
            "--runMode",
            "genomeGenerate",

            "--genomeDir",
            star_index
                .to_str()
                .unwrap(),

            "--genomeFastaFiles",
            fasta
                .to_str()
                .unwrap(),

            "--sjdbGTFfile",
            gtf
                .to_str()
                .unwrap(),

            // Tiny artificial genome.
            "--genomeSAindexNbases",
            "2",

            // Our reads are 60 bp.
            "--sjdbOverhang",
            "59",
        ])
        .status()
        .expect("failed to start STAR genomeGenerate");

    assert!(
        status.success(),
        "STAR genomeGenerate failed"
    );

    assert!(
        star_index.join("Genome").is_file(),
        "STAR finished successfully but no Genome file was created"
    );
}

fn ensure_splice_index(
    gtf: &Path,
    index_path: &Path,
) {
    if index_path.is_file() {
        return;
    }

    let index =
        AnnotationBuilder::new(1_000)
            .build_from_path(gtf)
            .expect("failed to build tiny splice index");

    index
        .save(index_path)
        .expect("failed to write tiny splice index");

    assert!(
        index_path.is_file(),
        "splice index was not created"
    );
}

fn assert_matrix_files(path: &Path) {
    for name in [
        "matrix.mtx.gz",
        "barcodes.tsv.gz",
        "features.tsv.gz",
    ] {
        assert!(
            path.join(name).is_file(),
            "missing expected matrix output: {}",
            path.join(name).display()
        );
    }
}

#[test]
fn integration_tiny_star() {
    // --------------------------------------------------------
    // 1. Test inputs
    // --------------------------------------------------------

    let fasta =
        test_path("tiny.fa");

    let gtf =
        test_path("tiny.gtf");

    let r1 =
        test_path("tiny_R1.fastq");

    let r2 =
        test_path("tiny_R2.fastq");

    for path in [
        &fasta,
        &gtf,
        &r1,
        &r2,
    ] {
        require_file(path);
    }

    // --------------------------------------------------------
    // 2. STAR index
    // --------------------------------------------------------

    let star_index =
        test_path("star_index");

    ensure_star_index(
        &fasta,
        &gtf,
        &star_index,
    );

    // --------------------------------------------------------
    // 3. gtf-splice-index
    // --------------------------------------------------------

    let splice_index_path =
        test_path("tiny.gtf.dat");

    ensure_splice_index(
        &gtf,
        &splice_index_path,
    );

    // Make sure it can actually be loaded again.
    let splice_index =
        SpliceIndex::load(
            &splice_index_path
        )
        .expect(
            "created splice index cannot be loaded"
        );

    assert_eq!(
        splice_index.genes.len(),
        2,
        "tiny annotation should contain exactly two genes"
    );

    assert_eq!(
        splice_index.transcripts.len(),
        2,
        "tiny annotation should contain exactly two transcripts"
    );

    // --------------------------------------------------------
    // 4. Run the actual Nelrune executable
    // --------------------------------------------------------

    let out =
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("target/test-integration-tiny");

    let _ =
        std::fs::remove_dir_all(&out);

    std::fs::create_dir_all(&out)
        .expect(
            "failed to create test output directory"
        );

    let output =
        Command::new(
            env!("CARGO_BIN_EXE_nelrune")
        )
        .args([
            "--r1",
            r1.to_str().unwrap(),

            "--r2",
            r2.to_str().unwrap(),

            "--primer-structure",
            "CELL:4+UMI:4",

            "--mapper",
            "star",

            "--mapper-index",
            star_index.to_str().unwrap(),

            "--mapper-threads",
            "2",

            "--index",
            splice_index_path
                .to_str()
                .unwrap(),

            "--require-strand",

            "--min-mapq",
            "0",

            // Two reads should survive cell filtering.
            "--min-cell-counts",
            "1",

            "--min-insert-len",
            "20",

            "--threads",
            "2",

            "--outpath",
            out.to_str().unwrap(),

            // Server behaviour deserves its own test.
            "--no-health-server",
        ])
        .output()
        .expect(
            "failed to start nelrune"
        );

    if !output.status.success() {
        panic!(
            "\nNelrune integration test failed.\n\
             status: {}\n\
             stdout:\n{}\n\
             stderr:\n{}\n",
            output.status,
            String::from_utf8_lossy(
                &output.stdout
            ),
            String::from_utf8_lossy(
                &output.stderr
            ),
        );
    }

    // --------------------------------------------------------
    // 5. Full output stack exists
    // --------------------------------------------------------

    let exonic =
        out.join("exonic");

    let intronic =
        out.join("intronic");

    assert_matrix_files(&exonic);
    assert_matrix_files(&intronic);

    assert!(
        out.join("nelrune.log").is_file(),
        "Nelrune log is missing"
    );

    assert!(
        out.join("nelrune-report.txt").is_file(),
        "Nelrune final report is missing"
    );

    // No SNP fixture in this test, so ref/alt are deliberately
    // not asserted here. That belongs in the SNP side-channel
    // integration test.

    // --------------------------------------------------------
    // 6. Read results back through Scdata
    // --------------------------------------------------------

    let gene_index =
        GeneFeatureIndex::new(
            &splice_index
        );

    let plus_gene =
        gene_index
            .feature_id("gene_plus")
            .expect(
                "gene_plus missing from feature index"
            );

    let minus_gene =
        gene_index
            .feature_id("gene_minus")
            .expect(
                "gene_minus missing from feature index"
            );

    let exonic_data =
        Scdata::read_matrix_market(
            &exonic,
            &gene_index,
        )
        .expect(
            "failed to reload exonic Nelrune matrix"
        );

    let intronic_data =
        Scdata::read_matrix_market(
            &intronic,
            &gene_index,
        )
        .expect(
            "failed to reload intronic Nelrune matrix"
        );

    // --------------------------------------------------------
    // Ground truth
    //
    // R2 read 1:
    //   exon1 -> exon2
    //   therefore gene_plus EXONIC
    //
    // R2 read 2:
    //   exon1 -> intron
    //   therefore gene_plus INTRONIC
    //
    // Both sequences are plus-strand genomic orientation.
    // Because --require-strand is enabled, gene_minus must
    // receive nothing despite sharing exactly the same exons.
    // --------------------------------------------------------

    assert_eq!(
        exonic_data.cell_ids().len(),
        1,
        "expected exactly one cell in exonic output"
    );

    assert_eq!(
        intronic_data.cell_ids().len(),
        1,
        "expected exactly one cell in intronic output"
    );

    let cell_id =
        *exonic_data
            .cell_ids()
            .iter()
            .next()
            .unwrap();

    assert!(
        intronic_data
            .cell_ids()
            .contains(&cell_id),
        "exonic and intronic reads should belong to the same cell"
    );

    let exonic_cell =
        exonic_data
            .get(&cell_id)
            .expect(
                "cell missing from exonic Scdata"
            );

    let intronic_cell =
        intronic_data
            .get(&cell_id)
            .expect(
                "cell missing from intronic Scdata"
            );

    assert_eq!(
        exonic_cell
            .total_umis_4_gene_id(
                &plus_gene
            ),
        1.0,
        "spliced read should produce exactly one gene_plus exonic UMI"
    );

    assert_eq!(
        intronic_cell
            .total_umis_4_gene_id(
                &plus_gene
            ),
        1.0,
        "unspliced read should produce exactly one gene_plus intronic UMI"
    );

    assert_eq!(
        exonic_cell
            .total_umis_4_gene_id(
                &minus_gene
            ),
        0.0,
        "plus-strand spliced read must not be assigned to gene_minus"
    );

    assert_eq!(
        intronic_cell
            .total_umis_4_gene_id(
                &minus_gene
            ),
        0.0,
        "plus-strand intronic read must not be assigned to gene_minus"
    );

    // Nice high-level sanity check too.
    assert_eq!(
        exonic_cell.total_umis(),
        1
    );

    assert_eq!(
        intronic_cell.total_umis(),
        1
    );
}
