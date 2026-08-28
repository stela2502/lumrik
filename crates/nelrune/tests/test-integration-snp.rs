use std::path::{Path, PathBuf};
use std::process::Command;

use bam_tide::index::GeneFeatureIndex;

use gtf_splice_index::{
    AnnotationBuilder,
    SpliceIndex,
};

use scdata::feature_index::FeatureIndex;
use scdata::Scdata;

use snp_index::{
    SnpIndex,
    VcfReadOptions,
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

    std::fs::create_dir_all(star_index)
        .expect(
            "failed to create STAR index directory"
        );

    let status =
        Command::new("STAR")
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

                // chrTiny is deliberately tiny.
                "--genomeSAindexNbases",
                "2",

                // Our test reads are 60 bp.
                "--sjdbOverhang",
                "59",
            ])
            .status()
            .expect(
                "failed to start STAR genomeGenerate"
            );

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
            .expect(
                "failed to build tiny splice index"
            );

    index
        .save(index_path)
        .expect(
            "failed to write tiny splice index"
        );

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
        let file =
            path.join(name);

        assert!(
            file.is_file(),
            "missing expected matrix output: {}",
            file.display()
        );
    }
}


#[test]
fn integration_tiny_star_snp() {
    // --------------------------------------------------------
    // 1. Test inputs
    // --------------------------------------------------------

    let fasta =
        test_path("tiny.fa");

    let gtf =
        test_path("tiny.gtf");

    let vcf =
        test_path("tiny.vcf");

    let r1 =
        test_path("tiny_R1_snp.fastq");

    let r2 =
        test_path("tiny_R2_snp.fastq");

    for path in [
        &fasta,
        &gtf,
        &vcf,
        &r1,
        &r2,
    ] {
        require_file(path);
    }

    assert!(
        star_available(),
        "STAR is required for this integration test but was not found in PATH"
    );


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
        PathBuf::from(
            env!("CARGO_MANIFEST_DIR")
        )
        .join(
            "target/test-integration-tiny-snp"
        );

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
            star_index
                .to_str()
                .unwrap(),

            "--mapper-threads",
            "2",

            "--index",
            splice_index_path
                .to_str()
                .unwrap(),

            // Needed by the SNP side channel.
            "--genome",
            fasta.to_str().unwrap(),

            "--vcf",
            vcf.to_str().unwrap(),

            "--require-strand",

            "--min-mapq",
            "0",

            // One cell with four distinct molecules.
            "--min-cell-counts",
            "1",

            "--min-insert-len",
            "20",

            "--threads",
            "2",

            "--outpath",
            out.to_str().unwrap(),

            // Health-server behaviour has its own responsibility;
            // don't make this integration test depend on a free port.
            "--no-health-server",
        ])
        .output()
        .expect(
            "failed to start nelrune"
        );

    if !output.status.success() {
        panic!(
            "\nNelrune SNP integration test failed.\n\
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

    let snp_ref =
        out.join("ref");

    let snp_alt =
        out.join("alt");

    assert_matrix_files(
        &exonic
    );

    assert_matrix_files(
        &intronic
    );

    assert_matrix_files(
        &snp_ref
    );

    assert_matrix_files(
        &snp_alt
    );

    assert!(
        out.join("nelrune.log").is_file(),
        "Nelrune log is missing"
    );

    assert!(
        out.join("nelrune-report.txt")
            .is_file(),
        "Nelrune final report is missing"
    );


    // --------------------------------------------------------
    // 6. Build exactly the feature indexes used to interpret
    //    the matrices.
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


    // --------------------------------------------------------
    // SNP index
    //
    // This must describe the same reference coordinate space
    // that Nelrune saw from STAR:
    //
    //     chrTiny, length 400
    //
    // tiny.vcf contains exactly:
    //
    //     chrTiny:135 A>C
    //
    // --------------------------------------------------------

    let snp_index =
        SnpIndex::from_vcf_path(
            &vcf,
            vec![
                "chrTiny".to_string()
            ],
            vec![
                400u32
            ],
            10_000,
            &VcfReadOptions::default(),
        )
        .expect(
            "failed to load tiny SNP index"
        );

    let tiny_snp =
        snp_index
            .feature_id("tiny_snp")
            .expect(
                "tiny_snp missing from SNP index"
            );


    // --------------------------------------------------------
    // 7. Reload gene matrices through Scdata
    // --------------------------------------------------------

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
    // Ground truth for gene quantification
    //
    // All four reads belong to CELL ACGT.
    //
    // spliced_REF:
    //     gene_plus exonic
    //
    // spliced_ALT:
    //     gene_plus exonic
    //
    // unspliced_REF:
    //     gene_plus intronic
    //
    // unspliced_ALT:
    //     gene_plus intronic
    //
    // All four have distinct UMIs.
    //
    // Because --require-strand is enabled, gene_minus must
    // receive nothing despite occupying exactly the same
    // genomic exon coordinates.
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

    let exonic_cell_ids =
        exonic_data.cell_ids();

    let cell_id =
        *exonic_cell_ids
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
        2.0,
        "REF + ALT spliced reads should produce two gene_plus exonic UMIs"
    );

    assert_eq!(
        intronic_cell
            .total_umis_4_gene_id(
                &plus_gene
            ),
        2.0,
        "REF + ALT unspliced reads should produce two gene_plus intronic UMIs"
    );

    assert_eq!(
        exonic_cell
            .total_umis_4_gene_id(
                &minus_gene
            ),
        0.0,
        "plus-strand spliced reads must not be assigned to gene_minus"
    );

    assert_eq!(
        intronic_cell
            .total_umis_4_gene_id(
                &minus_gene
            ),
        0.0,
        "plus-strand intronic reads must not be assigned to gene_minus"
    );

    assert_eq!(
        exonic_cell.total_umis(),
        2,
        "expected exactly two exonic UMIs"
    );

    assert_eq!(
        intronic_cell.total_umis(),
        2,
        "expected exactly two intronic UMIs"
    );


    // --------------------------------------------------------
    // 8. Reload SNP matrices through Scdata
    // --------------------------------------------------------

    let ref_data =
        Scdata::read_matrix_market(
            &snp_ref,
            &snp_index,
        )
        .expect(
            "failed to reload SNP REF matrix"
        );

    let alt_data =
        Scdata::read_matrix_market(
            &snp_alt,
            &snp_index,
        )
        .expect(
            "failed to reload SNP ALT matrix"
        );


    // --------------------------------------------------------
    // Ground truth for SNP side channel
    //
    // Position:
    //
    //     chrTiny:135 A>C
    //
    // Reads:
    //
    //     spliced_REF      A
    //     spliced_ALT      C
    //     unspliced_REF    A
    //     unspliced_ALT    C
    //
    // Each read contains exactly ONE SNP/mismatch.
    //
    // Therefore:
    //
    //     REF = 2 UMIs
    //     ALT = 2 UMIs
    //
    // --------------------------------------------------------

    assert_eq!(
        ref_data.cell_ids().len(),
        1,
        "expected exactly one cell in SNP REF output"
    );

    assert_eq!(
        alt_data.cell_ids().len(),
        1,
        "expected exactly one cell in SNP ALT output"
    );

    assert!(
        ref_data
            .cell_ids()
            .contains(&cell_id),
        "gene and SNP REF outputs should contain the same cell"
    );

    assert!(
        alt_data
            .cell_ids()
            .contains(&cell_id),
        "gene and SNP ALT outputs should contain the same cell"
    );

    let ref_cell =
        ref_data
            .get(&cell_id)
            .expect(
                "cell missing from SNP REF Scdata"
            );

    let alt_cell =
        alt_data
            .get(&cell_id)
            .expect(
                "cell missing from SNP ALT Scdata"
            );

    assert_eq!(
        ref_cell
            .total_umis_4_gene_id(
                &tiny_snp
            ),
        2.0,
        "expected exactly two REF-supporting UMIs for chrTiny:135 A>C"
    );

    assert_eq!(
        alt_cell
            .total_umis_4_gene_id(
                &tiny_snp
            ),
        2.0,
        "expected exactly two ALT-supporting UMIs for chrTiny:135 A>C"
    );

    assert_eq!(
        ref_cell.total_umis(),
        2,
        "SNP REF matrix should contain exactly two UMIs"
    );

    assert_eq!(
        alt_cell.total_umis(),
        2,
        "SNP ALT matrix should contain exactly two UMIs"
    );
}