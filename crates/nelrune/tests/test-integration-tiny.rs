use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use bam_tide::index::GeneFeatureIndex;
use gtf_splice_index::{AnnotationBuilder, SpliceIndex};
use scdata::{Scdata, feature_index::FeatureIndex};

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

fn ensure_star_index(fasta: &Path, gtf: &Path, star_index: &Path) {
    if star_index.join("Genome").is_file() {
        return;
    }

    assert!(
        star_available(),
        "STAR is required for this integration test but was not found in PATH"
    );

    std::fs::create_dir_all(star_index).expect("failed to create STAR index directory");

    let status = Command::new("STAR")
        .args([
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            star_index.to_str().unwrap(),
            "--genomeFastaFiles",
            fasta.to_str().unwrap(),
            "--sjdbGTFfile",
            gtf.to_str().unwrap(),
            // Tiny artificial genome.
            "--genomeSAindexNbases",
            "2",
            // Our reads are 60 bp.
            "--sjdbOverhang",
            "59",
        ])
        .status()
        .expect("failed to start STAR genomeGenerate");

    assert!(status.success(), "STAR genomeGenerate failed");

    assert!(
        star_index.join("Genome").is_file(),
        "STAR finished successfully but no Genome file was created"
    );
}

fn ensure_splice_index(gtf: &Path, index_path: &Path) {
    if index_path.is_file() {
        return;
    }

    let index = AnnotationBuilder::new(1_000)
        .build_from_path(gtf)
        .expect("failed to build tiny splice index");

    index
        .save(index_path)
        .expect("failed to write tiny splice index");

    assert!(index_path.is_file(), "splice index was not created");
}

fn assert_matrix_files(path: &Path) {
    for name in ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"] {
        assert!(
            path.join(name).is_file(),
            "missing expected matrix output: {}",
            path.join(name).display()
        );
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            other => panic!("unexpected base in tiny FASTQ fixture: {}", other as char),
        })
        .collect()
}

fn write_two_cell_diagonal_fixture(source_r1: &Path, source_r2: &Path, r1: &Path, r2: &Path) {
    let r1_text = fs::read_to_string(source_r1)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", source_r1.display()));
    let r2_text = fs::read_to_string(source_r2)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", source_r2.display()));

    let r1_lines: Vec<&str> = r1_text.lines().collect();
    let r2_lines: Vec<&str> = r2_text.lines().collect();

    assert_eq!(
        r1_lines.len(),
        8,
        "tiny R1 fixture must contain exactly two reads"
    );
    assert_eq!(
        r2_lines.len(),
        8,
        "tiny R2 fixture must contain exactly two reads"
    );

    let mut out_r1 = String::new();
    let mut out_r2 = String::new();

    // Cell 1 is the existing fixture: its two reads map to gene_plus, one
    // exonic and one intronic. Cell 2 gets the same molecules on the opposite
    // strand, so they map to gene_minus. This creates the deliberate 2 x 2
    // diagonal matrix in both outputs:
    //
    //              gene_plus  gene_minus
    //     cell 1       1          0
    //     cell 2       0          1
    //
    // The absent combinations are the regression target: a missing sparse
    // entry is a zero, not a missing cell or feature.
    for read in 0..2 {
        let i = read * 4;
        out_r1.push_str(&format!(
            "{}\n{}\n{}\n{}\n",
            r1_lines[i],
            r1_lines[i + 1],
            r1_lines[i + 2],
            r1_lines[i + 3]
        ));
        out_r2.push_str(&format!(
            "{}\n{}\n{}\n{}\n",
            r2_lines[i],
            r2_lines[i + 1],
            r2_lines[i + 2],
            r2_lines[i + 3]
        ));
    }

    for read in 0..2 {
        let i = read * 4;
        let source_cell_read = r1_lines[i + 1];
        assert!(
            source_cell_read.len() >= 4,
            "tiny R1 read is shorter than CELL:4"
        );
        let second_cell_read = format!("TGCA{}", &source_cell_read[4..]);
        let second_name = format!("{}-cell2", r1_lines[i]);

        out_r1.push_str(&format!(
            "{second_name}\n{second_cell_read}\n{}\n{}\n",
            r1_lines[i + 2],
            r1_lines[i + 3]
        ));

        let second_r2_name = format!("{}-cell2", r2_lines[i]);
        let second_r2_seq = reverse_complement(r2_lines[i + 1]);
        let second_r2_qual: String = r2_lines[i + 3].chars().rev().collect();
        out_r2.push_str(&format!(
            "{second_r2_name}\n{second_r2_seq}\n{}\n{second_r2_qual}\n",
            r2_lines[i + 2]
        ));
    }

    fs::write(r1, out_r1).unwrap_or_else(|err| panic!("failed to write {}: {err}", r1.display()));
    fs::write(r2, out_r2).unwrap_or_else(|err| panic!("failed to write {}: {err}", r2.display()));
}

#[test]
fn integration_tiny_star() {
    // --------------------------------------------------------
    // 1. Test inputs
    // --------------------------------------------------------

    let fasta = test_path("tiny.fa");

    let gtf = test_path("tiny.gtf");

    let source_r1 = test_path("tiny_R1.fastq");

    let source_r2 = test_path("tiny_R2.fastq");

    for path in [&fasta, &gtf, &source_r1, &source_r2] {
        require_file(path);
    }

    // --------------------------------------------------------
    // 2. STAR index
    // --------------------------------------------------------

    let star_index = test_path("star_index");

    ensure_star_index(&fasta, &gtf, &star_index);

    // --------------------------------------------------------
    // 3. gtf-splice-index
    // --------------------------------------------------------

    let splice_index_path = test_path("tiny.gtf.dat");

    ensure_splice_index(&gtf, &splice_index_path);

    // Make sure it can actually be loaded again.
    let splice_index =
        SpliceIndex::load(&splice_index_path).expect("created splice index cannot be loaded");

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

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test-integration-tiny");

    let _ = std::fs::remove_dir_all(&out);

    std::fs::create_dir_all(&out).expect("failed to create test output directory");

    let r1 = out.join("two_cell_R1.fastq");
    let r2 = out.join("two_cell_R2.fastq");
    write_two_cell_diagonal_fixture(&source_r1, &source_r2, &r1, &r2);

    let output = Command::new(env!("CARGO_BIN_EXE_nelrune"))
        .args([
            "--r1",
            r1.to_str().unwrap(),
            "--r2",
            r2.to_str().unwrap(),
            "--primer-structure",
            "TYPE:GEX+CELL:4+UMI:4",
            "--mapper",
            "star",
            "--mapper-index",
            star_index.to_str().unwrap(),
            "--mapper-threads",
            "2",
            "--index",
            splice_index_path.to_str().unwrap(),
            "--require-strand",
            "--min-mapq",
            "0",
            // Both synthetic cells have two molecules and must survive filtering.
            "--min-umi-count",
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
        .expect("failed to start nelrune");

    if !output.status.success() {
        panic!(
            "\nNelrune integration test failed.\n\
             status: {}\n\
             stdout:\n{}\n\
             stderr:\n{}\n",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr),
        );
    }

    // --------------------------------------------------------
    // 5. Full output stack exists
    // --------------------------------------------------------

    let exonic = out.join("filtered/exonic");
    let raw_exonic = out.join("raw/exonic");

    let intronic = out.join("filtered/intronic");
    let raw_intronic = out.join("raw/intronic");

    assert_matrix_files(&exonic);
    assert_matrix_files(&intronic);

    assert!(out.join("nelrune.log").is_file(), "Nelrune log is missing");

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

    let gene_index = GeneFeatureIndex::new(&splice_index);

    let plus_gene = gene_index
        .feature_id("gene_plus")
        .expect("gene_plus missing from feature index");

    let minus_gene = gene_index
        .feature_id("gene_minus")
        .expect("gene_minus missing from feature index");

    let exonic_data = Scdata::read_matrix_market(&exonic, &gene_index)
        .expect("failed to reload exonic Nelrune matrix");

    let intronic_data = Scdata::read_matrix_market(&intronic, &gene_index)
        .expect("failed to reload intronic Nelrune matrix");

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
        2,
        "expected both cells in exonic output"
    );

    assert_eq!(
        intronic_data.cell_ids().len(),
        2,
        "expected both cells in intronic output"
    );

    assert_eq!(
        exonic_data.cell_ids(),
        intronic_data.cell_ids(),
        "exonic and intronic outputs must retain the same two-cell axis"
    );

    let mut exonic_patterns = Vec::new();
    let mut intronic_patterns = Vec::new();

    for cell_id in exonic_data.cell_ids() {
        let exonic_cell = exonic_data
            .get(&cell_id)
            .expect("cell missing from exonic Scdata");
        let intronic_cell = intronic_data
            .get(&cell_id)
            .expect("cell missing from intronic Scdata");

        exonic_patterns.push((
            exonic_cell.total_umis_4_gene_id(&plus_gene),
            exonic_cell.total_umis_4_gene_id(&minus_gene),
        ));
        intronic_patterns.push((
            intronic_cell.total_umis_4_gene_id(&plus_gene),
            intronic_cell.total_umis_4_gene_id(&minus_gene),
        ));

        assert_eq!(exonic_cell.total_umis(), 1);
        assert_eq!(intronic_cell.total_umis(), 1);
    }

    exonic_patterns.sort_by(|a, b| a.partial_cmp(b).unwrap());
    intronic_patterns.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let expected = vec![(0.0, 1.0), (1.0, 0.0)];
    assert_eq!(
        exonic_patterns, expected,
        "exonic output must be a 2-cell x 2-gene diagonal matrix with one implicit zero per cell"
    );
    assert_eq!(
        intronic_patterns, expected,
        "intronic output must be a 2-cell x 2-gene diagonal matrix with one implicit zero per cell"
    );
}
