// test-gzipped-files.rs
use flate2::read::GzDecoder;
use mapping_info::MappingInfo;
use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};

use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

/// Minimal feature index for testing.
///
/// Feature IDs are the vector indices as u64:
/// Gene1 -> 0, Gene2 -> 1, Gene3 -> 2, Gene4 -> 3
struct TestFeatureIndex {
    names: Vec<String>,
    name_to_id: HashMap<String, u64>,
}

impl TestFeatureIndex {
    fn new(names: Vec<&str>) -> Self {
        let names: Vec<String> = names.into_iter().map(|s| s.to_string()).collect();
        let name_to_id = names
            .iter()
            .enumerate()
            .map(|(i, name)| (name.clone(), i as u64))
            .collect();

        Self { names, name_to_id }
    }

    fn feature_id(&self, name: &str) -> u64 {
        *self
            .name_to_id
            .get(name)
            .unwrap_or_else(|| panic!("unknown feature name: {name}"))
    }
}

impl FeatureIndex for TestFeatureIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        &self.names[feature_id as usize]
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let name = self.feature_name(feature_id);
        format!("{name}\t{name}\tGene Expression")
    }
    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.names.len() as u64).collect()
    }
    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }
}

fn test_out_dir(name: &str) -> PathBuf {
    PathBuf::from("tests").join("out").join(name)
}

fn read_gz_to_string(path: &PathBuf) -> String {
    let file =
        File::open(path).unwrap_or_else(|e| panic!("failed to open {}: {e}", path.display()));
    let mut decoder = GzDecoder::new(file);
    let mut text = String::new();
    decoder
        .read_to_string(&mut text)
        .unwrap_or_else(|e| panic!("failed to read gz {}: {e}", path.display()));
    text
}

#[test]
fn integer_10x_gz_files_match_expected_exactly() {
    let mut scdata = Scdata::new(1, MatrixValueType::Integer);
    let mut report = MappingInfo::new(None, 56.0, 10_000);

    let feature_index = TestFeatureIndex::new(vec!["Gene1", "Gene2", "Gene3", "Gene4"]);

    let gene1 = feature_index.feature_id("Gene1"); // 0
    let gene3 = feature_index.feature_id("Gene3"); // 2
    let gene4 = feature_index.feature_id("Gene4"); // 3

    // Cell 1: Gene1 x20, Gene4 x10
    for umi in 0..20 {
        assert!(
            scdata.try_insert(&1_u64, GeneUmiHash(gene1, umi as u64), 0.0, &mut report),
            "insert Gene1 umi {umi}"
        );
    }
    for umi in 20..30 {
        assert!(
            scdata.try_insert(&1_u64, GeneUmiHash(gene4, umi as u64), 0.0, &mut report),
            "insert Gene4 umi {umi}"
        );
    }

    // Cell 2: Gene3 x20, Gene1 x20
    for umi in 0..20 {
        assert!(
            scdata.try_insert(
                &13_452_355_u64,
                GeneUmiHash(gene3, umi as u64),
                0.0,
                &mut report
            ),
            "insert Gene3 umi {umi}"
        );
        assert!(
            scdata.try_insert(
                &13_452_355_u64,
                GeneUmiHash(gene1, umi as u64),
                0.0,
                &mut report
            ),
            "insert Gene1 umi {umi}"
        );
    }

    let out_dir = test_out_dir("integer_exact_file");
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir).unwrap();
    }

    scdata.finalize_for_export(0, &feature_index);
    scdata
        .write_sparse(&out_dir, &feature_index)
        .unwrap_or_else(|e| panic!("write_sparse failed: {e}"));

    let matrix_path = out_dir.join("matrix.mtx.gz");
    let features_path = out_dir.join("features.tsv.gz");
    let barcodes_path = out_dir.join("barcodes.tsv.gz");

    assert!(matrix_path.exists(), "matrix.mtx.gz not found");
    assert!(features_path.exists(), "features.tsv.gz not found");
    assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

    let actual_matrix = read_gz_to_string(&matrix_path);
    let actual_features = read_gz_to_string(&features_path);
    let actual_barcodes = read_gz_to_string(&barcodes_path);

    let expected_matrix = concat!(
        "%%MatrixMarket matrix coordinate integer general\n",
        "3 2 4\n",
        "1 1 20\n",
        "3 1 10\n",
        "1 2 20\n",
        "2 2 20\n",
    );

    let expected_features = concat!(
        "Gene1\tGene1\tGene Expression\n",
        "Gene3\tGene3\tGene Expression\n",
        "Gene4\tGene4\tGene Expression\n",
    );

    let expected_barcodes = concat!(
        "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n",
        "TAACACACCTATAAAAAAAAAAAAAAAAAAAA\n",
    );

    assert_eq!(
        actual_matrix, expected_matrix,
        "matrix.mtx.gz content mismatch"
    );
    assert_eq!(
        actual_features, expected_features,
        "features.tsv.gz content mismatch"
    );
    assert_eq!(
        actual_barcodes, expected_barcodes,
        "barcodes.tsv.gz content mismatch"
    );
}
