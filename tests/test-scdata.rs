use mapping_info::MappingInfo;
use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};

use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

/// Minimal feature index for integration tests.
///
/// Feature IDs are just the index in `names` as `u64`.
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

#[test]
fn singlecelldata_to_sparse_integer_roundtrip() {
    let mut celldata = Scdata::new(1, MatrixValueType::Integer);
    let mut report = MappingInfo::new(None, 56.0, 10_000);

    let feature_index = TestFeatureIndex::new(vec!["Gene1", "Gene2", "Gene3", "Gene4"]);

    let gene1 = feature_index.feature_id("Gene1");
    let gene3 = feature_index.feature_id("Gene3");
    let gene4 = feature_index.feature_id("Gene4");

    // Cell 1: Gene1 x20, Gene4 x10
    assert!(
        celldata.try_insert_value(&1_u64, GeneUmiHash(gene1, 0), 20.0, &mut report),
        "insert Gene1 umi 0"
    );

    assert!(
        !celldata.try_insert_value(&1_u64, GeneUmiHash(gene1, 0), 1.0, &mut report),
        "duplicate Gene1 umi 0 should fail"
    );

    assert!(
        celldata.try_insert_value(&1_u64, GeneUmiHash(gene4, 0), 10.0, &mut report),
        "insert Gene4 umi 0"
    );

    // Cell 2: Gene3 x30, Gene1 x20
    assert!(
        celldata.try_insert_value(&13_452_355_u64, GeneUmiHash(gene3, 0), 30.0, &mut report),
        "insert Gene3 umi 0"
    );
    assert!(
        celldata.try_insert_value(&13_452_355_u64, GeneUmiHash(gene1, 0), 20.0, &mut report),
        "insert Gene1 umi 0"
    );

    assert!(
        !celldata.try_insert(&13_452_355_u64, GeneUmiHash(gene1, 0), 1.0, &mut report),
        "duplicate Gene1 umi 0 should fail"
    );

    let out_dir = test_out_dir("integer");
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir).unwrap();
    }
    celldata.finalize_for_export(0, &feature_index);
    let write_result = celldata.write_sparse(&out_dir, &feature_index);
    assert!(
        write_result.is_ok(),
        "write_sparse failed: {write_result:?}"
    );

    let matrix_path = out_dir.join("matrix.mtx.gz");
    let features_path = out_dir.join("features.tsv.gz");
    let barcodes_path = out_dir.join("barcodes.tsv.gz");

    assert!(matrix_path.exists(), "matrix.mtx.gz not found");
    assert!(features_path.exists(), "features.tsv.gz not found");
    assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

    // Assumption:
    // read_matrix_market reconstructs feature IDs in the same order as written
    // from features.tsv.gz, so Gene1..Gene4 keep IDs 0..3.
    let scdata2 = Scdata::read_matrix_market(&out_dir, &feature_index).unwrap();

    // Cell 1: Gene1 x20, Gene4 x10
    assert_eq!(
        scdata2.get(&1_u64).unwrap().total_umis_4_gene_id(&gene1),
        20.0
    );
    assert_eq!(
        scdata2.get(&1_u64).unwrap().total_umis_4_gene_id(&gene3),
        0.0
    );
    assert_eq!(
        scdata2.get(&1_u64).unwrap().total_umis_4_gene_id(&gene4),
        10.0
    );

    // Cell 2: Gene3 x30, Gene1 x20
    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .total_umis_4_gene_id(&gene3),
        30.0
    );
    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .total_umis_4_gene_id(&gene1),
        20.0
    );
    assert_eq!(
        scdata2
            .get(&13_452_355_u64)
            .unwrap()
            .total_umis_4_gene_id(&gene4),
        0.0
    );
}

#[test]
fn singlecelldata_to_sparse_real_can_be_written_and_read() {
    let mut celldata = Scdata::new(1, MatrixValueType::Real);
    let mut report = MappingInfo::new(None, 56.0, 10_000);

    let feature_index = TestFeatureIndex::new(vec!["Gene1", "Gene2", "Gene3", "Gene4"]);

    let gene1 = feature_index.feature_id("Gene1");
    let gene3 = feature_index.feature_id("Gene3");
    let gene4 = feature_index.feature_id("Gene4");

    // Cell 1: Gene1 total 20.0, Gene4 total 10.0
    assert!(celldata.try_insert_value(&1_u64, GeneUmiHash(gene1, 0), 20.0, &mut report));
    assert!(celldata.try_insert_value(&1_u64, GeneUmiHash(gene4, 0), 10.0, &mut report));

    // Cell 2: Gene3 total 20.0, Gene1 total 20.0
    assert!(celldata.try_insert_value(&13_452_355_u64, GeneUmiHash(gene3, 0), 20.0, &mut report));
    assert!(celldata.try_insert_value(&13_452_355_u64, GeneUmiHash(gene1, 0), 20.0, &mut report));

    // Duplicate should still be rejected.
    assert!(
        !celldata.try_insert_value(&13_452_355_u64, GeneUmiHash(gene1, 0), 2.0, &mut report),
        "duplicate Gene1 umi 0 should fail"
    );

    let out_dir = test_out_dir("real");
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir).unwrap();
    }

    celldata.finalize_for_export(0, &feature_index);

    let write_result = celldata.write_sparse(&out_dir, &feature_index);
    assert!(
        write_result.is_ok(),
        "write_sparse failed: {write_result:?}"
    );

    let matrix_path = out_dir.join("matrix.mtx.gz");
    let features_path = out_dir.join("features.tsv.gz");
    let barcodes_path = out_dir.join("barcodes.tsv.gz");

    assert!(matrix_path.exists(), "matrix.mtx.gz not found");
    assert!(features_path.exists(), "features.tsv.gz not found");
    assert!(barcodes_path.exists(), "barcodes.tsv.gz not found");

    // Reading should at least succeed.
    let scdata2 = Scdata::read_matrix_market(&out_dir, &feature_index);
    assert!(
        scdata2.is_ok(),
        "read_matrix_market failed: {}",
        scdata2.unwrap()
    );
}
