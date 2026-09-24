use flate2::{Compression, write::GzEncoder};
use scdata::{FeatureIndex, load_mtx_feature_matrix, read_mtx_barcodes};
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn write_gz(path: &Path, text: &str) {
    let mut out = GzEncoder::new(File::create(path).unwrap(), Compression::default());
    out.write_all(text.as_bytes()).unwrap();
    out.finish().unwrap();
}

fn write_mex(dir: &Path, entries: &str) {
    write_gz(
        &dir.join("features.tsv.gz"),
        "g1\tGene1\tGene Expression\ng2\tGene2\tGene Expression\n",
    );
    write_gz(&dir.join("barcodes.tsv.gz"), "AAAA-1\nCCCC-1\n");
    write_gz(
        &dir.join("matrix.mtx.gz"),
        &format!(
            "%%MatrixMarket matrix coordinate integer general\n% test\n2 2 4\n{entries}"
        ),
    );
}

fn load_dense(dir: &Path) -> Vec<Vec<f32>> {
    let (data, index, _) = load_mtx_feature_matrix(dir, "Gene Expression", 1).unwrap();
    let matrix = data.as_sprs().unwrap();
    assert_eq!(index.ordered_feature_ids(), vec![0, 1]);
    matrix.to_dense().outer_iter().map(|row| row.to_vec()).collect()
}

#[test]
fn mtx_import_accepts_row_and_column_major_coordinates() {
    let row_major = tempfile::tempdir().unwrap();
    let column_major = tempfile::tempdir().unwrap();

    write_mex(row_major.path(), "1 1 1\n1 2 2\n2 1 3\n2 2 4\n");
    write_mex(column_major.path(), "1 1 1\n2 1 3\n1 2 2\n2 2 4\n");

    let expected = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
    assert_eq!(load_dense(row_major.path()), expected);
    assert_eq!(load_dense(column_major.path()), expected);
}

#[test]
fn mtx_import_accepts_numeric_vendor_cell_ids() {
    let dir = tempfile::tempdir().unwrap();
    write_gz(
        &dir.path().join("features.tsv.gz"),
        "g1\tGene1\tGene Expression\ng2\tGene2\tGene Expression\n",
    );
    write_gz(&dir.path().join("barcodes.tsv.gz"), "1764\n4422\n");
    write_gz(
        &dir.path().join("matrix.mtx.gz"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 7\n2 2 9\n",
    );

    let barcodes = read_mtx_barcodes(dir.path()).unwrap();
    assert_eq!(barcodes, vec![("1764".to_string(), 1764), ("4422".to_string(), 4422)]);

    let (data, _, _) = load_mtx_feature_matrix(dir.path(), "Gene Expression", 1).unwrap();
    assert_eq!(
        data.as_sprs().unwrap().to_dense().outer_iter().map(|row| row.to_vec()).collect::<Vec<_>>(),
        vec![vec![7.0, 0.0], vec![0.0, 9.0]],
    );
}

#[test]
fn mtx_import_preserves_arbitrary_cell_labels() {
    let dir = tempfile::tempdir().unwrap();
    write_gz(
        &dir.path().join("features.tsv.gz"),
        "g1\tGene1\tGene Expression\ng2\tGene2\tGene Expression\n",
    );
    write_gz(&dir.path().join("barcodes.tsv.gz"), "cell-A\nHorst\n");
    write_gz(
        &dir.path().join("matrix.mtx.gz"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 3\n2 2 5\n",
    );

    let first = read_mtx_barcodes(dir.path()).unwrap();
    let second = read_mtx_barcodes(dir.path()).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.iter().map(|(label, _)| label.as_str()).collect::<Vec<_>>(), vec!["cell-A", "Horst"]);

    let (data, _, _) = load_mtx_feature_matrix(dir.path(), "Gene Expression", 1).unwrap();
    assert_eq!(
        data.as_sprs().unwrap().to_dense().outer_iter().map(|row| row.to_vec()).collect::<Vec<_>>(),
        vec![vec![3.0, 0.0], vec![0.0, 5.0]],
    );
}
