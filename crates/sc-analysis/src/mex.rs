use crate::data::{CellAnnotations, SingleCellData};
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use scdata::MexFeatureIndex;
use sprs::TriMat;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub(crate) fn load_mex(dir: &Path) -> Result<SingleCellData> {
    // MEX cell columns are positional. Never decode/hash external barcodes into
    // scdata's molecule-level IDs: independent libraries reuse 10x sequences.
    let feature_index = MexFeatureIndex::from_dir(dir, "Gene Expression")?;
    let barcode_path = ["barcodes.tsv.gz", "barcodes.tsv"]
        .iter().map(|name| dir.join(name)).find(|path| path.is_file())
        .with_context(|| format!("missing barcodes.tsv[.gz] in {}", dir.display()))?;
    let file = File::open(&barcode_path)?;
    let reader: Box<dyn BufRead> = if barcode_path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let cells: Vec<String> = reader.lines().collect::<std::io::Result<Vec<_>>>()?
        .into_iter().map(|line| line.trim().to_owned()).collect();
    if cells.iter().any(String::is_empty) { bail!("empty cell label in {}", barcode_path.display()); }
    let matrix_path = ["matrix.mtx.gz", "matrix.mtx"]
        .iter().map(|name| dir.join(name)).find(|path| path.is_file())
        .with_context(|| format!("missing matrix.mtx[.gz] in {}", dir.display()))?;
    let source = SingleCellData::read_matrix_market(&matrix_path)
        .with_context(|| format!("reading {}", matrix_path.display()))?;
    if source.cols() != cells.len() {
        bail!("MEX has {} matrix columns but {} cell labels", source.cols(), cells.len());
    }
    let selected = feature_index.features();
    if selected.iter().any(|feature| feature.source_row >= source.rows()) {
        bail!("feature row exceeds matrix dimensions in {}", dir.display());
    }
    let mut triplets = TriMat::<f32>::with_capacity((selected.len(), cells.len()), source.nnz());
    for (target_row, feature) in selected.iter().enumerate() {
        if let Some(row) = source.outer_view(feature.source_row) {
            for (col, value) in row.iter() {
                triplets.add_triplet(target_row, col, *value);
            }
        }
    }
    let features = selected.iter().map(|f| f.name.clone()).collect();
    let annotations = CellAnnotations::from_columns(vec![("cell", cells.clone())], cells.len())?;
    SingleCellData::new(triplets.to_csr(), features, cells, annotations)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn preserves_reused_barcodes_across_samples_and_unsorted_coordinates() {
        let dir = std::env::temp_dir().join(format!("sc_analysis_mex_{}_{}", std::process::id(), line!()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("features.tsv"), "A\tA\tGene Expression\nB\tB\tGene Expression\n").unwrap();
        std::fs::write(dir.join("barcodes.tsv"), "ACGT-1_BP08\nACGT-1_DBA1_3\n").unwrap();
        let mut f = File::create(dir.join("matrix.mtx")).unwrap();
        write!(f, "%%MatrixMarket matrix coordinate integer general\n2 2 2\n2 2 7\n1 1 3\n").unwrap();
        let loaded = load_mex(&dir).unwrap();
        assert_eq!(loaded.cells, vec!["ACGT-1_BP08", "ACGT-1_DBA1_3"]);
        assert_eq!(loaded.matrix.get(0, 0), Some(&3.0));
        assert_eq!(loaded.matrix.get(1, 1), Some(&7.0));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
