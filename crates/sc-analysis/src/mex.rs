use crate::data::{CellAnnotations, SingleCellData};
use anyhow::{Context, Result};
use scdata::{FeatureIndex, load_mtx_feature_matrix, read_mtx_barcodes};
use std::path::Path;

pub(crate) fn load_mex(dir: &Path) -> Result<SingleCellData> {
    // scdata owns the external MEX contract. Its importer accepts standard
    // coordinate MatrixMarket input independent of whether entries arrive in
    // row-major, column-major, or otherwise valid coordinate order.
    let (counts, feature_index, _) = load_mtx_feature_matrix(dir, "Gene Expression", 1)
        .with_context(|| format!("loading Gene Expression MEX from {}", dir.display()))?;
    let barcodes = read_mtx_barcodes(dir)?;
    let cells = barcodes
        .into_iter()
        .map(|(barcode, _)| barcode)
        .collect::<Vec<_>>();
    let matrix = counts.as_sprs().map_err(anyhow::Error::msg)?;
    let features = feature_index
        .ordered_feature_ids()
        .into_iter()
        .map(|feature_id| feature_index.feature_name(feature_id).to_string())
        .collect::<Vec<_>>();
    let annotations = CellAnnotations::from_columns(vec![("cell", cells.clone())], cells.len())?;

    SingleCellData::new(matrix, features, cells, annotations)
}
