use crate::data::{CellAnnotations, SingleCellData};
use crate::qc::QcMetrics;
use anyhow::Result;
use sprs::{CsMat, TriMat};

pub(crate) fn normalize_surviving(
    raw: &SingleCellData,
    qc: &QcMetrics,
    scale: f64,
) -> Result<SingleCellData> {
    let kept = qc
        .keep_gene
        .iter()
        .enumerate()
        .filter_map(|(i, k)| k.then_some(i))
        .collect::<Vec<_>>();
    let mut tri = TriMat::<f32>::with_capacity((kept.len(), raw.n_cells()), raw.matrix.nnz());
    let mut features = Vec::with_capacity(kept.len());
    for (ng, &g) in kept.iter().enumerate() {
        features.push(raw.features[g].clone());
        if let Some(row) = raw.matrix.outer_view(g) {
            for (c, v) in row.iter() {
                if *v > 0.0 && qc.surviving[c] > 0.0 {
                    tri.add_triplet(
                        ng,
                        c,
                        (1.0 + scale * (*v as f64) / qc.surviving[c]).ln() as f32,
                    );
                }
            }
        }
    }
    let matrix: CsMat<f32> = tri.to_csr();
    let annotations =
        CellAnnotations::from_columns(vec![("cell", raw.cells.clone())], raw.n_cells())?;
    SingleCellData::new(matrix, features, raw.cells.clone(), annotations)
}
