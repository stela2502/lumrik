use anyhow::{Context, Result};
use linfa::prelude::*;
use linfa_reduction::Pca;
use ndarray::Array2;

pub(crate) fn pca(x: &Array2<f32>, components: usize) -> Result<(Array2<f32>, Vec<f64>)> {
    let max_k = x.nrows().saturating_sub(1).min(x.ncols().saturating_sub(1));
    let k = components.max(1).min(max_k.max(1));
    let dataset = DatasetBase::from(x.mapv(|v| v as f64));
    let model = Pca::params(k).fit(&dataset).context("fitting PCA")?;
    let coords = model.transform(dataset).records.mapv(|v| v as f32);
    let total = x.ncols().max(1) as f64; // input genes are z-scored, so total variance ~= number of genes
    let mut variance = Vec::with_capacity(coords.ncols());
    for j in 0..coords.ncols() {
        let m = (0..coords.nrows())
            .map(|i| coords[(i, j)] as f64)
            .sum::<f64>()
            / coords.nrows().max(1) as f64;
        let v = (0..coords.nrows())
            .map(|i| {
                let z = coords[(i, j)] as f64 - m;
                z * z
            })
            .sum::<f64>()
            / coords.nrows().max(1) as f64;
        variance.push(v / total);
    }
    Ok((coords, variance))
}
