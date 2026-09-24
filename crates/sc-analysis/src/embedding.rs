use anyhow::{Context, Result};
use manifolds_rs::{UmapParams, umap};
use ndarray::Array2;

/// Standard 2-D UMAP used only for visualisation.
///
/// Clustering is deliberately performed in PCA space; this embedding never
/// participates in population calling. Keeping that separation makes the
/// plotted geometry diagnostic rather than part of the biological decision.
pub(crate) fn umap_2d(pca: &Array2<f32>, neighbors: usize, epochs: usize) -> Result<Array2<f32>> {
    let n = pca.nrows();
    if n == 0 {
        return Ok(Array2::zeros((0, 2)));
    }
    if n < 3 {
        return Ok(Array2::zeros((n, 2)));
    }

    // manifolds-rs accepts a row-major slice tuple, which avoids coupling its
    // ndarray version to Lumrik's ndarray version.
    let row_major = pca.as_standard_layout().to_owned();
    let values = row_major
        .as_slice()
        .context("PCA coordinates are not contiguous")?;
    let mut params = UmapParams::<f32>::new_default_2d(Some(0.5), Some(1.0));
    params.k = neighbors.max(2).min(n - 1);
    params.optim_params.n_epochs = epochs.max(1);
    let coords = umap((values, n, pca.ncols()), None, &params, 42, 0)
        .map_err(|e| anyhow::anyhow!("UMAP failed: {e}"))?;
    if coords.len() != 2 || coords[0].len() != n || coords[1].len() != n {
        anyhow::bail!("UMAP returned unexpected shape");
    }
    let mut out = Array2::<f32>::zeros((n, 2));
    for i in 0..n {
        out[(i, 0)] = coords[0][i];
        out[(i, 1)] = coords[1][i];
    }
    Ok(out)
}
