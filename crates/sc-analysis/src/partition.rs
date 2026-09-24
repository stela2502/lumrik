use ndarray::Array2;

/// Deterministically partition PCA space into small axis-aligned micro-patches.
///
/// These labels are only seeds for the existing mean-expression merge.  They
/// deliberately carry no biological meaning by themselves.
pub(crate) fn spatial_patches(pca: &Array2<f32>, explained: &[f64]) -> Vec<usize> {
    let n = pca.nrows();
    if n == 0 {
        return Vec::new();
    }

    let max_patch_cells = 50usize.max((n + 99) / 100); // ceil(1% of retained cells)
    let dims = informative_dimensions(pca.ncols(), explained);
    let mut leaves = Vec::<Vec<usize>>::new();
    split_patch(
        pca,
        (0..n).collect(),
        &dims,
        max_patch_cells,
        &mut leaves,
    );

    let mut labels = vec![0usize; n];
    for (label, cells) in leaves.iter().enumerate() {
        for &cell in cells {
            labels[cell] = label;
        }
    }
    labels
}

/// Use only the leading PCs needed to explain 90% of the variance captured by
/// the computed PCA, with at least two dimensions when they exist.
fn informative_dimensions(ncols: usize, explained: &[f64]) -> Vec<usize> {
    if ncols == 0 {
        return Vec::new();
    }
    let available = ncols.min(explained.len());
    if available == 0 {
        return (0..ncols).collect();
    }
    let total = explained[..available].iter().copied().sum::<f64>();
    if total <= 0.0 {
        return (0..available).collect();
    }
    let target = total * 0.90;
    let mut cumulative = 0.0;
    let mut keep = 0usize;
    for &variance in &explained[..available] {
        cumulative += variance;
        keep += 1;
        if cumulative >= target {
            break;
        }
    }
    keep = keep.max(2.min(available));
    (0..keep).collect()
}

fn split_patch(
    pca: &Array2<f32>,
    mut cells: Vec<usize>,
    dims: &[usize],
    max_patch_cells: usize,
    leaves: &mut Vec<Vec<usize>>,
) {
    if cells.len() <= max_patch_cells || cells.len() < 2 || dims.is_empty() {
        leaves.push(cells);
        return;
    }

    let axis = highest_variance_axis(pca, &cells, dims);
    cells.sort_unstable_by(|&a, &b| {
        pca[(a, axis)]
            .total_cmp(&pca[(b, axis)])
            .then(a.cmp(&b))
    });
    let right = cells.split_off(cells.len() / 2);
    split_patch(pca, cells, dims, max_patch_cells, leaves);
    split_patch(pca, right, dims, max_patch_cells, leaves);
}

fn highest_variance_axis(pca: &Array2<f32>, cells: &[usize], dims: &[usize]) -> usize {
    let mut best_axis = dims[0];
    let mut best_variance = -1.0f64;
    for &axis in dims {
        let mean = cells
            .iter()
            .map(|&cell| pca[(cell, axis)] as f64)
            .sum::<f64>()
            / cells.len() as f64;
        let variance = cells
            .iter()
            .map(|&cell| {
                let delta = pca[(cell, axis)] as f64 - mean;
                delta * delta
            })
            .sum::<f64>()
            / cells.len() as f64;
        if variance > best_variance {
            best_variance = variance;
            best_axis = axis;
        }
    }
    best_axis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn patch_size_is_fifty_for_small_inputs() {
        let mut pca = Array2::<f32>::zeros((1076, 3));
        for i in 0..pca.nrows() {
            pca[(i, 0)] = i as f32;
            pca[(i, 1)] = (i % 17) as f32;
        }
        let labels = spatial_patches(&pca, &[0.6, 0.3, 0.1]);
        let k = labels.iter().copied().max().unwrap() + 1;
        let mut counts = vec![0usize; k];
        for label in labels {
            counts[label] += 1;
        }
        assert!(counts.iter().all(|&n| n <= 50));
    }

    #[test]
    fn patch_size_is_one_percent_for_large_inputs() {
        let n = 10_000;
        let mut pca = Array2::<f32>::zeros((n, 2));
        for i in 0..n {
            pca[(i, 0)] = i as f32;
        }
        let labels = spatial_patches(&pca, &[0.9, 0.1]);
        let k = labels.iter().copied().max().unwrap() + 1;
        let mut counts = vec![0usize; k];
        for label in labels {
            counts[label] += 1;
        }
        assert!(counts.iter().all(|&n| n <= 100));
    }
}
