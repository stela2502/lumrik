use crate::data::SingleCellData;
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
    let annotations = raw.annotations.clone();
    SingleCellData::new(matrix, features, raw.cells.clone(), annotations)
}


/// Fixed-depth stochastic rounding of UMI counts, followed by an exact-sum
/// correction. This follows FastWilcoxTest's NormalizeCells design, rather
/// than sampling individual molecules without replacement.
/// The RNG stream is independently keyed by cell index for reproducibility.
pub(crate) fn resample_surviving(
    raw: &SingleCellData,
    qc: &QcMetrics,
    target: usize,
    seed: u64,
) -> Result<SingleCellData> {
    anyhow::ensure!(target > 0, "resampling target must be positive");
    let kept: Vec<usize> = qc.keep_gene.iter().enumerate()
        .filter_map(|(g, &yes)| yes.then_some(g)).collect();
    let mut by_cell = vec![Vec::<(usize, u32)>::new(); raw.n_cells()];
    for (new_gene, &old_gene) in kept.iter().enumerate() {
        if let Some(row) = raw.matrix.outer_view(old_gene) {
            for (cell, &value) in row.iter() {
                if value > 0.0 {
                    anyhow::ensure!(value.is_finite() && value.fract() == 0.0 && value <= u32::MAX as f32,
                        "UMI resampling requires integer raw counts (gene {old_gene}, cell {cell}: {value})");
                    by_cell[cell].push((new_gene, value as u32));
                }
            }
        }
    }
    let mut tri = TriMat::<f32>::with_capacity((kept.len(), raw.n_cells()), raw.matrix.nnz());
    for (cell, genes) in by_cell.iter().enumerate() {
        let depth: u64 = genes.iter().map(|(_, v)| *v as u64).sum();
        anyhow::ensure!(depth >= target as u64, "cell {cell} has {depth} surviving UMIs, below required resampling depth {target}");
        if depth == target as u64 {
            for &(g, v) in genes { tri.add_triplet(g, cell, v as f32); }
            continue;
        }
        let mut state = seed ^ (cell as u64).wrapping_mul(0x9e3779b97f4a7c15);
        let mut rounded = Vec::with_capacity(genes.len());
        let mut sum = 0_i64;
        for &(g, count) in genes {
            let expected = count as f64 * target as f64 / depth as f64;
            let base = expected.floor() as i64;
            let frac = expected - base as f64;
            let draw = uniform01(&mut state);
            let value = base + i64::from(draw < frac);
            sum += value;
            // Positive priority: a rounded-up count is easier to decrement
            // when draw was close to frac; a rounded-down count is easier to
            // increment when draw was close to frac.
            rounded.push((g, value, frac - draw));
        }
        let delta = target as i64 - sum;
        if delta > 0 {
            rounded.sort_unstable_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)));
            for item in rounded.iter_mut().take(delta as usize) { item.1 += 1; }
        } else if delta < 0 {
            rounded.sort_unstable_by(|a, b| a.2.total_cmp(&b.2).then(a.0.cmp(&b.0)));
            let mut need = (-delta) as usize;
            for item in &mut rounded {
                if need == 0 { break; }
                if item.1 > 0 { item.1 -= 1; need -= 1; }
            }
            anyhow::ensure!(need == 0, "resampling correction failed for cell {cell}");
        }
        debug_assert_eq!(rounded.iter().map(|x| x.1).sum::<i64>(), target as i64);
        for (g, value, _) in rounded {
            if value > 0 { tri.add_triplet(g, cell, value as f32); }
        }
    }
    SingleCellData::new(tri.to_csr(), kept.iter().map(|&g| raw.features[g].clone()).collect(),
        raw.cells.clone(), raw.annotations.clone())
}

fn uniform01(state: &mut u64) -> f64 {
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut x = *state;
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    ((x ^ (x >> 31)) >> 11) as f64 * (1.0 / (1_u64 << 53) as f64)
}
