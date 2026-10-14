use sprs::CsMat;

#[derive(Debug, Clone)]
pub struct MergeStep {
    pub step: usize,
    pub from: usize,
    pub into: usize,
    pub correlation: f64,
    pub ward_delta_sse: Option<f64>,
    pub from_cells: usize,
    pub into_cells: usize,
    pub merged_cells: usize,
}

pub(crate) fn merge_correlated(
    x: &CsMat<f32>,
    selected_genes: &[usize],
    mut labels: Vec<usize>,
    cut: f64,
    target: Option<usize>,
    max_cluster_fraction: f64,
) -> anyhow::Result<(Vec<usize>, Vec<MergeStep>, f64)> {
    use rayon::prelude::*;

    anyhow::ensure!(max_cluster_fraction.is_finite() && (0.0..=1.0).contains(&max_cluster_fraction), "max_cluster_fraction must be in [0,1]");
    let max_cells = if max_cluster_fraction == 0.0 { usize::MAX } else { ((labels.len() as f64 * max_cluster_fraction).floor() as usize).max(1) };
    let ids = unique(&labels);
    let n_clusters = ids.len();
    anyhow::ensure!(target.is_none_or(|n| n > 0 && n <= n_clusters), "--target-clusters must be between 1 and {n_clusters}");
    let id_to_slot = ids.iter().enumerate().map(|(i, &id)| (id, i))
        .collect::<std::collections::HashMap<_, _>>();
    let mut cell_slots = labels.iter().map(|id| id_to_slot[id]).collect::<Vec<_>>();
    let mut sums = vec![vec![0.0f32; selected_genes.len()]; n_clusters];
    let mut counts = vec![0usize; n_clusters];
    for &slot in &cell_slots { counts[slot] += 1; }
    // The expression matrix is gene x cell (CSR). Visit it only once.
    for (slot, &gene) in selected_genes.iter().enumerate() {
        if let Some(row) = x.outer_view(gene) {
            for (cell, value) in row.iter() {
                sums[cell_slots[cell]][slot] += *value;
            }
        }
    }
    let mut means = sums.iter().zip(&counts).map(|(sum, &count)| {
        let denom = count.max(1) as f32;
        sum.iter().map(|v| *v / denom).collect::<Vec<_>>()
    }).collect::<Vec<_>>();
    let mut active = vec![true; n_clusters];
    // Cache pairwise correlations. Independent initial rows use Rayon.
    let mut correlations = (0..n_clusters).into_par_iter().map(|i| {
        (0..n_clusters).map(|j| {
            if i == j { -2.0 } else { pearson_f32_noalloc(&means[i], &means[j]) }
        }).collect::<Vec<f64>>()
    }).collect::<Vec<_>>();
    let mut history = Vec::new();
    let mut remaining_count = n_clusters;
    let highest_remaining = loop {
        // Globally highest rho, deterministic tie-breaking by cluster ID.
        let mut best = (-2.0f64, 0usize, 0usize);
        for from in 0..n_clusters {
            if !active[from] { continue; }
            for into in (from + 1)..n_clusters {
                if active[into] && counts[from] + counts[into] <= max_cells && correlations[from][into] > best.0 {
                    best = (correlations[from][into], from, into);
                }
            }
        }
        if remaining_count <= 1 || target.is_some_and(|n| remaining_count <= n)
            || (target.is_none() && best.0 < cut) {
            break if remaining_count <= 1 { f64::NAN } else { best.0 };
        }
        anyhow::ensure!(best.0 >= -1.0 || target.is_none(), "no eligible merge remains under --max-cluster-fraction; cannot reach --target-clusters");
        let (r, from, into) = best;
        let nf = counts[from];
        let ni = counts[into];
        if from < into {
            let (left, right) = sums.split_at_mut(into);
            for (dst, src) in right[0].iter_mut().zip(&left[from]) { *dst += *src; }
        } else {
            let (left, right) = sums.split_at_mut(from);
            for (dst, src) in left[into].iter_mut().zip(&right[0]) { *dst += *src; }
        }
        counts[into] += counts[from];
        counts[from] = 0;
        active[from] = false;
        remaining_count -= 1;
        let denom = counts[into] as f32;
        for (mean, sum) in means[into].iter_mut().zip(&sums[into]) {
            *mean = *sum / denom;
        }
        // Only the receiving cluster's profile changed.
        let updated = (0..n_clusters).into_par_iter().map(|j| {
            if active[j] && j != into {
                pearson_f32_noalloc(&means[into], &means[j])
            } else { -2.0 }
        }).collect::<Vec<_>>();
        for j in 0..n_clusters {
            correlations[into][j] = updated[j];
            correlations[j][into] = updated[j];
        }
        history.push(MergeStep {
            step: history.len() + 1,
            from: ids[from],
            into: ids[into],
            correlation: r,
            ward_delta_sse: None,
            from_cells: nf,
            into_cells: ni,
            merged_cells: nf + ni,
        });
        // Update cell assignments once after all merges, not once per merge.
        for slot in &mut cell_slots {
            if *slot == from { *slot = into; }
        }
    };
    for (label, &slot) in labels.iter_mut().zip(&cell_slots) {
        *label = ids[slot];
    }
    let remaining = unique(&labels);
    let map = remaining.iter().enumerate().map(|(i, &id)| (id, i))
        .collect::<std::collections::HashMap<_, _>>();
    for label in &mut labels { *label = map[label]; }
    Ok((labels, history, highest_remaining))
}

/// One conservative, non-recursive spatial refinement of merged clusters.
/// A split requires two sizeable groups separated by a large empty interval
/// along a PCA axis; elongated continuous populations remain intact.
pub(crate) fn split_disconnected_pca_groups(
    pca: &ndarray::Array2<f32>, labels: &[usize], min_side: usize,
) -> (Vec<usize>, usize) {
    let mut result = labels.to_vec();
    if labels.is_empty() || pca.ncols() == 0 { return (result, 0); }
    let mut next = labels.iter().copied().max().unwrap() + 1;
    let mut splits = 0;
    for id in unique(labels) {
        let cells: Vec<usize> = labels.iter().enumerate()
            .filter_map(|(i, &v)| (v == id).then_some(i)).collect();
        if cells.len() < 2 * min_side.max(2) { continue; }
        let mut best: Option<(f32, usize, Vec<usize>, usize)> = None;
        for axis in 0..pca.ncols() {
            let mut order = cells.clone();
            order.sort_unstable_by(|&a, &b| pca[(a, axis)].total_cmp(&pca[(b, axis)]).then(a.cmp(&b)));
            let range = pca[(*order.last().unwrap(), axis)] - pca[(order[0], axis)];
            if !range.is_finite() || range <= 0.0 { continue; }
            for cut in min_side.max(2)..=order.len() - min_side.max(2) {
                let gap = pca[(order[cut], axis)] - pca[(order[cut-1], axis)];
                let ratio = gap / range;
                if ratio >= 0.25 && best.as_ref().is_none_or(|(score, _, _, _)| ratio > *score) {
                    best = Some((ratio, axis, order.clone(), cut));
                }
            }
        }
        if let Some((_, _, order, cut)) = best {
            for &cell in &order[cut..] { result[cell] = next; }
            next += 1;
            splits += 1;
        }
    }
    let ids = unique(&result);
    let lookup: std::collections::HashMap<usize, usize> = ids.iter().enumerate()
        .map(|(new, &old)| (old, new)).collect();
    for label in &mut result { *label = lookup[label]; }
    (result, splits)
}

#[cfg(test)]
mod spatial_refinement_tests {
    use super::split_disconnected_pca_groups;
    use ndarray::Array2;
    #[test]
    fn splits_two_islands_but_not_continuous_cluster() {
        let islands = Array2::from_shape_vec((12, 1),
            (0..6).map(|i| i as f32 * 0.01)
                .chain((0..6).map(|i| 10.0 + i as f32 * 0.01)).collect()).unwrap();
        let (labels, splits) = split_disconnected_pca_groups(&islands, &[0; 12], 3);
        assert_eq!(splits, 1);
        assert_ne!(labels[0], labels[11]);
        let continuous = Array2::from_shape_vec((12, 1), (0..12).map(|i| i as f32).collect()).unwrap();
        assert_eq!(split_disconnected_pca_groups(&continuous, &[0; 12], 3).1, 0);
    }
}

// Pearson correlation in f64 without allocating converted vectors.
fn pearson_f32_noalloc(a: &[f32], b: &[f32]) -> f64 {
    if a.len() != b.len() || a.is_empty() { return 0.0; }
    let n = a.len() as f64;
    let ma = a.iter().map(|&v| v as f64).sum::<f64>() / n;
    let mb = b.iter().map(|&v| v as f64).sum::<f64>() / n;
    let (mut xy, mut aa, mut bb) = (0.0, 0.0, 0.0);
    for (&x, &y) in a.iter().zip(b) {
        let dx = x as f64 - ma;
        let dy = y as f64 - mb;
        xy += dx * dy;
        aa += dx * dx;
        bb += dy * dy;
    }
    let denom = (aa * bb).sqrt();
    if denom == 0.0 { 0.0 } else { xy / denom }
}

fn unique(labels: &[usize]) -> Vec<usize> {
    let mut v = labels.to_vec();
    v.sort_unstable();
    v.dedup();
    v
}

pub(crate) fn write_cluster_quality(
    stats: &std::path::Path, pca: &ndarray::Array2<f32>,
    hvg: &ndarray::Array2<f32>, labels: &[usize],
) -> anyhow::Result<()> {
    use std::io::Write;
    let k = labels.iter().copied().max().unwrap_or(0) + 1;
    let mut counts = vec![0usize; k];
    let mut pc_means = vec![vec![0.0f64; pca.ncols()]; k];
    let mut gene_means = vec![vec![0.0f32; hvg.ncols()]; k];
    for (i,&cl) in labels.iter().enumerate() {
        counts[cl] += 1;
        for d in 0..pca.ncols() { pc_means[cl][d] += pca[(i,d)] as f64; }
        for g in 0..hvg.ncols() { gene_means[cl][g] += hvg[(i,g)]; }
    }
    for cl in 0..k {
        let n = counts[cl].max(1) as f64;
        for v in &mut pc_means[cl] { *v /= n; }
        for v in &mut gene_means[cl] { *v /= n as f32; }
    }
    let mut pc_ss = vec![0.0f64; k];
    let mut rho_sum = vec![0.0f64; k];
    let mut rho_values = vec![Vec::<f64>::new(); k];
    for (i,&cl) in labels.iter().enumerate() {
        for d in 0..pca.ncols() {
            let delta = pca[(i,d)] as f64 - pc_means[cl][d];
            pc_ss[cl] += delta*delta;
        }
        let rho = pearson_f32_noalloc(hvg.row(i).as_slice().unwrap(), &gene_means[cl]);
        rho_sum[cl] += rho;
        rho_values[cl].push(rho);
    }
    let mut correlations = vec![vec![1.0f64; k]; k];
    let mut nearest = vec![-1.0f64; k];
    for a in 0..k { for b in (a+1)..k {
        let r = pearson_f32_noalloc(&gene_means[a], &gene_means[b]);
        correlations[a][b] = r; correlations[b][a] = r;
        nearest[a] = nearest[a].max(r); nearest[b] = nearest[b].max(r);
    }}
    std::fs::create_dir_all(stats)?;
    let mut out = std::io::BufWriter::new(std::fs::File::create(stats.join("cluster_quality.tsv"))?);
    writeln!(out, "cluster\tcells\tpca_variance_per_cell\tinternal_rho_mean\tinternal_rho_median\tinternal_rho_sd\tnearest_cluster_rho")?;
    for cl in 0..k {
        let n = counts[cl] as f64;
        let avg = rho_sum[cl] / n;
        rho_values[cl].sort_by(|a,b| a.total_cmp(b));
        let median = if counts[cl]%2 == 0 { (rho_values[cl][counts[cl]/2-1] + rho_values[cl][counts[cl]/2])*0.5 } else { rho_values[cl][counts[cl]/2] };
        let sd = (rho_values[cl].iter().map(|r| (r-avg)*(r-avg)).sum::<f64>()/n).sqrt();
        writeln!(out, "cluster_{cl:03}\t{}\t{:.6}\t{avg:.6}\t{median:.6}\t{sd:.6}\t{:.6}", counts[cl], pc_ss[cl]/n, nearest[cl])?;
    }
    let mut out = std::io::BufWriter::new(std::fs::File::create(stats.join("cluster_correlations.tsv"))?);
    write!(out, "cluster")?;
    for cl in 0..k { write!(out, "\tcluster_{cl:03}")?; }
    writeln!(out)?;
    for a in 0..k {
        write!(out, "cluster_{a:03}")?;
        for b in 0..k { write!(out, "\t{:.6}", correlations[a][b])?; }
        writeln!(out)?;
    }
    Ok(())
}

#[cfg(test)]
mod pearson_merge_tests {
    use super::merge_correlated;
    use sprs::TriMat;

    #[test]
    fn both_modes_use_full_expression_pearson() {
        let mut tri = TriMat::<f32>::new((3, 3));
        for (g, c, v) in [(0,0,3.),(1,0,1.),(0,1,4.),(1,1,1.),(2,2,5.)] {
            tri.add_triplet(g,c,v);
        }
        let matrix = tri.to_csr();
        let (labels, history, remaining) = merge_correlated(&matrix, &[0,1,2], vec![0,1,2], 0.0, Some(2), 0.0).unwrap();
        assert_eq!(history.len(), 1);
        assert_eq!(labels[0], labels[1]);
        assert_ne!(labels[0], labels[2]);
        assert!(remaining.is_finite());
        let (_, history, _) = merge_correlated(&matrix, &[0,1,2], vec![0,1,2], 0.9, None, 0.0).unwrap();
        assert_eq!(history.len(), 1);
        assert!(merge_correlated(&matrix, &[0,1,2], vec![0,1,2], 0.0, Some(4), 0.0).is_err());
    }
}
