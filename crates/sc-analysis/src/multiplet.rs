//! Cluster-mean mixture explanations in full linear expression space.
use crate::data::SingleCellData;
use crate::qc::QcMetrics;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Clone, Debug)]
pub(crate) struct Fit {
    pub(crate) explanation: String,
    pub(crate) components: Vec<(usize, f64)>,
    pub(crate) singlet_error: f64,
    pub(crate) mixture_error: f64,
    pub(crate) improvement: f64,
}

fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x,y)| x*y).sum() }
fn squared(a: &[f64]) -> f64 { dot(a,a) }
fn difference(a: &[f64], b: &[f64]) -> Vec<f64> { a.iter().zip(b).map(|(x,y)| x-y).collect() }
fn project(x: &[f64], a: &[f64], b: &[f64]) -> (f64, f64) {
    let ab = difference(a,b);
    let denom = squared(&ab);
    let w = if denom > 1e-12 { (dot(&difference(x,b), &ab)/denom).clamp(0.0,1.0) } else { 1.0 };
    let error = x.iter().zip(a).zip(b).map(|((&v,&aa),&bb)| (v-w*aa-(1.0-w)*bb).powi(2)).sum();
    (w,error)
}

/// Return one fit per cluster, using all surviving genes and original
/// surviving-UMI depths at the same linear 10,000 scale as user-facing data.
pub(crate) fn detect(raw: &SingleCellData, qc: &QcMetrics, labels: &[usize],
    min_parent_fraction: f64, min_improvement: f64) -> Vec<Fit> {
    let n = labels.len();
    if n == 0 { return Vec::new(); }
    let k = labels.iter().copied().max().unwrap_or(0) + 1;
    let mut counts = vec![0usize; k];
    for &c in labels { counts[c] += 1; }
    let mut means = vec![Vec::<f64>::new(); k];
    // Build full expression means in feature order, without a cell x gene
    // dense matrix. The normalization precedes the cluster average.
    for (g, row) in raw.matrix.outer_iterator().enumerate() {
        if !qc.keep_gene[g] { continue; }
        let mut sums = vec![0.0f64; k];
        for (cell, &v) in row.iter() {
            if qc.surviving[cell] > 0.0 {
                sums[labels[cell]] += 10_000.0 * v as f64 / qc.surviving[cell];
            }
        }
        for c in 0..k { means[c].push(sums[c] / counts[c].max(1) as f64); }
    }
    let min_cells = ((n as f64 * min_parent_fraction).ceil() as usize).max(1);
    let parents: Vec<usize> = (0..k).filter(|&c| counts[c] >= min_cells).collect();
    let mut cluster_fits = Vec::with_capacity(k);
    for target in 0..k {
        let unchanged = || Fit { explanation: format!("cluster_{target:03}"),
            components: Vec::new(), singlet_error: 0.0, mixture_error: 0.0, improvement: 0.0 };
        // Large clusters are reference populations, not targets for relabeling.
        if counts[target] == 0 || counts[target] >= min_cells || parents.len() < 2 {
            cluster_fits.push(unchanged());
            continue;
        }
        let x = &means[target];
        let baseline = parents.iter().map(|&c| squared(&difference(x, &means[c])))
            .fold(f64::INFINITY, f64::min);
        if !baseline.is_finite() || baseline <= 1e-12 {
            cluster_fits.push(unchanged());
            continue;
        }
        let mut best_error = baseline;
        let mut best = Vec::<(usize,f64)>::new();
        // Exhaustive pair search, followed by greedy forward selection up to
        // five populations. Re-optimize *all* weights after each addition.
        for (i, &a) in parents.iter().enumerate() {
            for &b in &parents[i+1..] {
                let (w, e) = project(x, &means[a], &means[b]);
                if w >= 0.15 && 1.0-w >= 0.15 && e < best_error {
                    best_error = e;
                    best = vec![(a,w), (b,1.0-w)];
                }
            }
        }
        while best.len() >= 2 && best.len() < 5 {
            let mut next = None;
            for &c in &parents {
                if best.iter().any(|(id,_)| *id == c) { continue; }
                let mut ids: Vec<usize> = best.iter().map(|(id,_)| *id).collect();
                ids.push(c);
                let (weights, err) = optimize_simplex(x, &means, &ids);
                if weights.iter().all(|&w| w >= 0.15) && err + 1e-10 < best_error {
                    next = Some((ids.into_iter().zip(weights).collect::<Vec<_>>(), err));
                    best_error = err;
                }
            }
            if let Some((fit, _)) = next { best = fit; } else { break; }
        }
        let improvement = (1.0 - best_error / baseline).max(0.0);
        if improvement < min_improvement { best.clear(); best_error = baseline; }
        let explanation = if best.is_empty() { format!("cluster_{target:03}") } else {
            format!("multiplet({})", best.iter().map(|(c,w)| format!("{w:.2} cluster_{c:03}"))
                .collect::<Vec<_>>().join(" + "))
        };
        cluster_fits.push(Fit { explanation, components: best,
            singlet_error: baseline, mixture_error: best_error, improvement });
    }
    cluster_fits
}

/// Coordinate-descent simplex fit. With <=5 parents and <=30 PCA dimensions,
/// this is tiny compared with the upstream PCA/UMAP workload.
fn optimize_simplex(x: &[f64], means: &[Vec<f64>], ids: &[usize]) -> (Vec<f64>, f64) {
    let m = ids.len();
    let mut weights = vec![1.0 / m as f64; m];
    for _ in 0..80 {
        let mut changed = 0.0f64;
        for i in 0..m {
            for j in i+1..m {
                let total = weights[i] + weights[j];
                if total <= 0.0 { continue; }
                let mut numerator = 0.0;
                let mut denom = 0.0;
                for t in 0..x.len() {
                    let diff = means[ids[i]][t] - means[ids[j]][t];
                    let mut rest = 0.0;
                    for h in 0..m {
                        if h != i && h != j { rest += weights[h] * means[ids[h]][t]; }
                    }
                    numerator += diff * (x[t] - rest - total * means[ids[j]][t]);
                    denom += diff * diff;
                }
                let new_i = if denom > 1e-12 { (numerator / denom).clamp(0.0, total) } else { weights[i] };
                changed = changed.max((new_i - weights[i]).abs());
                weights[i] = new_i;
                weights[j] = total - new_i;
            }
        }
        if changed < 1e-8 { break; }
    }
    let error = (0..x.len()).map(|t| {
        let predicted: f64 = ids.iter().enumerate().map(|(i,&id)| weights[i]*means[id][t]).sum();
        (x[t]-predicted).powi(2)
    }).sum();
    (weights,error)
}

pub(crate) fn write_clusters(path: &Path, labels: &[usize], fits: &[Fit]) -> anyhow::Result<()> {
    let k = labels.iter().copied().max().map_or(0, |v| v+1);
    let mut counts = vec![0usize; k];
    for &c in labels { counts[c] += 1; }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "cluster\tn_cells\texplanation\tcomponents\tweights\tsinglet_error\tmixture_error\trelative_improvement")?;
    for (c, f) in fits.iter().enumerate() {
        if counts[c] > 0 {
            let ids = f.components.iter().map(|(id,_)| format!("cluster_{id:03}")).collect::<Vec<_>>().join(",");
            let weights = f.components.iter().map(|(_,v)| format!("{v:.6}")).collect::<Vec<_>>().join(",");
            writeln!(w, "cluster_{c:03}\t{}\t{}\t{ids}\t{weights}\t{:.8}\t{:.8}\t{:.8}",
                counts[c], f.explanation, f.singlet_error, f.mixture_error, f.improvement)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mixture_projection() {
        let (w,e)=project(&[3.0,7.0], &[10.0,0.0], &[0.0,10.0]);
        assert!((w-0.3).abs()<1e-9);
        assert!(e<1e-9);
    }
    #[test]
    fn cluster_means_recover_linear_mixture_not_cell_fits() {
        use sprs::TriMat;
        let mut tri = TriMat::<f32>::new((2, 23));
        for c in 0..10 { tri.add_triplet(0, c, 100.0); }
        for c in 10..20 { tri.add_triplet(1, c, 100.0); }
        for c in 20..23 {
            tri.add_triplet(0, c, 30.0);
            tri.add_triplet(1, c, 70.0);
        }
        let data = SingleCellData::new(tri.to_csr(), vec!["G1".into(), "G2".into()],
            (0..23).map(|i| format!("cell{i}")).collect(),
            crate::data::CellAnnotations::from_columns(vec![], 23).unwrap()).unwrap();
        let qc = crate::qc::qc(&data, false);
        let labels: Vec<usize> = (0..23).map(|i| if i < 10 { 0 } else if i < 20 { 1 } else { 2 }).collect();
        let fits = detect(&data, &qc, &labels, 0.2, 0.35);
        assert_eq!(fits.len(), 3, "must return clusters, never cells");
        assert!(fits[0].components.is_empty() && fits[1].components.is_empty());
        assert_eq!(fits[2].components.len(), 2);
        assert!(fits[2].improvement > 0.99);
        assert!((fits[2].components[0].1 - 0.3).abs() < 1e-8);
    }
    #[test]
    fn cluster_report_has_no_cell_column() {
        let path = std::env::temp_dir().join(format!("lumrik_cluster_multiplet_{}.tsv", std::process::id()));
        let unchanged = Fit { explanation: "cluster_000".into(), components: vec![], singlet_error: 0.0, mixture_error: 0.0, improvement: 0.0 };
        write_clusters(&path, &[0,0,1], &[unchanged.clone(), Fit { explanation: "cluster_001".into(), ..unchanged }]).unwrap();
        let report = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(path).unwrap();
        assert_eq!(report.lines().count(), 3);
        assert!(report.starts_with("cluster\tn_cells\t"));
        assert!(!report.contains("cell\tcluster"));
    }
}
