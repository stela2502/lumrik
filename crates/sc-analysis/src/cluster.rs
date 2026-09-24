use crate::stats::pearson_f32;
use sprs::CsMat;

#[derive(Debug, Clone)]
pub struct MergeStep {
    pub step: usize,
    pub from: usize,
    pub into: usize,
    pub correlation: f64,
    pub from_cells: usize,
    pub into_cells: usize,
    pub merged_cells: usize,
}

pub(crate) fn merge_correlated(
    x: &CsMat<f32>,
    mut labels: Vec<usize>,
    cut: f64,
) -> (Vec<usize>, Vec<MergeStep>) {
    let mut history = Vec::new();
    loop {
        let ids = unique(&labels);
        if ids.len() < 2 {
            break;
        }
        let means = means(x, &labels, &ids);
        let mut chosen = None;
        // Match ScanpyAutoAnalyzer mergeClosest: first qualifying cluster in stable order,
        // merged into its best-correlated partner, then recompute all means.
        for (ai, &a) in ids.iter().enumerate() {
            let mut best = (-2.0f64, 0usize);
            for (bi, &b) in ids.iter().enumerate() {
                if ai == bi {
                    continue;
                }
                let r = pearson_f32(&means[ai], &means[bi]);
                if r > best.0 {
                    best = (r, b);
                }
            }
            if best.0 >= cut {
                chosen = Some((a, best.1, best.0));
                break;
            }
        }
        let Some((from, into, r)) = chosen else {
            break;
        };
        let nf = labels.iter().filter(|&&v| v == from).count();
        let ni = labels.iter().filter(|&&v| v == into).count();
        for v in &mut labels {
            if *v == from {
                *v = into;
            }
        }
        history.push(MergeStep {
            step: history.len() + 1,
            from,
            into,
            correlation: r,
            from_cells: nf,
            into_cells: ni,
            merged_cells: nf + ni,
        });
    }
    let ids = unique(&labels);
    let map = ids
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .collect::<std::collections::HashMap<_, _>>();
    for v in &mut labels {
        *v = map[v];
    }
    (labels, history)
}
fn unique(labels: &[usize]) -> Vec<usize> {
    let mut v = labels.to_vec();
    v.sort_unstable();
    v.dedup();
    v
}
fn means(x: &CsMat<f32>, labels: &[usize], ids: &[usize]) -> Vec<Vec<f32>> {
    let map = ids
        .iter()
        .enumerate()
        .map(|(i, &v)| (v, i))
        .collect::<std::collections::HashMap<_, _>>();
    let mut out = vec![vec![0.0; x.rows()]; ids.len()];
    let mut n = vec![0usize; ids.len()];
    for &lab in labels {
        n[map[&lab]] += 1;
    }
    for (gene, row) in x.outer_iterator().enumerate() {
        for (cell, v) in row.iter() {
            let g = map[&labels[cell]];
            out[g][gene] += *v;
        }
    }
    for g in 0..out.len() {
        let d = n[g].max(1) as f32;
        for v in &mut out[g] {
            *v /= d;
        }
    }
    out
}
