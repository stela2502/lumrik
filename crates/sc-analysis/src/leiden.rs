//! Independent PCA-graph Leiden baseline (not Scanpy fuzzy-simplicial weighting).
use anyhow::{Context, Result};
use leiden_rs::{GraphDataBuilder, Leiden, LeidenConfig, QualityType};
use ndarray::Array2;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

// Exact hierarchical kNN over spatial microclusters. Each tree node bounds all
// its descendant cells with a centroid and enclosing radius. Branches whose
// lower bound cannot beat the current kth neighbor are skipped, including
// branches in other mega-clusters. No global cell-pair distance matrix.
#[derive(Clone)]
struct Node {
    center: Vec<f64>,
    radius: f64,
    children: Option<(usize, usize)>,
    cells: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
struct Candidate {
    dist2: f64,
    id: usize,
}
impl PartialEq for Candidate {
    fn eq(&self, rhs: &Self) -> bool {
        self.dist2 == rhs.dist2 && self.id == rhs.id
    }
}
impl Eq for Candidate {}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, rhs: &Self) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}
impl Ord for Candidate {
    fn cmp(&self, rhs: &Self) -> Ordering {
        self.dist2.total_cmp(&rhs.dist2).then(self.id.cmp(&rhs.id))
    }
}

fn squared(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}
fn make_node(points: &[Vec<f64>], members: Vec<usize>, children: Option<(usize, usize)>) -> Node {
    let dims = points[0].len();
    let mut center = vec![0.0; dims];
    for &i in &members {
        for (v, x) in center.iter_mut().zip(&points[i]) {
            *v += *x;
        }
    }
    for v in &mut center {
        *v /= members.len() as f64;
    }
    let radius = members
        .iter()
        .map(|&i| squared(&center, &points[i]).sqrt())
        .fold(0.0f64, f64::max);
    Node {
        center,
        radius,
        children,
        cells: if children.is_some() {
            Vec::new()
        } else {
            members
        },
    }
}

// Build a balanced hierarchy of microcluster bounding balls. Internal nodes
// also act as the proposed mega-clusters, recursively up to the root.
fn build_tree(
    points: &[Vec<f64>],
    groups: &[Vec<usize>],
    ids: Vec<usize>,
    nodes: &mut Vec<Node>,
) -> usize {
    if ids.len() == 1 {
        let id = nodes.len();
        nodes.push(make_node(points, groups[ids[0]].clone(), None));
        return id;
    }
    let dims = points[0].len();
    let mut centers = Vec::with_capacity(ids.len());
    for &group in &ids {
        let mut mean = vec![0.0; dims];
        for &cell in &groups[group] {
            for (v, x) in mean.iter_mut().zip(&points[cell]) {
                *v += *x;
            }
        }
        for v in &mut mean {
            *v /= groups[group].len() as f64;
        }
        centers.push(mean);
    }
    let axis = (0..dims)
        .max_by(|&a, &b| {
            let variance = |axis: usize| {
                let mean = centers.iter().map(|c| c[axis]).sum::<f64>() / centers.len() as f64;
                centers
                    .iter()
                    .map(|c| (c[axis] - mean).powi(2))
                    .sum::<f64>()
            };
            variance(a).total_cmp(&variance(b))
        })
        .unwrap_or(0);
    let mut order: Vec<usize> = (0..ids.len()).collect();
    order.sort_unstable_by(|&a, &b| {
        centers[a][axis]
            .total_cmp(&centers[b][axis])
            .then(ids[a].cmp(&ids[b]))
    });
    let mid = order.len() / 2;
    let left_ids = order[..mid].iter().map(|&i| ids[i]).collect();
    let right_ids = order[mid..].iter().map(|&i| ids[i]).collect();
    let left = build_tree(points, groups, left_ids, nodes);
    let right = build_tree(points, groups, right_ids, nodes);
    let mut members = Vec::new();
    fn gather(node: usize, nodes: &[Node], into: &mut Vec<usize>) {
        if let Some((a, b)) = nodes[node].children {
            gather(a, nodes, into);
            gather(b, nodes, into);
        } else {
            into.extend_from_slice(&nodes[node].cells);
        }
    }
    gather(left, nodes, &mut members);
    gather(right, nodes, &mut members);
    let id = nodes.len();
    nodes.push(make_node(points, members, Some((left, right))));
    id
}

fn lower_bound2(x: &[f64], node: &Node) -> f64 {
    let gap = (squared(x, &node.center).sqrt() - node.radius).max(0.0);
    gap * gap
}
fn search(
    x_id: usize,
    node_id: usize,
    points: &[Vec<f64>],
    nodes: &[Node],
    k: usize,
    best: &mut BinaryHeap<Candidate>,
) {
    let node = &nodes[node_id];
    if best.len() == k
        && lower_bound2(&points[x_id], node) > best.peek().unwrap().dist2 * (1.0 + 1e-12) + 1e-12
    {
        return;
    }
    if let Some((a, b)) = node.children {
        let da = lower_bound2(&points[x_id], &nodes[a]);
        let db = lower_bound2(&points[x_id], &nodes[b]);
        let (first, second) = if da <= db { (a, b) } else { (b, a) };
        search(x_id, first, points, nodes, k, best);
        search(x_id, second, points, nodes, k, best);
    } else {
        for &j in &node.cells {
            if j == x_id {
                continue;
            }
            let candidate = Candidate {
                dist2: squared(&points[x_id], &points[j]),
                id: j,
            };
            if best.len() < k {
                best.push(candidate);
            } else if candidate < *best.peek().unwrap() {
                best.pop();
                best.push(candidate);
            }
        }
    }
}

fn hierarchical_knn(
    pca: &Array2<f32>,
    labels: &[usize],
    k: usize,
) -> Result<Vec<Vec<(usize, f64)>>> {
    let n = pca.nrows();
    anyhow::ensure!(
        labels.len() == n,
        "microcluster labels must match PCA cell count"
    );
    anyhow::ensure!(
        pca.iter().all(|v| v.is_finite()),
        "PCA contains nonfinite coordinates"
    );
    let points: Vec<Vec<f64>> = pca
        .outer_iter()
        .map(|r| r.iter().map(|&v| v as f64).collect())
        .collect();
    let mut by_label = BTreeMap::<usize, Vec<usize>>::new();
    for (i, &label) in labels.iter().enumerate() {
        by_label.entry(label).or_default().push(i);
    }
    let groups: Vec<Vec<usize>> = by_label.into_values().collect();
    let mut nodes = Vec::with_capacity(groups.len() * 2);
    let root = build_tree(&points, &groups, (0..groups.len()).collect(), &mut nodes);
    eprintln!(
        "  Leiden kNN: {} cells, {} microclusters, {} hierarchy nodes, k={k}",
        n,
        groups.len(),
        nodes.len()
    );
    Ok((0..n)
        .into_par_iter()
        .map(|i| {
            let mut best = BinaryHeap::with_capacity(k + 1);
            search(i, root, &points, &nodes, k, &mut best);
            let mut row: Vec<(usize, f64)> =
                best.into_iter().map(|v| (v.id, v.dist2.sqrt())).collect();
            row.sort_unstable_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            row
        })
        .collect())
}

/// Construct a sparse cell graph using exact, spatially pruned PCA kNN.
/// Weights are locally scaled, not Scanpy's fuzzy-simplicial connectivities.
pub(crate) fn leiden_pca(
    pca: &Array2<f32>,
    microclusters: &[usize],
    k: usize,
    resolution: f64,
) -> Result<Vec<usize>> {
    anyhow::ensure!(
        resolution.is_finite() && resolution > 0.0,
        "Leiden resolution must be finite and positive"
    );
    let n = pca.nrows();
    anyhow::ensure!(
        n > 0 && pca.ncols() > 0,
        "Leiden needs nonempty PCA coordinates"
    );
    if n == 1 {
        return Ok(vec![0]);
    }
    let k = k.clamp(1, n - 1);
    let neighbors = hierarchical_knn(pca, microclusters, k)?;
    eprintln!("  Leiden kNN complete; constructing weighted graph...");
    let scales: Vec<f64> = neighbors
        .iter()
        .map(|v| v.last().unwrap().1.max(1e-9))
        .collect();
    let mut edges = BTreeMap::<(usize, usize), f64>::new();
    for (i, row) in neighbors.iter().enumerate() {
        for &(j, distance) in row {
            let pair = (i.min(j), i.max(j));
            let weight = (-distance / (scales[i] * scales[j]).sqrt())
                .exp()
                .max(1e-12);
            edges
                .entry(pair)
                .and_modify(|w| *w = w.max(weight))
                .or_insert(weight);
        }
    }
    let mut builder = GraphDataBuilder::new(n);
    for ((i, j), w) in edges {
        builder.add_edge(i, j, w).context("Leiden graph edge")?;
    }
    let graph = builder.build().context("Leiden graph construction")?;
    let config = LeidenConfig {
        resolution,
        quality: QualityType::RBConfiguration,
        seed: Some(42),
        ..Default::default()
    };
    let result = Leiden::new(config)
        .run(&graph)
        .context("Leiden optimization")?;
    let raw = result.partition.as_slice();
    anyhow::ensure!(
        raw.len() == n,
        "Leiden returned an incorrect number of assignments"
    );
    let mut ids = BTreeSet::new();
    ids.extend(raw.iter().copied());
    let remap = ids
        .into_iter()
        .enumerate()
        .map(|(new, old)| (old, new))
        .collect::<BTreeMap<_, _>>();
    Ok(raw.iter().map(|id| remap[id]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hierarchical_neighbors_match_bruteforce() {
        let x = ndarray::Array2::from_shape_fn((137, 5), |(i, j)| {
            ((i * 17 + j * 31) % 101) as f32 / 9.0 + (i / 17) as f32
        });
        let labels: Vec<usize> = (0..137).map(|i| i / 11).collect();
        let actual = hierarchical_knn(&x, &labels, 7).unwrap();
        for i in 0..x.nrows() {
            let mut expected: Vec<(usize, f64)> = (0..x.nrows())
                .filter(|&j| j != i)
                .map(|j| {
                    let d2 = (0..x.ncols())
                        .map(|c| ((x[(i, c)] as f64) - (x[(j, c)] as f64)).powi(2))
                        .sum::<f64>();
                    (j, d2.sqrt())
                })
                .collect();
            expected.sort_unstable_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
            assert_eq!(actual[i], expected[..7], "cell {i}");
        }
    }
    #[test]
    fn two_separated_groups() {
        let x = ndarray::arr2(&[
            [0.0, 0.0],
            [0.1, 0.0],
            [0.0, 0.1],
            [100.0, 100.0],
            [100.1, 100.0],
            [100.0, 100.1],
        ]);
        let labels = leiden_pca(&x, &[0, 0, 0, 1, 1, 1], 2, 1.0).unwrap();
        assert_eq!(labels.len(), 6);
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[3], labels[4]);
        assert_ne!(labels[0], labels[3]);
    }
}
