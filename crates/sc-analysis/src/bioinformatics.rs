use crate::cluster::merge_correlated;
use crate::data::{CellAnnotations, SingleCellData};
use crate::embedding::umap_2d;
use crate::mex::load_mex;
use crate::normalize::normalize_surviving;
use crate::norn::{load_beacon_blocks, resolve_exonic};
use crate::partition::spatial_patches;
use crate::pca::pca;
use crate::qc::qc;
use crate::report::{write_all, write_plots};
use anyhow::{Result, bail};
use ndarray::Array2;
use sprs::TriMat;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

#[derive(Debug, Clone)]
pub struct AnalysisConfig {
    pub min_umi_count: Option<usize>,
    pub exclude_vdj_before_normalization: bool,
    pub variable_genes: usize,
    pub pca_components: usize,
    pub umap_components: usize,
    pub umap_neighbors: usize,
    pub umap_epochs: usize,
    pub overclusters: usize,
    pub kmeans_iterations: usize,
    pub merge_pearson: f64,
    pub holdout_fraction: f64,
    pub normalization_scale: f64,
    pub pseudo_samples: usize,
    pub output_dir: PathBuf,
}
impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            min_umi_count: None,
            exclude_vdj_before_normalization: false,
            variable_genes: 2000,
            pca_components: 30,
            umap_components: 10,
            umap_neighbors: 15,
            umap_epochs: 200,
            overclusters: 30,
            kmeans_iterations: 50,
            merge_pearson: 0.95,
            holdout_fraction: 0.10,
            normalization_scale: 10_000.0,
            pseudo_samples: 5,
            output_dir: PathBuf::from("sc_analysis"),
        }
    }
}
#[derive(Debug, Clone)]
pub struct AnalysisSummary {
    pub input_cells: usize,
    pub filtered_cells: usize,
    pub input_features: usize,
    pub surviving_features: usize,
    pub variable_genes: usize,
    pub initial_clusters: usize,
    pub merged_clusters: usize,
    pub cluster_sizes: Vec<usize>,
    pub holdout_cells: usize,
    pub mitochondrial_features: usize,
    pub ribosomal_features: usize,
    pub vdj_segment_features: usize,
    pub timings: Vec<AnalysisTiming>,
}

#[derive(Debug, Clone)]
pub struct AnalysisTiming {
    pub stage: String,
    pub elapsed: Duration,
}

fn format_elapsed(elapsed: Duration) -> String {
    let seconds = elapsed.as_secs_f64();
    if seconds >= 60.0 {
        format!(
            "{}m {:.1}s",
            (seconds / 60.0).floor() as u64,
            seconds % 60.0
        )
    } else {
        format!("{seconds:.2}s")
    }
}

fn finish_stage(timings: &mut Vec<AnalysisTiming>, stage: &str, started: Instant) {
    let elapsed = started.elapsed();
    eprintln!("Finished {stage} in {}", format_elapsed(elapsed));
    timings.push(AnalysisTiming {
        stage: stage.to_string(),
        elapsed,
    });
}

fn subset_cells(data: &SingleCellData, keep: &[usize]) -> Result<SingleCellData> {
    let mut old_to_new = vec![usize::MAX; data.n_cells()];
    for (new, &old) in keep.iter().enumerate() {
        old_to_new[old] = new;
    }

    let mut triplets = TriMat::<f32>::new((data.n_features(), keep.len()));
    for (gene, row) in data.matrix.outer_iterator().enumerate() {
        for (old_cell, value) in row.iter() {
            let new_cell = old_to_new[old_cell];
            if new_cell != usize::MAX {
                triplets.add_triplet(gene, new_cell, *value);
            }
        }
    }

    let cells = keep
        .iter()
        .map(|&i| data.cells[i].clone())
        .collect::<Vec<_>>();
    let annotations = CellAnnotations::from_columns(vec![("cell", cells.clone())], cells.len())?;
    SingleCellData::new(triplets.to_csr(), data.features.clone(), cells, annotations)
}

pub fn analyze_exon_matrix(
    path: impl AsRef<Path>,
    config: &AnalysisConfig,
) -> Result<(SingleCellData, AnalysisSummary)> {
    if !(0.0..1.0).contains(&config.holdout_fraction) {
        bail!("holdout_fraction must be in [0,1)");
    }
    if !(-1.0..=1.0).contains(&config.merge_pearson) {
        bail!("merge_pearson must be in [-1,1]");
    }
    let mut timings = Vec::new();

    eprintln!("Starting matrix loading...");
    let started = Instant::now();
    let exonic = resolve_exonic(path.as_ref())?;
    let raw = load_mex(&exonic)?;
    finish_stage(&mut timings, "matrix loading", started);
    let input_cells = raw.n_cells();
    let n_genes = raw.n_features();
    if input_cells < 4 || n_genes < 2 {
        bail!("matrix too small: {n_genes} genes x {input_cells} cells");
    }

    // sc-analyse deliberately does not perform automatic cell calling.  This
    // optional threshold is an explicit exploratory filter so a permissive
    // Nelrune quantification can be re-cut quickly without re-quantifying.
    eprintln!("Starting UMI filtering...");
    let started = Instant::now();
    let raw = if let Some(min_umi_count) = config.min_umi_count {
        let raw_qc = qc(&raw, false);
        let keep = raw_qc
            .total
            .iter()
            .enumerate()
            .filter_map(|(i, &umi)| (umi >= min_umi_count as f64).then_some(i))
            .collect::<Vec<_>>();
        subset_cells(&raw, &keep)?
    } else {
        raw
    };
    finish_stage(&mut timings, "UMI filtering", started);
    let n_cells = raw.n_cells();
    if n_cells < 4 {
        bail!("min_umi_count retained only {n_cells} cells; at least 4 are required");
    }
    eprintln!("Starting QC and normalization...");
    let started = Instant::now();
    let q = qc(&raw, config.exclude_vdj_before_normalization);
    let mito_features = raw
        .features
        .iter()
        .filter(|x| crate::qc::classify_gene(x).0)
        .count();
    let ribo_features = raw
        .features
        .iter()
        .filter(|x| crate::qc::classify_gene(x).1)
        .count();
    let vdj_segment_features = raw
        .features
        .iter()
        .filter(|x| crate::qc::classify_gene(x).2)
        .count();
    let norm = normalize_surviving(&raw, &q, config.normalization_scale)?;
    finish_stage(&mut timings, "QC and normalization", started);
    eprintln!("Starting variable-gene preparation...");
    let started = Instant::now();
    let ng = norm.n_features();
    let mut variances = Vec::with_capacity(ng);
    for (g, row) in norm.matrix.outer_iterator().enumerate() {
        if crate::qc::is_vdj_segment_gene(&norm.features[g]) {
            continue;
        }
        let mut s = 0.0;
        let mut s2 = 0.0;
        for (_, v) in row.iter() {
            let x = *v as f64;
            s += x;
            s2 += x * x;
        }
        let m = s / n_cells as f64;
        variances.push(((s2 / n_cells as f64 - m * m).max(0.0), g));
    }
    variances.sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    let selected = variances
        .iter()
        .take(config.variable_genes)
        .map(|x| x.1)
        .collect::<Vec<_>>();
    let mut x = Array2::<f32>::zeros((n_cells, selected.len()));
    for (j, &g) in selected.iter().enumerate() {
        if let Some(row) = norm.matrix.outer_view(g) {
            for (c, v) in row.iter() {
                x[(c, j)] = *v;
            }
        }
        let m = (0..n_cells).map(|c| x[(c, j)] as f64).sum::<f64>() / n_cells as f64;
        let sd = ((0..n_cells)
            .map(|c| {
                let z = x[(c, j)] as f64 - m;
                z * z
            })
            .sum::<f64>()
            / n_cells as f64)
            .sqrt();
        if sd > 0.0 {
            for c in 0..n_cells {
                x[(c, j)] = ((x[(c, j)] as f64 - m) / sd) as f32;
            }
        }
    }
    finish_stage(&mut timings, "variable-gene preparation", started);

    eprintln!("Starting PCA...");
    let started = Instant::now();
    let (pca_coords, pca_var) = pca(&x, config.pca_components)?;
    finish_stage(&mut timings, "PCA", started);
    // Initial groups are deliberately geometric micro-patches, not biological
    // clusters.  Recursively split PCA space until each patch contains at most
    // max(50 cells, 1% of the retained dataset), then let the existing
    // mean-expression Pearson merge determine the final biological groups.
    eprintln!("Starting spatial partitioning...");
    let started = Instant::now();
    let initial = spatial_patches(&pca_coords, &pca_var);
    finish_stage(&mut timings, "spatial partitioning", started);
    let initial_k = initial.iter().copied().max().unwrap_or(0) + 1;
    eprintln!("Starting mean-profile merge...");
    let started = Instant::now();
    let (labels, history) = merge_correlated(&norm.matrix, initial.clone(), config.merge_pearson);
    finish_stage(&mut timings, "mean-profile merge", started);
    let merged_k = labels.iter().copied().max().unwrap_or(0) + 1;
    let mut cluster_sizes = vec![0usize; merged_k];
    for &cluster in &labels {
        cluster_sizes[cluster] += 1;
    }
    eprintln!("Starting UMAP...");
    let started = Instant::now();
    let viz = umap_2d(&pca_coords, config.umap_neighbors, config.umap_epochs)?;
    finish_stage(&mut timings, "UMAP", started);
    eprintln!("Starting Beacon integration...");
    let started = Instant::now();
    let beacon = load_beacon_blocks(&exonic, &norm.cells)?;
    finish_stage(&mut timings, "Beacon integration", started);
    eprintln!("Starting pseudo-deep holdout...");
    let started = Instant::now();
    let mut by = HashMap::<usize, Vec<usize>>::new();
    for (c, &cl) in labels.iter().enumerate() {
        by.entry(cl).or_default().push(c);
    }
    let mut holdout = vec![false; n_cells];
    for cells in by.values_mut() {
        cells.sort_unstable_by(|&a, &b| q.surviving[b].total_cmp(&q.surviving[a]).then(a.cmp(&b)));
        let n = ((cells.len() as f64 * config.holdout_fraction).ceil() as usize)
            .max(1)
            .min(cells.len().saturating_sub(1));
        for &c in cells.iter().take(n) {
            holdout[c] = true;
        }
    }
    finish_stage(&mut timings, "pseudo-deep holdout", started);

    eprintln!("Starting output writing...");
    let started = Instant::now();
    write_all(
        &config.output_dir,
        &norm,
        &q,
        &labels,
        &initial,
        &holdout,
        &pca_coords,
        &pca_var,
        &viz,
        &history,
        &beacon,
        config.pseudo_samples,
    )?;
    finish_stage(&mut timings, "output writing", started);

    eprintln!("Starting plot generation...");
    let started = Instant::now();
    write_plots(
        &config.output_dir,
        &norm,
        &q,
        &labels,
        &pca_var,
        &viz,
        &beacon,
    )?;
    finish_stage(&mut timings, "plot generation", started);

    let annotations = CellAnnotations::from_columns(
        vec![
            ("cell", norm.cells.clone()),
            (
                "cluster",
                labels.iter().map(|x| format!("cluster_{x:03}")).collect(),
            ),
            (
                "tech",
                holdout
                    .iter()
                    .map(|&x| {
                        if x {
                            "pseudo_deep".into()
                        } else {
                            "training".into()
                        }
                    })
                    .collect(),
            ),
            (
                "umi_count",
                q.total.iter().map(|x| format!("{x:.0}")).collect(),
            ),
            (
                "surviving_umi",
                q.surviving.iter().map(|x| format!("{x:.0}")).collect(),
            ),
            (
                "mitochondrial_fraction",
                (0..n_cells)
                    .map(|i| {
                        format!(
                            "{:.6}",
                            if q.total[i] > 0.0 {
                                q.mito[i] / q.total[i]
                            } else {
                                0.0
                            }
                        )
                    })
                    .collect(),
            ),
            (
                "ribosomal_fraction",
                (0..n_cells)
                    .map(|i| {
                        format!(
                            "{:.6}",
                            if q.total[i] > 0.0 {
                                q.ribo[i] / q.total[i]
                            } else {
                                0.0
                            }
                        )
                    })
                    .collect(),
            ),
            (
                "vdj_fraction",
                (0..n_cells)
                    .map(|i| {
                        format!(
                            "{:.6}",
                            if q.total[i] > 0.0 {
                                q.vdj[i] / q.total[i]
                            } else {
                                0.0
                            }
                        )
                    })
                    .collect(),
            ),
        ],
        n_cells,
    )?;
    let data = SingleCellData::new(norm.matrix, norm.features, norm.cells, annotations)?;
    Ok((
        data,
        AnalysisSummary {
            input_cells,
            filtered_cells: n_cells,
            input_features: n_genes,
            surviving_features: ng,
            variable_genes: selected.len(),
            initial_clusters: initial_k,
            merged_clusters: merged_k,
            cluster_sizes,
            holdout_cells: holdout.iter().filter(|&&x| x).count(),
            mitochondrial_features: mito_features,
            ribosomal_features: ribo_features,
            vdj_segment_features,
            timings,
        },
    ))
}
