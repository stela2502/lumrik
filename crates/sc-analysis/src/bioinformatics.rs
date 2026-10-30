use crate::cluster::{merge_correlated, write_cluster_quality};
use crate::data::{CellAnnotations, SingleCellData};
use crate::embedding::umap_2d;
use crate::mex::load_mex;
use crate::normalize::{normalize_surviving, resample_surviving};
use crate::norn::{load_beacon_blocks, load_clonomap, resolve_exonic};
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
    pub target_gene_list: Option<PathBuf>,
    pub exclude_vdj_before_normalization: bool,
    pub variable_genes: usize,
    pub pca_components: usize,
    pub umap_neighbors: usize,
    pub leiden_resolution: f64,
    pub leiden_neighbors: usize,
    pub umap_epochs: usize,
    pub spatial_min_cells: usize,
    pub spatial_fraction: f64,
    pub merge_pearson: f64,
    pub ward_fraction_max: f64,
    pub max_cluster_fraction: f64,
    pub target_clusters: Option<usize>,
    pub holdout_fraction: f64,
    pub pseudo_samples: usize,
    pub multiplet_min_parent_fraction: f64,
    pub multiplet_min_improvement: f64,
    pub output_dir: PathBuf,
}
impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            min_umi_count: Some(800),
            target_gene_list: None,
            exclude_vdj_before_normalization: false,
            variable_genes: 2000,
            pca_components: 30,
            umap_neighbors: 15,
            leiden_resolution: 1.0,
            leiden_neighbors: 15,
            umap_epochs: 200,
            spatial_min_cells: 50,
            spatial_fraction: 0.01,
            merge_pearson: 0.95,
            ward_fraction_max: 0.2,
            max_cluster_fraction: 0.10,
            target_clusters: None,
            holdout_fraction: 0.10,
            pseudo_samples: 5,
            multiplet_min_parent_fraction: 0.01,
            multiplet_min_improvement: 0.35,
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
    let annotations = data.annotations.subset(keep)?;
    SingleCellData::new(triplets.to_csr(), data.features.clone(), cells, annotations)
}

/// Restrict the raw feature-by-cell matrix before QC or any normalization.
fn filter_target_genes(data: SingleCellData, path: &Path) -> Result<SingleCellData> {
    use std::collections::HashSet;
    let contents = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("cannot read target gene list {}: {e}", path.display()))?;
    let targets: HashSet<String> = contents
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| line.to_string())
        .collect();
    anyhow::ensure!(
        !targets.is_empty(),
        "target gene list {} is empty",
        path.display()
    );
    let keep: Vec<usize> = data
        .features
        .iter()
        .enumerate()
        .filter_map(|(i, name)| targets.contains(name).then_some(i))
        .collect();
    anyhow::ensure!(
        !keep.is_empty(),
        "none of {} target genes in {} match matrix feature names (check symbols/IDs and species)",
        targets.len(),
        path.display()
    );
    eprintln!(
        "  Target gene filter: {} unique requested, {} matched matrix features, {} excluded",
        targets.len(),
        keep.len(),
        data.n_features() - keep.len()
    );
    let mut triplets = TriMat::<f32>::new((keep.len(), data.n_cells()));
    for (new_gene, &old_gene) in keep.iter().enumerate() {
        if let Some(row) = data.matrix.outer_view(old_gene) {
            for (cell, value) in row.iter() {
                triplets.add_triplet(new_gene, cell, *value);
            }
        }
    }
    let features = keep.iter().map(|&i| data.features[i].clone()).collect();
    SingleCellData::new(triplets.to_csr(), features, data.cells, data.annotations)
}

fn load_input(path: &Path) -> Result<SingleCellData> {
    let is_h5ad = path.extension().is_some_and(|ext| ext == "h5ad");
    let resolved = if is_h5ad {
        path.to_path_buf()
    } else {
        resolve_exonic(path)?
    };
    if is_h5ad {
        SingleCellData::from_h5ad(&resolved)
    } else {
        load_mex(&resolved)
    }
}

/// Merge independent samples by feature identifier, preserving their original barcodes
/// in `cellid` and using `<sample>_<cellid>` as experiment-wide cell labels.
pub fn merge_sample_inputs(paths: &[PathBuf], names: &[String]) -> Result<SingleCellData> {
    use std::collections::{HashMap, HashSet};
    anyhow::ensure!(
        !paths.is_empty() && paths.len() == names.len(),
        "--data and --sample-names require equal nonzero lengths"
    );
    let mut seen_names = HashSet::new();
    for name in names {
        anyhow::ensure!(
            !name.is_empty() && !name.contains(['\t', '\n']),
            "invalid sample name {name:?}"
        );
        anyhow::ensure!(seen_names.insert(name), "duplicate sample name {name:?}");
    }
    let mut datasets = Vec::new();
    let mut features = Vec::new();
    let mut lookup = HashMap::<String, usize>::new();
    let mut cells = Vec::new();
    let mut samples = Vec::new();
    let mut barcodes = Vec::new();
    let mut seen_cells = HashSet::new();
    for (path, name) in paths.iter().zip(names) {
        let data = load_input(path)?;
        let mut local = HashSet::new();
        for cell in &data.cells {
            anyhow::ensure!(
                local.insert(cell),
                "duplicate barcode {cell:?} in sample {name:?}"
            );
            let qualified = format!("{name}_{cell}");
            anyhow::ensure!(
                seen_cells.insert(qualified.clone()),
                "ambiguous combined cell identifier {qualified:?}"
            );
            cells.push(qualified);
            samples.push(name.clone());
            barcodes.push(cell.clone());
        }
        for feature in &data.features {
            if !lookup.contains_key(feature) {
                let idx = features.len();
                lookup.insert(feature.clone(), idx);
                features.push(feature.clone());
            }
        }
        datasets.push(data);
    }
    let nnz = datasets.iter().map(|d| d.matrix.nnz()).sum();
    let mut triplets = TriMat::<f32>::with_capacity((features.len(), cells.len()), nnz);
    let mut offset = 0;
    for data in datasets {
        for (row_idx, row) in data.matrix.outer_iterator().enumerate() {
            let target = lookup[&data.features[row_idx]];
            for (col, value) in row.iter() {
                triplets.add_triplet(target, offset + col, *value);
            }
        }
        offset += data.n_cells();
    }
    let annotations = crate::data::CellAnnotations::from_columns(
        vec![
            ("cell", cells.clone()),
            ("sample", samples),
            ("cellid", barcodes),
        ],
        cells.len(),
    )?;
    SingleCellData::new(triplets.to_csr(), features, cells, annotations)
}

pub fn analyze_sample_matrices(
    paths: &[PathBuf],
    names: &[String],
    config: &AnalysisConfig,
) -> Result<(SingleCellData, AnalysisSummary)> {
    let merged = merge_sample_inputs(paths, names)?;
    analyze_loaded_matrix(merged, config, None)
}

pub fn analyze_exon_matrix(
    path: impl AsRef<Path>,
    config: &AnalysisConfig,
) -> Result<(SingleCellData, AnalysisSummary)> {
    let raw = load_input(path.as_ref())?;
    analyze_loaded_matrix(raw, config, Some(path.as_ref()))
}

fn analyze_loaded_matrix(
    raw: SingleCellData,
    config: &AnalysisConfig,
    source: Option<&Path>,
) -> Result<(SingleCellData, AnalysisSummary)> {
    // The target panel is the FIRST filter, before any QC, cell filtering, or normalization.
    let raw = if let Some(path) = &config.target_gene_list {
        filter_target_genes(raw, path)?
    } else {
        raw
    };
    if !(0.0..1.0).contains(&config.holdout_fraction) {
        bail!("holdout_fraction must be in [0,1)");
    }
    anyhow::ensure!(
        config.ward_fraction_max.is_finite() && (0.0..=1.0).contains(&config.ward_fraction_max),
        "ward_fraction_max must be in [0,1]"
    );
    if !(-1.0..=1.0).contains(&config.merge_pearson) {
        bail!("merge_pearson must be in [-1,1]");
    }
    anyhow::ensure!(
        config.min_umi_count.is_some_and(|n| n > 0),
        "min_umi_count must be positive"
    );
    anyhow::ensure!(
        config.max_cluster_fraction.is_finite()
            && (0.0..=1.0).contains(&config.max_cluster_fraction),
        "max_cluster_fraction must be in [0,1]"
    );
    let mut timings = Vec::new();

    eprintln!("Starting matrix loading...");
    let started = Instant::now();
    finish_stage(&mut timings, "matrix loading", started);
    let input_cells = raw.n_cells();
    let n_genes = raw.n_features();
    if input_cells < 4 || n_genes < 2 {
        bail!("matrix too small: {n_genes} genes x {input_cells} cells");
    }

    // Filter on surviving genes at the exact depth required for resampling.
    // This is an explicit analysis eligibility threshold, not cell calling.
    eprintln!("Starting UMI filtering...");
    let started = Instant::now();
    // Eligibility uses the same surviving genes that enter fixed-depth resampling.
    let raw_qc = qc(&raw, config.exclude_vdj_before_normalization);
    let threshold = config.min_umi_count.expect("validated minimum UMI depth");
    let mut keep = Vec::with_capacity(input_cells);
    let mut removed_umi = Vec::new();
    let mut retained_umi = Vec::new();
    for (i, &umi) in raw_qc.surviving.iter().enumerate() {
        if umi >= threshold as f64 {
            keep.push(i);
            retained_umi.push(umi);
        } else {
            removed_umi.push(umi);
        }
    }
    let retained_count = keep.len();
    let removed_count = input_cells - retained_count;
    let median = |values: &mut Vec<f64>| -> Option<f64> {
        if values.is_empty() {
            return None;
        }
        values.sort_by(f64::total_cmp);
        let middle = values.len() / 2;
        Some(if values.len() % 2 == 0 {
            (values[middle - 1] + values[middle]) / 2.0
        } else {
            values[middle]
        })
    };
    let retained_median = median(&mut retained_umi);
    let removed_median = median(&mut removed_umi);
    let median_text =
        |value: Option<f64>| value.map_or_else(|| "NA".to_string(), |v| format!("{v:.1}"));
    eprintln!("  Minimum surviving UMIs/cell and resampling depth: {threshold}");
    eprintln!("  Input cells:             {input_cells}");
    eprintln!("  Retained cells:          {retained_count}");
    eprintln!("  Removed cells:           {removed_count}");
    eprintln!(
        "  Retained fraction:       {:.2}%",
        100.0 * retained_count as f64 / input_cells as f64
    );
    eprintln!(
        "  Median UMIs retained:    {}",
        median_text(retained_median)
    );
    eprintln!("  Median UMIs removed:     {}", median_text(removed_median));
    // Write diagnostics before expensive downstream analysis.
    // No need to rebuild the matrix if every cell is retained.
    let qc_dir = config.output_dir.join("qc");
    std::fs::create_dir_all(&qc_dir)?;
    let summary = format!(
        "metric\tvalue\nminimum_surviving_umi\t{}\ninput_cells\t{input_cells}\nretained_cells\t{retained_count}\nremoved_cells\t{removed_count}\nretained_fraction\t{:.6}\nmedian_umi_retained\t{}\nmedian_umi_removed\t{}\n",
        threshold,
        retained_count as f64 / input_cells as f64,
        median_text(retained_median),
        median_text(removed_median),
    );
    std::fs::write(qc_dir.join("umi_filter_summary.tsv"), summary)?;
    let raw = if removed_count == 0 {
        raw
    } else {
        subset_cells(&raw, &keep)?
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
    anyhow::ensure!(
        q.surviving.iter().all(|&v| v >= threshold as f64),
        "retained cell has fewer surviving UMIs than the resampling depth"
    );
    let norm = normalize_surviving(&raw, &q, 10_000.0)?;
    let sampled = resample_surviving(&raw, &q, threshold, 42)?;
    // Keep exact-depth counts for exports; use log1p counts for geometry and Pearson merging.
    let mut sampled_log = sampled.matrix.clone();
    for value in sampled_log.data_mut() {
        *value = value.ln_1p();
    }
    eprintln!("  PCA/UMAP/Pearson: log1p of {threshold} resampled surviving UMIs/cell, seed 42");
    eprintln!("  User-facing expression: log1p(10,000 * surviving_gene_umi / surviving_cell_umi)");
    finish_stage(&mut timings, "QC and normalization", started);
    eprintln!("Starting variable-gene preparation...");
    let started = Instant::now();
    let ng = sampled.n_features();
    let mut variances = Vec::with_capacity(ng);
    for (g, row) in sampled_log.outer_iterator().enumerate() {
        if crate::qc::is_vdj_segment_gene(&sampled.features[g]) {
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
        if let Some(row) = sampled_log.outer_view(g) {
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
    // max(configured minimum cells, configured fraction of retained cells).
    // Leiden graph communities then determine the biological groups.
    eprintln!("Starting spatial partitioning...");
    let started = Instant::now();
    let initial = spatial_patches(
        &pca_coords,
        &pca_var,
        config.spatial_min_cells,
        config.spatial_fraction,
    );
    finish_stage(&mut timings, "spatial partitioning", started);
    let initial_k = initial.iter().copied().max().unwrap_or(0) + 1;
    let mut initial_sizes = vec![0usize; initial_k];
    for &cluster in &initial {
        initial_sizes[cluster] += 1;
    }
    let initial_max = initial_sizes.iter().copied().max().unwrap_or(0);
    eprintln!(
        "  Initial spatial clusters: {initial_k} (largest: {initial_max} / {n_cells} cells, {:.1}%)",
        100.0 * initial_max as f64 / n_cells.max(1) as f64
    );
    // Spatial patches are a search index only, never biological cluster labels.
    eprintln!("Starting Leiden PCA graph and community detection...");
    let started = Instant::now();
    let leiden_initial = crate::leiden::leiden_pca(
        &pca_coords,
        &initial,
        config.leiden_neighbors,
        config.leiden_resolution,
    )?;
    finish_stage(&mut timings, "Leiden graph and communities", started);
    let leiden_k = leiden_initial.iter().copied().max().unwrap_or(0) + 1;
    eprintln!("  Leiden initial communities: {leiden_k}");
    eprintln!("Starting Leiden mean-profile Pearson merge...");
    let started = Instant::now();
    let all_genes = (0..norm.matrix.rows()).collect::<Vec<_>>();
    let (labels, history, highest_remaining) = merge_correlated(
        &norm.matrix,
        &all_genes,
        leiden_initial.clone(),
        config.merge_pearson,
        config.target_clusters,
        config.max_cluster_fraction,
    )?;
    finish_stage(&mut timings, "Leiden Pearson merge", started);
    let merged_k = labels.iter().copied().max().unwrap_or(0) + 1;
    let mut cluster_sizes = vec![0usize; merged_k];
    for &cluster in &labels {
        cluster_sizes[cluster] += 1;
    }
    let largest = cluster_sizes.iter().copied().max().unwrap_or(0);
    eprintln!(
        "  Leiden communities: {leiden_k} -> {merged_k} ({} Pearson merges)",
        history.len()
    );
    eprintln!("  Highest remaining Pearson rho: {highest_remaining:.6}");
    std::fs::create_dir_all(&config.output_dir)?;
    std::fs::write(
        config.output_dir.join("cluster_merge_stopping.tsv"),
        format!(
            "mode\tfinal_clusters\thighest_remaining_rho\n{}\t{}\t{:.8}\n",
            if config.target_clusters.is_some() {
                "target"
            } else {
                "threshold"
            },
            merged_k,
            highest_remaining
        ),
    )?;
    eprintln!(
        "  Largest final cluster: {largest} / {n_cells} cells ({:.1}%)",
        100.0 * largest as f64 / n_cells.max(1) as f64
    );
    eprintln!("  Final cluster sizes:");
    for (id, &count) in cluster_sizes.iter().enumerate() {
        eprintln!(
            "    {id:03}: {count} ({:.1}%)",
            100.0 * count as f64 / n_cells.max(1) as f64
        );
    }
    eprintln!("Starting multiplet mixture fitting...");
    let started = Instant::now();
    let multiplet_fits = crate::multiplet::detect(
        &raw,
        &q,
        &labels,
        config.multiplet_min_parent_fraction,
        config.multiplet_min_improvement,
    );
    let n_multiplets = multiplet_fits
        .iter()
        .filter(|f| !f.components.is_empty())
        .count();
    eprintln!(
        "  Multiplet explanations: {n_multiplets} clusters (full-gene linear-expression cluster-mean mixture model)"
    );
    finish_stage(&mut timings, "multiplet fitting", started);
    eprintln!("Starting UMAP...");
    let started = Instant::now();
    let viz = umap_2d(&pca_coords, config.umap_neighbors, config.umap_epochs)?;
    finish_stage(&mut timings, "UMAP", started);
    eprintln!("Starting Beacon integration...");
    let started = Instant::now();
    let beacon =
        if let Some(path) = source.filter(|p| p.extension().is_none_or(|ext| ext != "h5ad")) {
            load_beacon_blocks(&resolve_exonic(path)?, &norm.cells)?
        } else {
            Vec::new()
        };
    finish_stage(&mut timings, "Beacon integration", started);
    eprintln!("Starting ClonoMap integration...");
    let started = Instant::now();
    let clonomap =
        if let Some(path) = source.filter(|p| !p.extension().is_some_and(|ext| ext == "h5ad")) {
            load_clonomap(path, &norm.cells)?
        } else {
            None
        };
    finish_stage(&mut timings, "ClonoMap integration", started);
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
    let cell_multiplet_fits: Vec<_> = labels.iter().map(|&c| multiplet_fits[c].clone()).collect();
    write_all(
        &config.output_dir,
        &norm,
        &sampled,
        &q,
        &labels,
        &leiden_initial,
        &initial,
        &holdout,
        &pca_coords,
        &pca_var,
        &viz,
        &history,
        &beacon,
        clonomap.as_ref(),
        config.pseudo_samples,
        &cell_multiplet_fits,
        &selected,
    )?;
    crate::multiplet::write_clusters(
        &config.output_dir.join("cluster_multiplet_explanations.tsv"),
        &labels,
        &multiplet_fits,
    )?;
    write_cluster_quality(&config.output_dir.join("stats"), &pca_coords, &x, &labels)?;
    write_cluster_umi_stats(
        &config.output_dir.join("stats/cluster_umi_distribution.tsv"),
        &q,
        &labels,
    )?;
    finish_stage(&mut timings, "output writing", started);

    eprintln!("Starting plot generation...");
    let started = Instant::now();
    write_plots(
        &config.output_dir,
        &norm,
        &q,
        &labels,
        &leiden_initial,
        &pca_coords,
        &pca_var,
        &viz,
        &beacon,
    )?;
    crate::report::write_leiden_diagnostics(
        &config.output_dir,
        &pca_coords,
        &viz,
        &leiden_initial,
        &labels,
    )?;
    finish_stage(&mut timings, "plot generation", started);

    let mut annotation_columns = vec![
        (
            "multiplet_explanation",
            cell_multiplet_fits
                .iter()
                .map(|f| f.explanation.clone())
                .collect(),
        ),
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
    ];
    if let Some(c) = &clonomap {
        for (header, values) in c.headers.iter().zip(&c.columns) {
            annotation_columns.push((header.as_str(), values.clone()));
        }
    }
    let annotations = CellAnnotations::from_columns(annotation_columns, n_cells)?;
    let data = SingleCellData::new(norm.matrix, norm.features, norm.cells, annotations)?;
    Ok((
        data,
        AnalysisSummary {
            input_cells,
            filtered_cells: n_cells,
            input_features: n_genes,
            surviving_features: ng,
            variable_genes: selected.len(),
            initial_clusters: leiden_k,
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

fn write_cluster_umi_stats(path: &Path, q: &crate::qc::QcMetrics, labels: &[usize]) -> Result<()> {
    let n = labels.iter().copied().max().unwrap_or(0) + 1;
    let mut groups = vec![Vec::<usize>::new(); n];
    for (cell, &cluster) in labels.iter().enumerate() {
        groups[cluster].push(cell);
    }
    let mut out = String::from("cluster\tn_cells\tmetric\tmin\tq1\tmedian\tq3\tmax\n");
    for (cluster, cells) in groups.iter().enumerate() {
        for (name, values) in [("total_umi", &q.total), ("surviving_umi", &q.surviving)] {
            let mut sorted: Vec<f64> = cells.iter().map(|&i| values[i]).collect();
            sorted.sort_by(f64::total_cmp);
            if sorted.is_empty() {
                continue;
            }
            let quantile = |p: f64| -> f64 {
                let index = p * (sorted.len() - 1) as f64;
                let lo = index.floor() as usize;
                let hi = index.ceil() as usize;
                sorted[lo] + (sorted[hi] - sorted[lo]) * (index - lo as f64)
            };
            out.push_str(&format!(
                "cluster_{cluster:03}\t{}\t{name}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\n",
                cells.len(),
                quantile(0.0),
                quantile(0.25),
                quantile(0.5),
                quantile(0.75),
                quantile(1.0)
            ));
        }
    }
    std::fs::write(path, out)?;
    Ok(())
}
