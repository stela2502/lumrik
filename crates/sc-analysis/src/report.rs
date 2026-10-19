use crate::cluster::MergeStep;
use crate::data::SingleCellData;
use crate::norn::{BeaconBlock, ClonomapBlock};
use crate::qc::QcMetrics;
use crate::stats::{bh_adjust, mann_whitney};
use anyhow::Result;
use ndarray::Array2;
use plotters::prelude::*;
use plotters::style::TextStyle;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use rayon::prelude::*;
use sc_anndata::{AnnData, Column, Matrix, Table};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

pub(crate) fn write_all(
    out: &Path,
    norm: &SingleCellData,
    sampled: &SingleCellData,
    qc: &QcMetrics,
    labels: &[usize],
    initial: &[usize],
    spatial_raw: &[usize],
    holdout: &[bool],
    pca: &Array2<f32>,
    variance: &[f64],
    umap: &Array2<f32>,
    history: &[MergeStep],
    beacon: &[BeaconBlock],
    clonomap: Option<&ClonomapBlock>,
    pseudo_groups: usize,
    multiplet_fits: &[crate::multiplet::Fit],
) -> Result<()> {
    fs::create_dir_all(out)?;
    fs::create_dir_all(out.join("stats"))?;
    write_cells(
        out, norm, qc, labels, initial, holdout, umap, beacon, clonomap,
    )?;
    write_clusters(out, labels, holdout)?;
    write_beacon_summary(out, labels, beacon)?;
    write_qc_summary(out, qc)?;
    write_pca(out, norm, pca, variance)?;
    write_history(out, history)?;
    // Publish marker summaries before the potentially slow H5AD export.
    write_stats(
        &out.join("stats"),
        norm,
        labels,
        pseudo_groups,
        Some(initial),
        Some(out),
    )?;
    write_spatial_metadata(out, norm, labels, initial, spatial_raw)?;
    write_h5ad(
        out,
        norm,
        sampled,
        qc,
        labels,
        initial,
        spatial_raw,
        holdout,
        pca,
        umap,
        beacon,
        clonomap,
        multiplet_fits,
    )?;
    Ok(())
}

fn write_spatial_metadata(
    out: &Path,
    d: &SingleCellData,
    merged: &[usize],
    leiden: &[usize],
    raw: &[usize],
) -> Result<()> {
    anyhow::ensure!(
        [merged.len(), leiden.len(), raw.len()]
            .iter()
            .all(|&n| n == d.n_cells()),
        "spatial metadata cell count mismatch"
    );
    let mut w = BufWriter::new(File::create(out.join("spatial_cluster_crosswalk.tsv"))?);
    writeln!(w, "cell\tspatial_patch_raw\tleiden_initial\tleiden_merged")?;
    for i in 0..d.n_cells() {
        writeln!(
            w,
            "{}\tpatch_{:04}\tcluster_{:03}\tcluster_{:03}",
            d.cells[i], raw[i], leiden[i], merged[i]
        )?;
    }
    fs::write(
        out.join("SPATIAL_ANNOTATIONS.md"),
        concat!(
            "# Spatial and Leiden annotations\n\n",
            "`analysis.h5ad` stores cell-aligned `obs` columns:\n\n",
            "- `spatial_patch_raw`: original PCA-space geometric partition used for kNN acceleration.\n",
            "- `leiden_initial`: Leiden community before Pearson merging (also legacy `overcluster`).\n",
            "- `leiden_merged`: final community after Pearson merging (also legacy `cluster`).\n\n",
            "Patch IDs are arbitrary identifiers, not physical coordinates, ordered developmental stages, or biological clusters. ",
            "Patches are computed from PCA feature space, not the displayed 2D UMAP. ",
            "The same patch can contain cells from different Leiden communities; use both columns to subset.\n\n",
            "`spatial_cluster_crosswalk.tsv` duplicates these annotations with cell IDs for easy joins. ",
            "Use the H5AD expression matrix for per-cell gene detection and conditional intensity; ",
            "cluster membership alone cannot establish apoptosis, TP53 activity, or biological trajectories.\n"
        ),
    )?;
    Ok(())
}

fn write_h5ad(
    out: &Path,
    d: &SingleCellData,
    sampled: &SingleCellData,
    q: &QcMetrics,
    labels: &[usize],
    initial: &[usize],
    spatial_raw: &[usize],
    holdout: &[bool],
    pca: &Array2<f32>,
    umap: &Array2<f32>,
    beacon: &[BeaconBlock],
    clonomap: Option<&ClonomapBlock>,
    multiplet_fits: &[crate::multiplet::Fit],
) -> Result<()> {
    anyhow::ensure!(
        d.features == sampled.features && d.cells == sampled.cells,
        "H5AD layers have different cell/gene ordering"
    );
    anyhow::ensure!(
        multiplet_fits.len() == d.n_cells(),
        "multiplet annotation length mismatch"
    );
    anyhow::ensure!(
        labels.len() == d.n_cells()
            && initial.len() == d.n_cells()
            && spatial_raw.len() == d.n_cells(),
        "H5AD cluster annotations / cell count mismatch"
    );
    // sc-analysis stores feature x cell CSR. AnnData's invariant is
    // observation(cell) x variable(feature), so transpose exactly once here.
    let x = d.matrix.transpose_view().to_csr();
    let mut obs = Table::new("_index", d.cells.clone());
    // Preserve source sample/condition annotations, without overriding generated columns.
    let generated = [
        "total_umi",
        "detected_genes",
        "mitochondrial_umi",
        "ribosomal_umi",
        "vdj_umi",
        "surviving_umi",
        "overcluster",
        "cluster",
        "tech",
        "spatial_patch_raw",
        "leiden_initial",
        "leiden_merged",
    ];
    for name in d.annotations.headers() {
        if name == "cell" || generated.contains(&name.as_str()) {
            continue;
        }
        if let Some(values) = d.annotations.column(name) {
            obs.add_column(name.clone(), Column::String(values.to_vec()))?;
        }
    }
    obs.add_column("total_umi", Column::F64(q.total.clone()))?;
    obs.add_column(
        "detected_genes",
        Column::I64(q.detected.iter().map(|&x| x as i64).collect()),
    )?;
    obs.add_column("mitochondrial_umi", Column::F64(q.mito.clone()))?;
    obs.add_column("ribosomal_umi", Column::F64(q.ribo.clone()))?;
    obs.add_column("vdj_umi", Column::F64(q.vdj.clone()))?;
    obs.add_column("surviving_umi", Column::F64(q.surviving.clone()))?;
    obs.add_column(
        "overcluster",
        Column::I64(initial.iter().map(|&x| x as i64).collect()),
    )?;
    obs.add_column(
        "spatial_patch_raw",
        Column::String(
            spatial_raw
                .iter()
                .map(|&x| format!("patch_{x:04}"))
                .collect(),
        ),
    )?;
    obs.add_column(
        "leiden_initial",
        Column::String(initial.iter().map(|&x| format!("cluster_{x:03}")).collect()),
    )?;
    obs.add_column(
        "leiden_merged",
        Column::String(labels.iter().map(|&x| format!("cluster_{x:03}")).collect()),
    )?;
    obs.add_column(
        "cluster",
        Column::String(labels.iter().map(|x| format!("cluster_{x:03}")).collect()),
    )?;
    obs.add_column(
        "multiplet_explanation",
        Column::String(
            multiplet_fits
                .iter()
                .map(|f| f.explanation.clone())
                .collect(),
        ),
    )?;
    obs.add_column(
        "tech",
        Column::String(
            holdout
                .iter()
                .map(|&x| if x { "pseudo_deep" } else { "training" }.to_string())
                .collect(),
        ),
    )?;
    for b in beacon {
        obs.add_column(
            format!("{}_assignment", b.name),
            Column::String(b.assignment.clone()),
        )?;
        obs.add_column(
            format!("{}_called_features", b.name),
            Column::String(b.called_features.clone()),
        )?;
        obs.add_column(
            format!("{}_n_called", b.name),
            Column::String(b.n_called.clone()),
        )?;
        obs.add_column(
            format!("{}_best_feature", b.name),
            Column::String(b.best_feature.clone()),
        )?;
        obs.add_column(
            format!("{}_best_log_odds", b.name),
            Column::String(b.best_log_odds.clone()),
        )?;
    }
    if let Some(c) = clonomap {
        for (header, values) in c.headers.iter().zip(&c.columns) {
            obs.add_column(header.clone(), Column::String(values.clone()))?;
        }
    }
    let var = Table::new("_index", d.features.clone());
    let mut obsm = BTreeMap::new();
    obsm.insert("X_pca".to_string(), Matrix::DenseF32(pca.clone()));
    obsm.insert("X_umap".to_string(), Matrix::DenseF32(umap.clone()));
    let mut layers = BTreeMap::new();
    layers.insert(
        "resampled_umi".to_string(),
        Matrix::CsrF32(sampled.matrix.transpose_view().to_csr()),
    );
    AnnData {
        x: Matrix::CsrF32(x),
        obs,
        var,
        obsm,
        layers,
    }
    .write_h5ad(out.join("analysis.h5ad"))?;
    Ok(())
}

pub(crate) fn write_plots(
    out: &Path,
    norm: &SingleCellData,
    qc: &QcMetrics,
    labels: &[usize],
    initial: &[usize],
    pca: &Array2<f32>,
    variance: &[f64],
    umap: &Array2<f32>,
    beacon: &[BeaconBlock],
) -> Result<()> {
    anyhow::ensure!(
        pca.nrows() == umap.nrows() && initial.len() == umap.nrows(),
        "initial patch / PCA / UMAP cell ordering mismatch"
    );
    anyhow::ensure!(
        pca.ncols() >= 2 && umap.ncols() >= 2,
        "diagnostic plots require at least two coordinates"
    );
    let dir = out.join("plots");
    fs::create_dir_all(&dir)?;
    plots(&dir, norm, qc, labels, variance, umap, beacon)?;
    // Keep the traditional labeled cluster legend. Initial labels are Leiden,
    // not the spatial microclusters used only for neighbor indexing.
    scatter_cat(
        &dir.join("umap_leiden_initial.svg"),
        "UMAP — initial Leiden communities",
        umap,
        initial,
    )?;
    scatter_cat(
        &dir.join("umap_leiden_merged.svg"),
        "UMAP — merged Leiden communities",
        umap,
        labels,
    )?;
    if let Some(samples) = norm.annotations.column("sample") {
        scatter_labels(
            &dir.join("umap_samples.svg"),
            "UMAP — samples",
            umap,
            samples,
        )?;
    }
    Ok(())
}

fn frac(a: f64, b: f64) -> f64 {
    if b > 0.0 { a / b } else { 0.0 }
}
fn write_cells(
    out: &Path,
    d: &SingleCellData,
    q: &QcMetrics,
    l: &[usize],
    initial: &[usize],
    h: &[bool],
    u: &Array2<f32>,
    beacon: &[BeaconBlock],
    clonomap: Option<&ClonomapBlock>,
) -> Result<()> {
    let plasma = marker_score(d, &["Jchain", "Mzb1", "Sdc1", "Xbp1", "Prdm1"]);
    let b_cell = marker_score(d, &["Cd79a", "Cd79b", "Ms4a1", "Cd37", "H2-Aa"]);
    let mut w = BufWriter::new(File::create(out.join("cells.tsv"))?);
    write!(
        w,
        "cell\ttotal_umi\tdetected_genes\tmitochondrial_umi\tmitochondrial_fraction\tribosomal_umi\tribosomal_fraction\tvdj_umi\tvdj_fraction\tsurviving_umi\tplasma_score\tb_cell_score\tovercluster\tcluster\ttech\tumap_1\tumap_2"
    )?;
    if d.annotations.column("sample").is_some() {
        write!(w, "\tsample\tcellid")?;
    }
    for b in beacon {
        write!(
            w,
            "\t{}_assignment\t{}_called_features\t{}_n_called\t{}_best_feature\t{}_best_log_odds",
            b.name, b.name, b.name, b.name, b.name
        )?;
    }
    if let Some(c) = clonomap {
        for header in &c.headers {
            write!(w, "\t{header}")?;
        }
    }
    writeln!(w)?;
    for i in 0..d.n_cells() {
        write!(
            w,
            "{}\t{:.0}\t{}\t{:.0}\t{:.6}\t{:.0}\t{:.6}\t{:.0}\t{:.6}\t{:.0}\t{:.8}\t{:.8}\t{}\tcluster_{:03}\t{}\t{:.6}\t{:.6}",
            d.cells[i],
            q.total[i],
            q.detected[i],
            q.mito[i],
            frac(q.mito[i], q.total[i]),
            q.ribo[i],
            frac(q.ribo[i], q.total[i]),
            q.vdj[i],
            frac(q.vdj[i], q.total[i]),
            q.surviving[i],
            plasma[i],
            b_cell[i],
            initial[i],
            l[i],
            if h[i] { "pseudo_deep" } else { "training" },
            u[(i, 0)],
            u[(i, 1)]
        )?;
        if let (Some(samples), Some(cellids)) = (
            d.annotations.column("sample"),
            d.annotations.column("cellid"),
        ) {
            write!(w, "\t{}\t{}", samples[i], cellids[i])?;
        }
        for b in beacon {
            write!(
                w,
                "\t{}\t{}\t{}\t{}\t{}",
                b.assignment[i],
                b.called_features[i],
                b.n_called[i],
                b.best_feature[i],
                b.best_log_odds[i]
            )?;
        }
        if let Some(c) = clonomap {
            for column in &c.columns {
                write!(w, "\t{}", column[i])?;
            }
        }
        writeln!(w)?;
    }
    Ok(())
}

fn write_clusters(out: &Path, l: &[usize], h: &[bool]) -> Result<()> {
    let k = l.iter().copied().max().unwrap_or(0) + 1;
    let mut w = BufWriter::new(File::create(out.join("clusters.tsv"))?);
    writeln!(
        w,
        "cluster\tcells\tfraction_cells\ttraining_cells\tpseudo_deep_cells"
    )?;
    for cl in 0..k {
        let cells = (0..l.len()).filter(|&i| l[i] == cl).collect::<Vec<_>>();
        let held = cells.iter().filter(|&&i| h[i]).count();
        writeln!(
            w,
            "cluster_{cl:03}\t{}\t{:.8}\t{}\t{}",
            cells.len(),
            cells.len() as f64 / l.len().max(1) as f64,
            cells.len() - held,
            held
        )?;
    }
    Ok(())
}
fn write_beacon_summary(out: &Path, l: &[usize], beacon: &[BeaconBlock]) -> Result<()> {
    if beacon.is_empty() {
        return Ok(());
    }
    let mut w = BufWriter::new(File::create(out.join("beacon_summary.tsv"))?);
    writeln!(w, "feature_type\tcluster\tassignment\tcells")?;
    for b in beacon {
        let mut counts = std::collections::BTreeMap::<(usize, String), usize>::new();
        for (i, label) in b.label.iter().enumerate() {
            *counts.entry((l[i], label.clone())).or_default() += 1;
        }
        for ((cl, label), n) in counts {
            writeln!(w, "{}\tcluster_{:03}\t{}\t{}", b.name, cl, label, n)?;
        }
    }
    Ok(())
}
fn write_qc_summary(out: &Path, q: &QcMetrics) -> Result<()> {
    let n = q.total.len().max(1) as f64;
    let total = q.total.iter().sum::<f64>();
    let mito = q.mito.iter().sum::<f64>();
    let ribo = q.ribo.iter().sum::<f64>();
    let vdj = q.vdj.iter().sum::<f64>();
    let surviving = q.surviving.iter().sum::<f64>();
    let mut w = BufWriter::new(File::create(out.join("qc_summary.tsv"))?);
    writeln!(w, "metric\tvalue")?;
    writeln!(w, "cells\t{}", q.total.len())?;
    writeln!(w, "total_umi\t{total:.0}")?;
    writeln!(w, "mitochondrial_umi\t{mito:.0}")?;
    writeln!(w, "ribosomal_umi\t{ribo:.0}")?;
    writeln!(w, "vdj_umi\t{vdj:.0}")?;
    writeln!(w, "surviving_umi\t{surviving:.0}")?;
    writeln!(
        w,
        "mitochondrial_fraction\t{:.8}",
        if total > 0.0 { mito / total } else { 0.0 }
    )?;
    writeln!(
        w,
        "ribosomal_fraction\t{:.8}",
        if total > 0.0 { ribo / total } else { 0.0 }
    )?;
    writeln!(
        w,
        "vdj_fraction\t{:.8}",
        if total > 0.0 { vdj / total } else { 0.0 }
    )?;
    writeln!(w, "mean_total_umi_per_cell\t{:.4}", total / n)?;
    writeln!(w, "mean_surviving_umi_per_cell\t{:.4}", surviving / n)?;
    Ok(())
}
fn write_pca(out: &Path, d: &SingleCellData, p: &Array2<f32>, v: &[f64]) -> Result<()> {
    let mut w = BufWriter::new(File::create(out.join("pca.tsv"))?);
    write!(w, "cell")?;
    for j in 0..p.ncols() {
        write!(w, "\tPC{}", j + 1)?;
    }
    writeln!(w)?;
    for i in 0..p.nrows() {
        write!(w, "{}", d.cells[i])?;
        for j in 0..p.ncols() {
            write!(w, "\t{:.6}", p[(i, j)])?;
        }
        writeln!(w)?;
    }
    let mut x = BufWriter::new(File::create(out.join("pca_variance.tsv"))?);
    writeln!(
        x,
        "component\texplained_variance_fraction\tcumulative_fraction"
    )?;
    let mut c = 0.0;
    for (i, &z) in v.iter().enumerate() {
        c += z;
        writeln!(x, "{}\t{:.8}\t{:.8}", i + 1, z, c)?;
    }
    Ok(())
}
fn write_history(out: &Path, h: &[MergeStep]) -> Result<()> {
    let mut w = BufWriter::new(File::create(out.join("cluster_merge_history.tsv"))?);
    writeln!(
        w,
        "step\tfrom\tinto\tpearson\tfrom_cells\tinto_cells\tmerged_cells\tward_delta_sse\tpca_centroid_distance\tumap10_centroid_distance"
    )?;
    for x in h {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.8}\t{}\t{}\t{}\t{}\t{}\t{}",
            x.step,
            x.from,
            x.into,
            x.correlation,
            x.from_cells,
            x.into_cells,
            x.merged_cells,
            x.ward_delta_sse
                .map(|v| format!("{v:.8}"))
                .unwrap_or_default(),
            x.pca_centroid_distance
                .map(|v| format!("{v:.8}"))
                .unwrap_or_default(),
            x.umap10_centroid_distance
                .map(|v| format!("{v:.8}"))
                .unwrap_or_default()
        )?;
    }
    Ok(())
}
/// Mean-expression pseudo-bulks. Spatial seeds are never mixed across merged clusters.
/// When spatial labels are unavailable (external H5AD), use deterministic
/// round-robin grouping until PCA coordinates can be supplied explicitly.
fn pseudo_samples(
    d: &SingleCellData,
    cells: &[usize],
    n: usize,
    spatial: Option<&[usize]>,
) -> Vec<Vec<f64>> {
    use std::collections::BTreeMap;
    if cells.is_empty() {
        return Vec::new();
    }
    let groups: Vec<Vec<usize>> = if let Some(patches) = spatial {
        let mut by = BTreeMap::<usize, Vec<usize>>::new();
        for &cell in cells {
            by.entry(patches[cell]).or_default().push(cell);
        }
        by.into_values().collect()
    } else {
        let k = n.max(2).min(cells.len());
        let mut by = vec![Vec::new(); k];
        for (i, &cell) in cells.iter().enumerate() {
            by[i % k].push(cell);
        }
        by
    };
    let mut slot = vec![usize::MAX; d.n_cells()];
    for (i, cells) in groups.iter().enumerate() {
        for &cell in cells {
            slot[cell] = i;
        }
    }
    let mut out = vec![vec![0.0; d.n_features()]; groups.len()];
    for (g, row) in d.matrix.outer_iterator().enumerate() {
        for (c, value) in row.iter() {
            let i = slot[c];
            if i != usize::MAX {
                out[i][g] += *value as f64;
            }
        }
    }
    for (profile, cells) in out.iter_mut().zip(&groups) {
        let denom = cells.len() as f64;
        for value in profile {
            *value /= denom;
        }
    }
    out
}
/// Run the standard sc-analysis marker statistics for existing cell annotations,
/// without reclustering. The output schema is identical to the standard run.
pub fn write_group_stats(
    dir: &Path,
    data: &SingleCellData,
    grouping: &str,
    pseudo_groups: usize,
) -> Result<()> {
    let groups = data.annotations.groups(grouping)?;
    anyhow::ensure!(
        groups.len() >= 2,
        "marker statistics require at least two populations"
    );
    let mut labels = vec![usize::MAX; data.n_cells()];
    let mut names = Vec::new();
    for (idx, (name, cells)) in groups.iter().enumerate() {
        names.push(name.clone());
        for &cell in cells {
            labels[cell] = idx;
        }
    }
    anyhow::ensure!(
        labels.iter().all(|&x| x != usize::MAX),
        "missing cluster assignment"
    );
    fs::create_dir_all(dir)?;
    write_stats(dir, data, &labels, pseudo_groups, None, None)?;
    let mut mapping = BufWriter::new(File::create(dir.join("cluster_labels.tsv"))?);
    writeln!(mapping, "cluster\t{grouping}")?;
    for (idx, name) in names.iter().enumerate() {
        writeln!(mapping, "cluster_{idx:03}\t{name}")?;
    }
    Ok(())
}
// Gene, pseudobulk means, detection fractions, fold change, U, p, FDR,
// and observed mean expression conditional on detection (cluster/rest).
type MarkerRow = (usize, f64, f64, f64, f64, f64, f64, f64, f64, f64, f64);

fn general_score(r: &MarkerRow) -> f64 {
    let specificity = (r.3 - r.4).max(0.0);
    let intensity = ((1.0 + r.9) / (1.0 + r.10)).ln().max(0.0);
    specificity * (1.0 + intensity) * (1.0 + r.9).ln()
}

fn sparse_score(r: &MarkerRow) -> f64 {
    // Prefer strong positive cells, but require a nonzero observed fraction.
    // Softly reward enrichment and penalize ubiquitous genes.
    if r.3 <= 0.0 {
        return 0.0;
    }
    let intensity_gate = r.9 / (1.0 + r.9);
    (1.0 + r.9).ln() * intensity_gate * r.3.sqrt() * (1.0 - r.3).max(0.0) * (r.3 - r.4).max(0.0)
}

fn low_signal_score(r: &MarkerRow) -> f64 {
    // Prefer detectable changes in sparse/weak genes, without an arbitrary
    // detection cutoff. This is a descriptive prioritization, NOT a p-value.
    if r.3 == 0.0 && r.4 == 0.0 {
        return 0.0;
    }
    // A strong signal on either side belongs primarily in the strong-expression
    // ranking. Cubic attenuation prevents detection differences from dominating
    // this ranking for genes with substantial positive-cell expression.
    let weak_weight = 1.0 / (1.0 + r.9.max(r.10)).powi(3);
    let detection_change = (r.3 - r.4).abs();
    let intensity_change = ((1.0 + r.9) / (1.0 + r.10)).ln().abs();
    weak_weight * (detection_change + 0.1 * intensity_change * r.3.min(r.4).sqrt())
}

const STATS_README: &str = r#"# sc-analysis marker statistics

This folder contains cluster-versus-rest comparisons. Clusters are Leiden
communities after the configured merging step (or annotation groups when
using group statistics). The same cell may not belong to two clusters.

## Files

- `cluster_NNN.tsv`: **original inferential ranking**, sorted by adjusted
  Mann–Whitney p-value, with pseudobulk mean expression, detection fractions,
  fold change, U statistic, raw p-value, and Benjamini–Hochberg FDR.
- `cluster_NNN_general.tsv`: descriptive **broad marker** ranking.
- `cluster_NNN_sparse_strong.tsv`: descriptive **minority/high-expression**
  ranking. A rare high-count gene is not discarded simply for being rare.
- `cluster_NNN_low_expression.tsv`: descriptive **low-expression change**
  ranking, useful for screening weakly expressed regulators including TFs.
  **This is NOT a TF-only list**: no TF annotation database is applied.
- `marker_summary_general.tsv`, `marker_summary_sparse_strong.tsv`, and
  `marker_summary_low_expression.tsv`: one row per cluster, with top 50
  up/down gene names and aligned scores (comma-separated). These summaries
  use the exact same ranked rows as their corresponding cluster tables.
  A positive/negative direction is determined by the legacy log2 fold change;
  the low-expression score itself is two-sided.
- `marker_summary.tsv` and `top20_up_down.tsv`: legacy significance-ranked
  summaries. Despite the name, the latter currently includes up to 50 per
  direction; these are NOT the new descriptive rankings.
- `cluster_labels.tsv` (annotation-group analysis only): cluster-to-group map.

## Expression measurements

`pct_cluster` / `pct_rest` are fractions of cells with observed expression > 0.
`positive_mean_cluster` / `positive_mean_rest` are the arithmetic means of
observed positive values only, zero if no cells are positive.
`mean_cluster` / `mean_rest` in the legacy table are means of pseudobulk
samples, not conditional positive-cell means. Their exact units depend on the
input normalization. The descriptive files retain both measurements.

The two distinct sparse patterns are **rare/strong expression** (a minority
of cells with high observed counts) and **rare/weak expression** (few cells
with low observed counts). Neither by itself proves a distinct cell state.

## Descriptive scores (NOT statistical tests)

All scores are exploratory, dimensionless and depend on the input expression
scale. They are not probabilities, significance levels or cross-dataset
comparable effect sizes. Rankings are cluster-specific.

Let d = pct_cluster, b = pct_rest, x = positive_mean_cluster,
y = positive_mean_rest. Natural logarithms are used.

- General: `max(d-b,0) * (1+max(ln((1+x)/(1+y)),0)) * ln(1+x)`.
- Sparse/strong: `ln(1+x) * [x/(1+x)] * sqrt(d) * max(1-d,0) * max(d-b,0)`.
  The additional intensity gate suppresses weak positive-cell signals.
- Low-expression: `[abs(d-b) + 0.1*abs(ln((1+x)/(1+y)))*sqrt(min(d,b))]
  / [1 + max(x,y)]^3`. The stronger attenuation prevents high-expression
  genes from ranking well solely because their detection rates differ.

The intensity weights are continuous rather than hard classification cutoffs.
A gene can still appear in both lists, but the rankings now have a stronger
preference for different positive-cell expression regimes.

These formulas are initial prioritization heuristics, not validated marker
models. Scores of zero may arise from no expression or no positive enrichment.
For low-expression candidates, inspect the detection and positive-mean columns
and original p-values together; a change in detection can reflect sequencing
depth rather than regulation. The low-expression ranking is two-sided.

## Interpretation and limitations

Sparse high expression may reflect a minority phenotype, doublets or ambient
RNA. Sparse low expression may reflect weak transcription or sampling/dropout.
Depth, batch, condition, and library complexity can confound comparisons.
Strong changes in erythroid, TP53 or apoptosis genes do not establish DBA
rescue, differentiation arrest or cell death without additional evidence.

For CRISPR experiments, compare candidate markers with guide assignments and
replicates before making causal claims. TF identity requires external
annotation; this output does not infer TF status or regulatory targets.

The descriptive scores are not used to recompute the legacy p-values or FDR.
Do not interpret their rank as inferential evidence.
"#;

pub(crate) fn write_stats(
    dir: &Path,
    d: &SingleCellData,
    l: &[usize],
    groups: usize,
    spatial: Option<&[usize]>,
    summary_root: Option<&Path>,
) -> Result<()> {
    if let Some(patches) = spatial {
        anyhow::ensure!(
            patches.len() == d.n_cells(),
            "PCA patch/cell count mismatch"
        );
    }
    fs::create_dir_all(dir)?;
    fs::write(dir.join("README.md"), STATS_README)?;
    let k = l.iter().copied().max().unwrap_or(0) + 1;
    let mut top = BufWriter::new(File::create(dir.join("top20_up_down.tsv"))?);
    let mut summary = BufWriter::new(File::create(dir.join("marker_summary.tsv"))?);
    let mut descriptive_summaries = ["general", "sparse_strong", "low_expression"]
        .map(|name| {
            File::create(dir.join(format!("marker_summary_{name}.tsv"))).map(BufWriter::new)
        })
        .into_iter()
        .collect::<std::io::Result<Vec<_>>>()?;
    for writer in &mut descriptive_summaries {
        writeln!(
            writer,
            "cluster\tn_cells\tfraction_cells\ttop50_up\ttop50_up_scores\ttop50_down\ttop50_down_scores"
        )?;
    }
    // Create the discoverable root-level summary immediately, not after all stats.
    let mut root_summary = summary_root
        .map(|root| File::create(root.join("marker_summary.tsv")))
        .transpose()?;
    writeln!(
        top,
        "cluster\tdirection\trank\tgene\tmean_cluster\tmean_rest\tpct_cluster\tpct_rest\tlog2_fold_change\twilcoxon_u\tp_value\tp_adj"
    )?;
    writeln!(
        summary,
        "cluster\tn_cells\tfraction_cells\ttop20_up\ttop20_down"
    )?;
    summary.flush()?;
    if let Some(root) = root_summary.as_mut() {
        writeln!(
            root,
            "cluster\tn_cells\tfraction_cells\ttop20_up\ttop20_down"
        )?;
        root.flush()?;
    }
    for cl in 0..k {
        let inside = (0..d.n_cells()).filter(|&i| l[i] == cl).collect::<Vec<_>>();
        let outside = (0..d.n_cells()).filter(|&i| l[i] != cl).collect::<Vec<_>>();
        if inside.is_empty() || outside.is_empty() {
            continue;
        }
        eprintln!(
            "  marker stats: cluster {}/{} ({} cells)",
            cl + 1,
            k,
            inside.len()
        );
        let a = pseudo_samples(d, &inside, groups, spatial);
        let b = pseudo_samples(d, &outside, groups, spatial);
        // One sparse row traversal per gene; avoid O(cells) random CSR lookups.
        let mut membership = vec![false; d.n_cells()];
        for &cell in &inside {
            membership[cell] = true;
        }
        let raw: Vec<_> = (0..d.n_features())
            .into_par_iter()
            .map(|g| {
                let av = a.iter().map(|x| x[g]).collect::<Vec<_>>();
                let bv = b.iter().map(|x| x[g]).collect::<Vec<_>>();
                let ma = av.iter().sum::<f64>() / av.len() as f64;
                let mb = bv.iter().sum::<f64>() / bv.len() as f64;
                let (u, p) = mann_whitney(&av, &bv);
                let mut in_detected = 0usize;
                let mut out_detected = 0usize;
                let mut in_sum = 0.0_f64;
                let mut out_sum = 0.0_f64;
                if let Some(row) = d.matrix.outer_view(g) {
                    for (cell, value) in row.iter() {
                        if *value > 0.0 {
                            if membership[cell] {
                                in_detected += 1;
                                in_sum += *value as f64;
                            } else {
                                out_detected += 1;
                                out_sum += *value as f64;
                            }
                        }
                    }
                }
                let da = in_detected as f64 / inside.len() as f64;
                let db = out_detected as f64 / outside.len() as f64;
                let fc = ((ma + 1e-6) / (mb + 1e-6)).log2();
                let pos_in = if in_detected > 0 {
                    in_sum / in_detected as f64
                } else {
                    0.0
                };
                let pos_out = if out_detected > 0 {
                    out_sum / out_detected as f64
                } else {
                    0.0
                };
                (g, ma, mb, da, db, fc, u, p, pos_in, pos_out)
            })
            .collect();
        let adj = bh_adjust(&raw.iter().map(|x| x.7).collect::<Vec<_>>());
        let mut rows = raw
            .into_iter()
            .enumerate()
            .map(|(i, r)| (r.0, r.1, r.2, r.3, r.4, r.5, r.6, r.7, adj[i], r.8, r.9))
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| {
            a.8.total_cmp(&b.8)
                .then_with(|| a.7.total_cmp(&b.7))
                .then_with(|| b.5.abs().total_cmp(&a.5.abs()))
                .then_with(|| d.features[a.0].cmp(&d.features[b.0]))
        });
        // Orthogonal descriptive rankings. These do not replace the inferential table.
        let mut general = rows.clone();
        general.sort_by(|a, b| {
            general_score(b)
                .total_cmp(&general_score(a))
                .then_with(|| a.8.total_cmp(&b.8))
                .then_with(|| d.features[a.0].cmp(&d.features[b.0]))
        });
        let mut sparse = rows.clone();
        sparse.sort_by(|a, b| {
            sparse_score(b)
                .total_cmp(&sparse_score(a))
                .then_with(|| d.features[a.0].cmp(&d.features[b.0]))
        });
        let mut low = rows.clone();
        low.sort_by(|a, b| {
            low_signal_score(b)
                .total_cmp(&low_signal_score(a))
                .then_with(|| d.features[a.0].cmp(&d.features[b.0]))
        });
        for (summary_writer, (name, ranked, score)) in descriptive_summaries.iter_mut().zip([
            ("general", &general, general_score as fn(&MarkerRow) -> f64),
            (
                "sparse_strong",
                &sparse,
                sparse_score as fn(&MarkerRow) -> f64,
            ),
            (
                "low_expression",
                &low,
                low_signal_score as fn(&MarkerRow) -> f64,
            ),
        ]) {
            let select = |up: bool| -> (String, String) {
                let selected = ranked
                    .iter()
                    .filter(|r| score(r) > 0.0 && if up { r.5 > 0.0 } else { r.5 < 0.0 })
                    .take(50)
                    .collect::<Vec<_>>();
                (
                    selected
                        .iter()
                        .map(|r| d.features[r.0].as_str())
                        .collect::<Vec<_>>()
                        .join(","),
                    selected
                        .iter()
                        .map(|r| format!("{:.8}", score(r)))
                        .collect::<Vec<_>>()
                        .join(","),
                )
            };
            let (up_genes, up_scores) = select(true);
            let (down_genes, down_scores) = select(false);
            writeln!(
                summary_writer,
                "cluster_{cl:03}\t{}\t{:.8}\t{up_genes}\t{up_scores}\t{down_genes}\t{down_scores}",
                inside.len(),
                inside.len() as f64 / d.n_cells().max(1) as f64
            )?;
            let mut file = BufWriter::new(File::create(
                dir.join(format!("cluster_{cl:03}_{name}.tsv")),
            )?);
            writeln!(
                file,
                "rank\tgene\tscore\tpct_cluster\tpct_rest\tpositive_mean_cluster\tpositive_mean_rest\tmean_cluster\tmean_rest\tlog2_fold_change\tp_adj"
            )?;
            for (rank, r) in ranked.iter().enumerate() {
                writeln!(
                    file,
                    "{}\t{}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.6e}",
                    rank + 1,
                    d.features[r.0],
                    score(r),
                    r.3,
                    r.4,
                    r.9,
                    r.10,
                    r.1,
                    r.2,
                    r.5,
                    r.8
                )?;
            }
        }
        let mut w = BufWriter::new(File::create(dir.join(format!("cluster_{cl:03}.tsv")))?);
        writeln!(
            w,
            "rank\tgene\tmean_cluster\tmean_rest\tpct_cluster\tpct_rest\tlog2_fold_change\twilcoxon_u\tp_value\tp_adj"
        )?;
        for (rank, r) in rows.iter().enumerate() {
            writeln!(
                w,
                "{}\t{}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.4}\t{:.6e}\t{:.6e}",
                rank + 1,
                d.features[r.0],
                r.1,
                r.2,
                r.3,
                r.4,
                r.5,
                r.6,
                r.7,
                r.8
            )?;
        }
        let up = rows
            .iter()
            .filter(|r| r.5 > 0.0)
            .take(50)
            .collect::<Vec<_>>();
        let down = rows
            .iter()
            .filter(|r| r.5 < 0.0)
            .take(50)
            .collect::<Vec<_>>();
        for (direction, set) in [("up", &up), ("down", &down)] {
            for (rank, r) in set.iter().enumerate() {
                writeln!(
                    top,
                    "cluster_{cl:03}\t{direction}\t{}\t{}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.8}\t{:.4}\t{:.6e}\t{:.6e}",
                    rank + 1,
                    d.features[r.0],
                    r.1,
                    r.2,
                    r.3,
                    r.4,
                    r.5,
                    r.6,
                    r.7,
                    r.8
                )?;
            }
        }
        let up_names = up
            .iter()
            .map(|r| d.features[r.0].as_str())
            .collect::<Vec<_>>()
            .join(",");
        let down_names = down
            .iter()
            .map(|r| d.features[r.0].as_str())
            .collect::<Vec<_>>()
            .join(",");
        let summary_line = format!(
            "cluster_{cl:03}\t{}\t{:.8}\t{up_names}\t{down_names}",
            inside.len(),
            inside.len() as f64 / d.n_cells().max(1) as f64
        );
        writeln!(summary, "{summary_line}")?;
        summary.flush()?;
        if let Some(root) = root_summary.as_mut() {
            writeln!(root, "{summary_line}")?;
            root.flush()?;
        }
    }
    Ok(())
}
fn marker_score(d: &SingleCellData, genes: &[&str]) -> Vec<f64> {
    let mut score = vec![0.0; d.n_cells()];
    for gene in genes {
        let Some(g) = d
            .features
            .iter()
            .position(|feature| feature.eq_ignore_ascii_case(gene))
        else {
            continue;
        };
        if let Some(row) = d.matrix.outer_view(g) {
            for (cell, value) in row.iter() {
                score[cell] += *value as f64;
            }
        }
    }
    score
}

fn plots(
    dir: &Path,
    d: &SingleCellData,
    q: &QcMetrics,
    l: &[usize],
    var: &[f64],
    u: &Array2<f32>,
    beacon: &[BeaconBlock],
) -> Result<()> {
    bar_plot(&dir.join("pca_variance.svg"), "PCA explained variance", var)?;
    let plasma = marker_score(d, &["Jchain", "Mzb1", "Sdc1", "Xbp1", "Prdm1"]);
    let b_cell = marker_score(d, &["Cd79a", "Cd79b", "Ms4a1", "Cd37", "H2-Aa"]);
    scatter_num(
        &dir.join("umap_plasma_score.svg"),
        "UMAP — plasma / plasmablast score",
        u,
        &plasma,
    )?;
    scatter_num(
        &dir.join("umap_b_cell_score.svg"),
        "UMAP — B-cell score",
        u,
        &b_cell,
    )?;
    scatter_cat(
        &dir.join("umap_clusters.svg"),
        "UMAP — final clusters",
        u,
        l,
    )?;
    scatter_num(
        &dir.join("umap_total_umi.svg"),
        "UMAP — total UMI",
        u,
        &q.total,
    )?;
    scatter_num(
        &dir.join("umap_surviving_umi.svg"),
        "UMAP — surviving UMI",
        u,
        &q.surviving,
    )?;
    let mito = (0..q.total.len())
        .map(|i| frac(q.mito[i], q.total[i]))
        .collect::<Vec<_>>();
    let ribo = (0..q.total.len())
        .map(|i| frac(q.ribo[i], q.total[i]))
        .collect::<Vec<_>>();
    scatter_num(
        &dir.join("umap_mitochondrial.svg"),
        "UMAP — mitochondrial fraction",
        u,
        &mito,
    )?;
    scatter_num(
        &dir.join("umap_ribosomal.svg"),
        "UMAP — ribosomal fraction",
        u,
        &ribo,
    )?;
    let vdj = (0..q.total.len())
        .map(|i| frac(q.vdj[i], q.total[i]))
        .collect::<Vec<_>>();
    scatter_num(&dir.join("umap_vdj.svg"), "UMAP — V/D/J fraction", u, &vdj)?;
    hist(
        &dir.join("qc_mitochondrial.svg"),
        "Mitochondrial fraction",
        &mito,
    )?;
    hist(&dir.join("qc_ribosomal.svg"), "Ribosomal fraction", &ribo)?;
    hist(&dir.join("qc_umi.svg"), "Surviving UMI", &q.surviving)?;
    umi_rank_plot(&dir.join("qc_umi_rank.svg"), &q.total)?;
    for b in beacon {
        scatter_labels(
            &dir.join(format!("umap_{}.svg", b.name)),
            &format!("UMAP — {} Beacon assignment", b.name),
            u,
            &b.label,
        )?;
    }
    Ok(())
}
fn range(u: &Array2<f32>, j: usize) -> (f32, f32) {
    let mut lo = f32::INFINITY;
    let mut hi = f32::NEG_INFINITY;
    for i in 0..u.nrows() {
        lo = lo.min(u[(i, j)]);
        hi = hi.max(u[(i, j)]);
    }
    if lo == hi {
        (lo - 1.0, hi + 1.0)
    } else {
        (lo, hi)
    }
}

// Same ordered-colour-pair / four-axis encoding as ClonoMap's rooted tree.
// A larger palette is safe here because spatial diagnostics are not constrained
// to ClonoMap's six high-contrast lineage colours.
const PATCH_PALETTE: [RGBColor; 12] = [
    RGBColor(0, 0, 0),
    RGBColor(211, 211, 211),
    RGBColor(230, 159, 0),
    RGBColor(0, 114, 178),
    RGBColor(213, 94, 0),
    RGBColor(204, 121, 167),
    RGBColor(0, 158, 115),
    RGBColor(240, 228, 66),
    RGBColor(86, 180, 233),
    RGBColor(128, 64, 0),
    RGBColor(120, 60, 180),
    RGBColor(80, 130, 40),
];

// Distinct, high-contrast colour pairs come FIRST; only then do we rotate.
// Pair black with light colours and light grey with dark ones.
// Interleave families so adjacent IDs do not all share the same base colour.
fn patch_symbol(id: usize) -> (RGBColor, RGBColor, usize) {
    const LIGHT: [usize; 5] = [2, 5, 7, 8, 11];
    const DARK: [usize; 5] = [3, 4, 6, 9, 10];
    let mut pairs = Vec::new();
    for i in 0..LIGHT.len() {
        pairs.push((PATCH_PALETTE[0], PATCH_PALETTE[LIGHT[i]]));
        pairs.push((PATCH_PALETTE[1], PATCH_PALETTE[DARK[i]]));
    }
    pairs.push((PATCH_PALETTE[0], PATCH_PALETTE[1]));
    // Other ordered pairs follow, with the same luminance-aware contrast rule.
    fn luminance(c: RGBColor) -> f64 {
        let channel = |x: u8| {
            let v = x as f64 / 255.0;
            if v <= 0.04045 {
                v / 12.92
            } else {
                ((v + 0.055) / 1.055).powf(2.4)
            }
        };
        0.2126 * channel(c.0) + 0.7152 * channel(c.1) + 0.0722 * channel(c.2)
    }
    let mut others = Vec::new();
    for a in 0..PATCH_PALETTE.len() {
        for b in (a + 1)..PATCH_PALETTE.len() {
            let ca = PATCH_PALETTE[a];
            let cb = PATCH_PALETTE[b];
            if pairs
                .iter()
                .any(|&(x, y)| (x == ca && y == cb) || (x == cb && y == ca))
            {
                continue;
            }
            let (lo, hi) = if luminance(ca) < luminance(cb) {
                (ca, cb)
            } else {
                (cb, ca)
            };
            let contrast = (luminance(hi) + 0.05) / (luminance(lo) + 0.05);
            if contrast >= 3.0 {
                others.push((lo, hi));
            }
        }
    }
    others.sort_by(|&(a, b), &(c, d)| {
        let contrast = |x: RGBColor, y: RGBColor| (luminance(y) + 0.05) / (luminance(x) + 0.05);
        contrast(c, d).total_cmp(&contrast(a, b))
    });
    pairs.extend(others);
    let pair_count = pairs.len();
    let pair_index = id % pair_count;
    let cycle = id / pair_count;
    let (a, b) = pairs[pair_index];
    let axis = (cycle / 2) % 4;
    if cycle % 2 == 0 {
        (a, b, axis)
    } else {
        (b, a, axis)
    }
}

fn split_half(radius: i32, axis: usize, second: bool) -> Vec<(i32, i32)> {
    let tilt = match axis {
        0 => 90.0_f64,
        1 => 45.0,
        2 => 0.0,
        _ => -45.0,
    }
    .to_radians();
    let offset = if second { std::f64::consts::PI } else { 0.0 };
    let mut points = vec![(0, 0)];
    for step in 0..=18 {
        let a = tilt + offset + std::f64::consts::PI * step as f64 / 18.0;
        points.push((
            (radius as f64 * a.cos()).round() as i32,
            (radius as f64 * a.sin()).round() as i32,
        ));
    }
    points
}

// Keep colour assignments deterministic across initial and merged Leiden plots.
// The primary cue is the ordered colour pair. Rotate the split only when
// another cluster has a potentially confusable pair (including reversed pairs).
fn leiden_symbol(id: usize) -> (RGBColor, RGBColor, usize) {
    fn distance(a: RGBColor, b: RGBColor) -> f64 {
        let dr = (a.0 as f64 - b.0 as f64) / 255.0;
        let dg = (a.1 as f64 - b.1 as f64) / 255.0;
        let db = (a.2 as f64 - b.2 as f64) / 255.0;
        (dr * dr + dg * dg + db * db).sqrt()
    }
    fn similar(a: (RGBColor, RGBColor), b: (RGBColor, RGBColor)) -> bool {
        let forward = distance(a.0, b.0) + distance(a.1, b.1);
        let reversed = distance(a.0, b.1) + distance(a.1, b.0);
        forward.min(reversed) < 0.72
    }
    let mut assigned: Vec<(RGBColor, RGBColor, usize)> = Vec::with_capacity(id + 1);
    for current in 0..=id {
        let (a, b, _) = patch_symbol(current);
        let mut occupied = [false; 4];
        for &(pa, pb, axis) in &assigned {
            if similar((a, b), (pa, pb)) {
                occupied[axis] = true;
            }
        }
        let axis = occupied
            .iter()
            .position(|&taken| !taken)
            .unwrap_or(current % 4);
        assigned.push((a, b, axis));
    }
    assigned[id]
}

fn scatter_cat(path: &Path, title: &str, u: &Array2<f32>, v: &[usize]) -> Result<()> {
    let mut counts = BTreeMap::<usize, usize>::new();
    for &cl in v {
        *counts.entry(cl).or_default() += 1;
    }
    // Lexical (numeric ID) legend order is stable across runs and plots.
    // Every entry shows the exact same two-colour, oriented symbol as the cells.
    // Reserve extra canvas width for the legend; never shrink the 1000x900 plot.
    // Overflow legends wrap into additional columns rather than increasing plot height.
    const LEGEND_ROWS: usize = 27;
    const LEGEND_COLUMN_WIDTH: u32 = 290;
    let legend_columns = counts.len().div_ceil(LEGEND_ROWS).max(1) as u32;
    let root = SVGBackend::new(path, (1000 + LEGEND_COLUMN_WIDTH * legend_columns, 900))
        .into_drawing_area();
    root.fill(&WHITE)?;
    let (plot, legend) = root.split_horizontally(1000);
    let (x0, x1) = range(u, 0);
    let (y0, y1) = range(u, 1);
    let mut c = ChartBuilder::on(&plot)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .build_cartesian_2d(x0..x1, y0..y1)?;

    let symbols: BTreeMap<usize, _> = counts.keys().map(|&id| (id, leiden_symbol(id))).collect();
    c.draw_series((0..u.nrows()).map(|i| {
        let &(a, b, axis) = &symbols[&v[i]];
        EmptyElement::at((u[(i, 0)], u[(i, 1)]))
            + Polygon::new(split_half(4, axis, false), a.filled())
            + Polygon::new(split_half(4, axis, true), b.filled())
    }))?;

    for (row, (&cl, &count)) in counts.iter().enumerate() {
        let (a, b, axis) = symbols[&cl];
        let x = 20 + (row / LEGEND_ROWS) as i32 * LEGEND_COLUMN_WIDTH as i32;
        let y = 70 + (row % LEGEND_ROWS) as i32 * 30;
        legend.draw(&Polygon::new(
            split_half(9, axis, false)
                .into_iter()
                .map(|(dx, dy)| (dx + x, dy + y))
                .collect::<Vec<_>>(),
            a.filled(),
        ))?;
        legend.draw(&Polygon::new(
            split_half(9, axis, true)
                .into_iter()
                .map(|(dx, dy)| (dx + x, dy + y))
                .collect::<Vec<_>>(),
            b.filled(),
        ))?;
        legend.draw(&Circle::new((x, y), 9, BLACK.stroke_width(1)))?;
        legend.draw(&Text::new(
            format!(
                "cluster_{cl:03} — {count} ({:.1}%)",
                100.0 * count as f64 / v.len().max(1) as f64
            ),
            (x + 18, y),
            TextStyle::from(("sans-serif", 16).into_font()).pos(Pos::new(HPos::Left, VPos::Center)),
        ))?;
    }
    root.present()?;
    Ok(())
}

fn scatter_labels(path: &Path, title: &str, u: &Array2<f32>, labels: &[String]) -> Result<()> {
    let mut names = labels.to_vec();
    names.sort();
    names.dedup();
    let root = SVGBackend::new(path, (1250, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    let (x0, x1) = range(u, 0);
    let (y0, y1) = range(u, 1);
    let mut c = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .right_y_label_area_size(180)
        .build_cartesian_2d(x0..x1, y0..y1)?;
    for (name_idx, name) in names.iter().enumerate() {
        let color = Palette99::pick(name_idx).mix(1.0);
        c.draw_series(
            (0..u.nrows())
                .filter(|&i| labels[i] == *name)
                .map(|i| Circle::new((u[(i, 0)], u[(i, 1)]), 2, color.filled())),
        )?
        .label(name.clone())
        .legend(move |(x, y)| Circle::new((x, y), 4, color.filled()));
    }
    c.configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.8))
        .draw()?;
    root.present()?;
    Ok(())
}
fn yellow_red_blue(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let lerp =
        |a: u8, b: u8, u: f64| -> u8 { (a as f64 + (b as f64 - a as f64) * u).round() as u8 };

    if t <= 0.5 {
        let u = t * 2.0;
        RGBColor(255, lerp(220, 0, u), 0)
    } else {
        let u = (t - 0.5) * 2.0;
        RGBColor(lerp(255, 0, u), 0, lerp(0, 178, u))
    }
}

fn legend_value(value: f64, integer_values: bool) -> String {
    if integer_values {
        format!("{value:.0}")
    } else {
        format!("{value:.3}")
    }
}

fn scatter_num(path: &Path, title: &str, u: &Array2<f32>, v: &[f64]) -> Result<()> {
    let root = SVGBackend::new(path, (1250, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    let (plot, legend) = root.split_horizontally(1080);
    let (x0, x1) = range(u, 0);
    let (y0, y1) = range(u, 1);
    let lo = v.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let integer_values = v.iter().all(|value| value.fract().abs() < 1e-9);
    let mut c = ChartBuilder::on(&plot)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .build_cartesian_2d(x0..x1, y0..y1)?;
    c.draw_series((0..u.nrows()).map(|i| {
        let t = if hi > lo {
            (v[i] - lo) / (hi - lo)
        } else {
            0.0
        };
        Circle::new((u[(i, 0)], u[(i, 1)]), 2, yellow_red_blue(t).filled())
    }))?;

    let bar_x0 = 35;
    let bar_x1 = 75;
    let bar_y0 = 140;
    let bar_y1 = 700;
    let steps = 100;
    for step in 0..steps {
        let t0 = step as f64 / steps as f64;
        let y_top = bar_y1 - ((step + 1) * (bar_y1 - bar_y0) / steps) as i32;
        let y_bottom = bar_y1 - (step * (bar_y1 - bar_y0) / steps) as i32;
        legend.draw(&Rectangle::new(
            [(bar_x0, y_top), (bar_x1, y_bottom)],
            yellow_red_blue(t0).filled(),
        ))?;
    }
    legend.draw(&Text::new(
        legend_value(hi, integer_values),
        (85, bar_y0 + 5),
        ("sans-serif", 16).into_font(),
    ))?;
    legend.draw(&Text::new(
        legend_value(lo, integer_values),
        (85, bar_y1),
        ("sans-serif", 16).into_font(),
    ))?;
    legend.draw(&Text::new(
        "high",
        (35, bar_y0 - 20),
        ("sans-serif", 16).into_font(),
    ))?;
    legend.draw(&Text::new(
        "low",
        (35, bar_y1 + 25),
        ("sans-serif", 16).into_font(),
    ))?;
    root.present()?;
    Ok(())
}

fn umi_rank_plot(path: &Path, umi: &[f64]) -> Result<()> {
    let mut ranked = umi.to_vec();
    ranked.sort_by(|a, b| b.total_cmp(a));
    let points = ranked
        .iter()
        .enumerate()
        .filter(|(_, umi)| **umi > 0.0)
        .map(|(i, umi)| (((i + 1) as f64).log10(), umi.log10()))
        .collect::<Vec<_>>();
    if points.is_empty() {
        return Ok(());
    }
    let x_max = points.last().map(|x| x.0).unwrap_or(1.0).max(1.0);
    let y_max = points.iter().map(|x| x.1).fold(0.0_f64, f64::max).max(1.0);
    let root = SVGBackend::new(path, (1000, 700)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Barcode rank — total UMI", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0_f64..x_max, 0.0_f64..y_max)?;
    chart
        .configure_mesh()
        .x_desc("log10 cell rank")
        .y_desc("log10 total UMI")
        .draw()?;
    chart.draw_series(LineSeries::new(points, &BLUE))?;
    root.present()?;
    Ok(())
}

fn hist(path: &Path, title: &str, v: &[f64]) -> Result<()> {
    let lo = v.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = v.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let bins = 50usize;
    let step = ((hi - lo) / bins as f64).max(1e-9);
    let mut n = vec![0u32; bins];
    for &x in v {
        let b = (((x - lo) / step) as usize).min(bins - 1);
        n[b] += 1;
    }
    let root = SVGBackend::new(path, (1000, 700)).into_drawing_area();
    root.fill(&WHITE)?;
    let ymax = *n.iter().max().unwrap_or(&1);
    let mut c = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .build_cartesian_2d(0usize..bins, 0u32..ymax)?;
    c.configure_mesh().draw()?;
    c.draw_series(
        n.iter()
            .enumerate()
            .map(|(i, &y)| Rectangle::new([(i, 0), (i + 1, y)], BLUE.filled())),
    )?;
    root.present()?;
    Ok(())
}
fn bar_plot(path: &Path, title: &str, v: &[f64]) -> Result<()> {
    let root = SVGBackend::new(path, (1000, 700)).into_drawing_area();
    root.fill(&WHITE)?;
    let ymax = v.iter().copied().fold(0.0, f64::max).max(1e-9);
    let mut c = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .build_cartesian_2d(0usize..v.len(), 0.0..ymax)?;
    c.configure_mesh().draw()?;
    c.draw_series(
        v.iter()
            .enumerate()
            .map(|(i, &y)| Rectangle::new([(i, 0.0), (i + 1, y)], BLUE.filled())),
    )?;
    root.present()?;
    Ok(())
}

/// Minimal PCA/UMAP Leiden diagnostics with the original labeled legend.
pub(crate) fn write_leiden_diagnostics(
    out: &Path,
    pca: &Array2<f32>,
    umap: &Array2<f32>,
    initial: &[usize],
    merged: &[usize],
) -> Result<()> {
    let dir = out.join("plots");
    fs::create_dir_all(&dir)?;
    scatter_cat(
        &dir.join("pca_leiden_initial.svg"),
        "PCA — initial Leiden communities",
        pca,
        initial,
    )?;
    scatter_cat(
        &dir.join("umap_direct_leiden_initial.svg"),
        "UMAP — initial Leiden communities",
        umap,
        initial,
    )?;
    scatter_cat(
        &dir.join("umap_direct_leiden_merged.svg"),
        "UMAP — merged Leiden communities",
        umap,
        merged,
    )?;
    Ok(())
}
