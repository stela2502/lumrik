use crate::cluster::MergeStep;
use crate::data::SingleCellData;
use crate::norn::BeaconBlock;
use crate::qc::QcMetrics;
use crate::stats::{bh_adjust, mann_whitney};
use anyhow::Result;
use flate2::{Compression, write::GzEncoder};
use ndarray::Array2;
use plotters::prelude::*;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

pub(crate) fn write_all(
    out: &Path,
    norm: &SingleCellData,
    qc: &QcMetrics,
    labels: &[usize],
    initial: &[usize],
    holdout: &[bool],
    pca: &Array2<f32>,
    variance: &[f64],
    umap: &Array2<f32>,
    history: &[MergeStep],
    beacon: &[BeaconBlock],
    pseudo_groups: usize,
) -> Result<()> {
    fs::create_dir_all(out)?;
    fs::create_dir_all(out.join("normalized"))?;
    fs::create_dir_all(out.join("stats"))?;
    write_mex(&out.join("normalized"), norm)?;
    write_cells(out, norm, qc, labels, initial, holdout, umap, beacon)?;
    write_clusters(out, labels, holdout)?;
    write_beacon_summary(out, labels, beacon)?;
    write_qc_summary(out, qc)?;
    write_pca(out, norm, pca, variance)?;
    write_history(out, history)?;
    write_stats(&out.join("stats"), norm, labels, pseudo_groups)?;
    Ok(())
}

pub(crate) fn write_plots(
    out: &Path,
    norm: &SingleCellData,
    qc: &QcMetrics,
    labels: &[usize],
    variance: &[f64],
    umap: &Array2<f32>,
    beacon: &[BeaconBlock],
) -> Result<()> {
    let dir = out.join("plots");
    fs::create_dir_all(&dir)?;
    plots(&dir, norm, qc, labels, variance, umap, beacon)
}
fn gz(path: &Path) -> Result<GzEncoder<File>> {
    Ok(GzEncoder::new(File::create(path)?, Compression::default()))
}
fn write_mex(dir: &Path, d: &SingleCellData) -> Result<()> {
    let mut m = gz(&dir.join("matrix.mtx.gz"))?;
    writeln!(m, "%%MatrixMarket matrix coordinate real general")?;
    writeln!(m, "{} {} {}", d.n_features(), d.n_cells(), d.matrix.nnz())?;
    let csc = d.matrix.to_csc();
    for (c, col) in csc.outer_iterator().enumerate() {
        for (r, v) in col.iter() {
            writeln!(m, "{} {} {:.8}", r + 1, c + 1, v)?;
        }
    }
    m.finish()?;
    let mut f = gz(&dir.join("features.tsv.gz"))?;
    for x in &d.features {
        writeln!(f, "{0}\t{0}\tGene Expression", x)?;
    }
    f.finish()?;
    let mut b = gz(&dir.join("barcodes.tsv.gz"))?;
    for x in &d.cells {
        writeln!(b, "{x}")?;
    }
    b.finish()?;
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
) -> Result<()> {
    let plasma = marker_score(d, &["Jchain", "Mzb1", "Sdc1", "Xbp1", "Prdm1"]);
    let b_cell = marker_score(d, &["Cd79a", "Cd79b", "Ms4a1", "Cd37", "H2-Aa"]);
    let mut w = BufWriter::new(File::create(out.join("cells.tsv"))?);
    write!(
        w,
        "cell\ttotal_umi\tdetected_genes\tmitochondrial_umi\tmitochondrial_fraction\tribosomal_umi\tribosomal_fraction\tvdj_umi\tvdj_fraction\tsurviving_umi\tplasma_score\tb_cell_score\tovercluster\tcluster\ttech\tumap_1\tumap_2"
    )?;
    for b in beacon {
        write!(
            w,
            "\t{}_assignment\t{}_called_features\t{}_n_called\t{}_best_feature\t{}_best_log_odds",
            b.name, b.name, b.name, b.name, b.name
        )?;
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
        "step\tfrom\tinto\tpearson\tfrom_cells\tinto_cells\tmerged_cells"
    )?;
    for x in h {
        writeln!(
            w,
            "{}\t{}\t{}\t{:.8}\t{}\t{}\t{}",
            x.step, x.from, x.into, x.correlation, x.from_cells, x.into_cells, x.merged_cells
        )?;
    }
    Ok(())
}
fn pseudo_samples(d: &SingleCellData, cells: &[usize], n: usize) -> Vec<Vec<f64>> {
    let k = n.max(2).min(cells.len().max(1));
    let mut out = vec![vec![0.0; d.n_features()]; k];
    let mut cnt = vec![0usize; k];
    let mut slot = vec![usize::MAX; d.n_cells()];
    for (i, &c) in cells.iter().enumerate() {
        slot[c] = i % k;
        cnt[i % k] += 1;
    }
    for (g, row) in d.matrix.outer_iterator().enumerate() {
        for (c, v) in row.iter() {
            let s = slot[c];
            if s != usize::MAX {
                out[s][g] += *v as f64;
            }
        }
    }
    for s in 0..k {
        let z = cnt[s].max(1) as f64;
        for x in &mut out[s] {
            *x /= z;
        }
    }
    out
}
fn write_stats(dir: &Path, d: &SingleCellData, l: &[usize], groups: usize) -> Result<()> {
    let k = l.iter().copied().max().unwrap_or(0) + 1;
    let mut top = BufWriter::new(File::create(dir.join("top20_up_down.tsv"))?);
    let mut summary = BufWriter::new(File::create(dir.join("marker_summary.tsv"))?);
    writeln!(
        top,
        "cluster\tdirection\trank\tgene\tmean_cluster\tmean_rest\tpct_cluster\tpct_rest\tlog2_fold_change\twilcoxon_u\tp_value\tp_adj"
    )?;
    writeln!(
        summary,
        "cluster\tn_cells\tfraction_cells\ttop20_up\ttop20_down"
    )?;
    for cl in 0..k {
        let inside = (0..d.n_cells()).filter(|&i| l[i] == cl).collect::<Vec<_>>();
        let outside = (0..d.n_cells()).filter(|&i| l[i] != cl).collect::<Vec<_>>();
        if inside.is_empty() || outside.is_empty() {
            continue;
        }
        let a = pseudo_samples(d, &inside, groups);
        let b = pseudo_samples(d, &outside, groups);
        let mut raw = Vec::with_capacity(d.n_features());
        for g in 0..d.n_features() {
            let av = a.iter().map(|x| x[g]).collect::<Vec<_>>();
            let bv = b.iter().map(|x| x[g]).collect::<Vec<_>>();
            let ma = av.iter().sum::<f64>() / av.len() as f64;
            let mb = bv.iter().sum::<f64>() / bv.len() as f64;
            let (u, p) = mann_whitney(&av, &bv);
            let da = inside
                .iter()
                .filter(|&&c| d.matrix.get(g, c).copied().unwrap_or(0.0) > 0.0)
                .count() as f64
                / inside.len() as f64;
            let db = outside
                .iter()
                .filter(|&&c| d.matrix.get(g, c).copied().unwrap_or(0.0) > 0.0)
                .count() as f64
                / outside.len() as f64;
            let fc = ((ma + 1e-6) / (mb + 1e-6)).log2();
            raw.push((g, ma, mb, da, db, fc, u, p));
        }
        let adj = bh_adjust(&raw.iter().map(|x| x.7).collect::<Vec<_>>());
        let mut rows = raw
            .into_iter()
            .enumerate()
            .map(|(i, r)| (r.0, r.1, r.2, r.3, r.4, r.5, r.6, r.7, adj[i]))
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| {
            a.8.total_cmp(&b.8)
                .then_with(|| b.5.abs().total_cmp(&a.5.abs()))
                .then_with(|| d.features[a.0].cmp(&d.features[b.0]))
        });
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
            .take(20)
            .collect::<Vec<_>>();
        let down = rows
            .iter()
            .filter(|r| r.5 < 0.0)
            .take(20)
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
        writeln!(
            summary,
            "cluster_{cl:03}\t{}\t{:.8}\t{up_names}\t{down_names}",
            inside.len(),
            inside.len() as f64 / d.n_cells().max(1) as f64
        )?;
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
fn scatter_cat(path: &Path, title: &str, u: &Array2<f32>, v: &[usize]) -> Result<()> {
    let k = v.iter().copied().max().unwrap_or(0) + 1;
    let mut counts = vec![0usize; k];
    for &cl in v {
        counts[cl] += 1;
    }

    let root = SVGBackend::new(path, (1250, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    let (plot, legend) = root.split_horizontally(1000);
    let (x0, x1) = range(u, 0);
    let (y0, y1) = range(u, 1);
    let mut c = ChartBuilder::on(&plot)
        .caption(title, ("sans-serif", 30))
        .margin(20)
        .build_cartesian_2d(x0..x1, y0..y1)?;

    for cl in 0..k {
        let color = Palette99::pick(cl).mix(1.0);
        c.draw_series(
            (0..u.nrows())
                .filter(|&i| v[i] == cl)
                .map(|i| Circle::new((u[(i, 0)], u[(i, 1)]), 2, color.filled())),
        )?;

        let y = 70 + cl as i32 * 30;
        legend.draw(&Circle::new((20, y), 5, color.filled()))?;
        legend.draw(&Text::new(
            format!(
                "cluster_{cl:03} — {} ({:.1}%)",
                counts[cl],
                100.0 * counts[cl] as f64 / v.len().max(1) as f64
            ),
            (35, y + 5),
            ("sans-serif", 16).into_font(),
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
