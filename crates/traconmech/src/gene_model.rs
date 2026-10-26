//! Portable cluster centroids. Apply never updates a trained model.
use anyhow::{Context, Result, ensure};
use sc_analysis::SingleCellData;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Serialize, Deserialize)]
struct Model {
    format: String,
    version: u32,
    gene: String,
    grouping: String,
    features: Vec<String>,
    clusters: Vec<Centroid>,
    #[serde(default)]
    modules: BTreeMap<String, usize>,
}
#[derive(Serialize, Deserialize)]
struct Centroid {
    name: String,
    n_cells: usize,
    log1p_mean: Vec<f64>,
}
pub fn create(
    data: &SingleCellData,
    cluster: &str,
    gene: &str,
    model_path: &Path,
    out: &Path,
) -> Result<()> {
    std::fs::create_dir_all(out)?;
    if let Some(parent) = model_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    ensure!(
        !data.cells.is_empty() && !data.features.is_empty(),
        "empty training dataset"
    );
    let groups = data.annotations.groups(cluster)?;
    ensure!(!groups.is_empty(), "no clusters in {}", cluster);
    // Each centroid is independent. Keep the per-gene accumulation order unchanged.
    let clusters: Vec<Centroid> = groups.into_par_iter().filter_map(|(name, indices)| {
        if name.trim().is_empty() { return None; }
        let mut selected = vec![false; data.cells.len()];
        for &i in &indices { selected[i] = true; }
        let mut means = vec![0.0; data.features.len()];
        for (g, row) in data.matrix.outer_iterator().enumerate() {
            for (c, v) in row.iter() {
                if selected[c] { means[g] += f64::from(*v); }
            }
        }
        for m in &mut means { *m /= indices.len() as f64; }
        Some(Centroid { name, n_cells: indices.len(), log1p_mean: means })
    }).collect();
    ensure!(!clusters.is_empty(), "no nonempty cluster labels");
    let modules = read_training_modules(&out.join("module_membership.tsv"))?;
    ensure!(
        !modules.is_empty(),
        "training produced no effect modules; use --effect-cluster during create"
    );
    let model = Model {
        modules,
        format: "traconmech.gene-model.centroids".into(),
        version: 1,
        gene: gene.to_owned(),
        grouping: cluster.to_owned(),
        features: data.features.clone(),
        clusters,
    };
    let file =
        File::create(model_path).with_context(|| format!("creating {}", model_path.display()))?;
    serde_json::to_writer(BufWriter::new(file), &model)?;
    println!(
        "Saved {} clusters and {} genes to {}",
        model.clusters.len(),
        model.features.len(),
        model_path.display()
    );
    let mut report = BufWriter::new(File::create(out.join("gene_model_training_clusters.tsv"))?);
    writeln!(report, "cluster\tn_cells\tgenes")?;
    for c in &model.clusters {
        writeln!(
            report,
            "{}\t{}\t{}",
            c.name,
            c.n_cells,
            model.features.len()
        )?;
    }
    Ok(())
}
pub fn validate(model_path: &Path, gene: &str, data: &SingleCellData) -> Result<()> {
    let model = read_model(model_path)?;
    ensure!(
        model.gene == gene,
        "model gene {} does not match requested {gene}",
        model.gene
    );
    ensure!(!data.cells.is_empty(), "empty external dataset");
    Ok(())
}
fn read_model(path: &Path) -> Result<Model> {
    let model: Model = serde_json::from_reader(
        File::open(path).with_context(|| format!("opening model {}", path.display()))?,
    )?;
    ensure!(
        model.format == "traconmech.gene-model.centroids" && model.version == 1,
        "unsupported gene model format/version"
    );
    ensure!(
        !model.clusters.is_empty() && !model.features.is_empty(),
        "empty model"
    );
    for c in &model.clusters {
        ensure!(
            c.log1p_mean.len() == model.features.len(),
            "invalid centroid length"
        );
    }
    Ok(model)
}

fn read_training_modules(path: &Path) -> Result<BTreeMap<String, usize>> {
    use std::io::{BufRead, BufReader};
    let mut lines = BufReader::new(
        File::open(path).with_context(|| format!("reading training modules {}", path.display()))?,
    )
    .lines();
    let header = lines.next().context("empty module membership table")??;
    let columns: Vec<_> = header.split('\t').collect();
    ensure!(
        columns.first() == Some(&"gene"),
        "invalid module membership header"
    );
    let ids: Vec<usize> = columns
        .iter()
        .skip(1)
        .map(|s| {
            s.strip_prefix("module_")
                .context("invalid module column")?
                .parse::<usize>()
                .map_err(Into::into)
        })
        .collect::<Result<_>>()?;
    let mut modules = BTreeMap::new();
    for line in lines {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        ensure!(
            fields.len() == ids.len() + 1,
            "invalid module membership row: {line}"
        );
        let assigned: Vec<_> = fields
            .iter()
            .skip(1)
            .zip(&ids)
            .filter_map(|(v, id)| (*v == "1").then_some(*id))
            .collect();
        ensure!(
            assigned.len() == 1,
            "gene {} must belong to exactly one module",
            fields[0]
        );
        ensure!(
            modules.insert(fields[0].to_owned(), assigned[0]).is_none(),
            "duplicate module gene {}",
            fields[0]
        );
    }
    Ok(modules)
}

/// Plot the external expression of the ORIGINAL training modules; never refit.
pub fn plot_frozen_modules(
    data: &SingleCellData,
    model_path: &Path,
    cluster: &str,
    condition: &str,
    cases: &str,
    controls: &str,
    out: &Path,
    label_angle: Option<f64>,
    count_label_layout: &str,
) -> Result<()> {
    let model = read_model(model_path)?;
    ensure!(
        !model.modules.is_empty(),
        "model has no frozen modules; recreate with --effect-cluster"
    );
    let clusters = data
        .annotations
        .column(cluster)
        .with_context(|| format!("missing cluster annotation {cluster}"))?;
    let conditions = data
        .annotations
        .column(condition)
        .with_context(|| format!("missing condition annotation {condition}"))?;
    let case_set: std::collections::HashSet<&str> = cases.split(',').map(str::trim).filter(|s| !s.is_empty()).collect();
    let control_set: std::collections::HashSet<&str> = controls.split(',').map(str::trim).filter(|s| !s.is_empty()).collect();
    let mut groups: BTreeMap<String, (Vec<usize>, Vec<usize>)> = BTreeMap::new();
    for (i, (cl, co)) in clusters.iter().zip(conditions).enumerate() {
        let pair = groups.entry(cl.clone()).or_default();
        if case_set.contains(co.as_str()) {
            pair.0.push(i);
        }
        if control_set.contains(co.as_str()) {
            pair.1.push(i);
        }
    }
    groups.retain(|_, (a, b)| !a.is_empty() || !b.is_empty());
    ensure!(!groups.is_empty(), "no case/control cells to plot");
    let mut ordered: Vec<String> = groups.keys().cloned().collect();
    if let (Some(names), Some(orders)) = (
        data.annotations.column(cluster),
        data.annotations.column("developmental_order"),
    ) {
        let mut rank: HashMap<&str, f64> = HashMap::new();
        for (name, order) in names.iter().zip(orders) {
            if let Ok(v) = order.parse::<f64>() {
                if v.is_finite() {
                    rank.entry(name.as_str())
                        .and_modify(|old| *old = old.min(v))
                        .or_insert(v);
                }
            }
        }
        if !rank.is_empty() {
            ordered.sort_by(|a, b| {
                rank.get(a.as_str())
                    .copied()
                    .unwrap_or(f64::INFINITY)
                    .total_cmp(&rank.get(b.as_str()).copied().unwrap_or(f64::INFINITY))
                    .then_with(|| a.cmp(b))
            });
        }
    } else if let Some(levels) = data.annotations.categorical_levels(cluster) {
        let rank: HashMap<&str, usize> = levels
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        ordered.sort_by_key(|s| rank.get(s.as_str()).copied().unwrap_or(usize::MAX));
    }
    let references: BTreeMap<String, String> =
        ordered.iter().map(|s| (s.clone(), s.clone())).collect();
    crate::module_plots::write_module_plots(
        out,
        data,
        &model.modules,
        &ordered,
        &groups,
        &groups,
        &references,
        cases,
        controls,
        label_angle,
        count_label_layout,
    )?;
    let mut report = BufWriter::new(File::create(out.join("frozen_module_coverage.tsv"))?);
    writeln!(report, "module\ttraining_genes\texternal_genes")?;
    let mut totals: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
    for (gene, id) in &model.modules {
        let pair = totals.entry(*id).or_default();
        pair.0 += 1;
        if data.feature_index(gene).is_some() {
            pair.1 += 1;
        }
    }
    for (id, (total, present)) in totals {
        writeln!(report, "{id}\t{total}\t{present}")?;
    }
    Ok(())
}

pub fn apply(data: &SingleCellData, model_path: &Path, condition: &str, out: &Path) -> Result<()> {
    let model = read_model(model_path)?;
    let labels = data
        .annotations
        .column(condition)
        .with_context(|| format!("missing condition annotation {condition}"))?;
    ensure!(
        labels.len() == data.cells.len(),
        "condition annotation length mismatch"
    );
    let lookup: HashMap<&str, usize> = data
        .features
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let shared: Vec<(usize, usize)> = model
        .features
        .iter()
        .enumerate()
        .filter_map(|(i, s)| lookup.get(s.as_str()).map(|&j| (i, j)))
        .collect();
    ensure!(
        shared.len() >= 10,
        "only {} overlapping genes; need at least 10",
        shared.len()
    );
    let mut seen = std::collections::HashSet::new();
    ensure!(
        model.features.iter().all(|g| seen.insert(g)),
        "duplicate genes in model"
    );
    let coverage = shared.len() as f64 / model.features.len() as f64;
    ensure!(
        coverage >= 0.5,
        "only {:.1}% of model genes overlap external data (minimum 50%)",
        coverage * 100.0
    );
    std::fs::create_dir_all(out)?;
    let quantified = crate::module_score::score(data, &model.modules)?;
    let module_ids = &quantified.ids;
    let module_columns: Vec<String> = module_ids.iter().map(|id| format!("module_{id}")).collect();
    let mut output = BufWriter::new(File::create(out.join("gene_model_scores.tsv"))?);
    write!(output, "cell")?;
    for name in &module_columns { write!(output, "\t{name}")?; }
    writeln!(output, "\tbest_module\tcondition")?;
    let mut summaries: BTreeMap<(String, usize), (usize, f64, f64)> = BTreeMap::new();
    for (cell_index, cell) in data.cells.iter().enumerate() {
        write!(output, "{cell}")?;
        let mut best: Option<(usize, f64)> = None;
        for (j, &id) in module_ids.iter().enumerate() {
            let score = quantified.scores[cell_index][j];
            write!(output, "\t{score:.6}")?;
            if best.map_or(true, |(_, previous)| score > previous) {
                best = Some((j, score));
            }
            let entry = summaries.entry((labels[cell_index].clone(), id)).or_default();
            entry.0 += 1;
            entry.1 += score;
            entry.2 += score * score;
        }
        writeln!(output, "\t{}\t{}", best.map(|(j, _)| module_columns[j].as_str()).unwrap_or("NA"), labels[cell_index])?;
    }
    let mut raw = BufWriter::new(File::create(out.join("gene_model_module_means.tsv"))?);
    write!(raw, "cell")?;
    for name in &module_columns { write!(raw, "\t{name}")?; }
    writeln!(raw)?;
    for (cell, values) in data.cells.iter().zip(&quantified.module_means) {
        write!(raw, "{cell}")?;
        for value in values { write!(raw, "\t{value:.6}")?; }
        writeln!(raw)?;
    }
    let mut summary = BufWriter::new(File::create(out.join("gene_model_condition_scores.tsv"))?);
    writeln!(summary, "condition_column\tcondition\tmodule\tn_cells\tmean_module_score\tsd_module_score")?;
    for ((label, id), (n, sum, sum_sq)) in summaries {
        let mean = sum / n as f64;
        let sd = (sum_sq / n as f64 - mean * mean).max(0.0).sqrt();
        writeln!(summary, "{condition}\t{label}\tmodule_{id}\t{n}\t{mean:.6}\t{sd:.6}")?;
    }
    let mut module_coverage = BufWriter::new(File::create(out.join("gene_model_module_coverage.tsv"))?);
    writeln!(module_coverage, "module\tmodel_genes\tobserved_genes")?;
    for (id, (total, present)) in quantified.coverage {
        writeln!(module_coverage, "module_{id}\t{total}\t{present}")?;
    }
    let mut info = BTreeMap::new();
    info.insert("model_genes", model.features.len());
    info.insert("overlapping_genes", shared.len());
    info.insert("cells", data.cells.len());
    serde_json::to_writer_pretty(File::create(out.join("gene_model_coverage.json"))?, &info)?;
    println!("Quantified {} frozen modules in {} cells ({:.1}% model gene overlap): {}", module_ids.len(), data.cells.len(), coverage * 100.0, out.display());
    Ok(())
}
