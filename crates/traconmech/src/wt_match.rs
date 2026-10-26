//! Select an adequate within-Louvain WT reference, otherwise find the closest
//! other Louvain population and compare against its WT cells only.
use crate::experiment::require_annotation;
use anyhow::{Context, Result, bail};
use sc_analysis::SingleCellData;
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

fn selected(s: &str) -> HashSet<&str> {
    s.split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect()
}
pub(crate) fn profile(data: &SingleCellData, cells: &[usize]) -> Vec<f64> {
    let mut mask = vec![false; data.cells.len()];
    for &i in cells {
        mask[i] = true;
    }
    data.matrix
        .outer_iterator()
        .map(|row| {
            let sum: f64 = row
                .iter()
                .filter(|(i, _)| mask[*i])
                .map(|(_, v)| f64::from(*v))
                .sum();
            sum / cells.len() as f64
        })
        .collect()
}
pub(crate) fn correlation(a: &[f64], b: &[f64]) -> Option<f64> {
    let a: Vec<_> = a.iter().map(|x| x.max(0.).ln_1p()).collect();
    let b: Vec<_> = b.iter().map(|x| x.max(0.).ln_1p()).collect();
    // Pearson on transformed expression, with means recomputed after transform.
    pearson_raw(&a, &b)
}
fn pearson_raw(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.len() < 2 {
        return None;
    }
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let (mut ab, mut aa, mut bb) = (0., 0., 0.);
    for (&x, &y) in a.iter().zip(b) {
        let (x, y) = (x - ma, y - mb);
        ab += x * y;
        aa += x * x;
        bb += y * y;
    }
    (aa > 0. && bb > 0.).then(|| ab / (aa * bb).sqrt())
}
pub fn run(
    data: &SingleCellData,
    out: &Path,
    cluster_col: &str,
    condition_col: &str,
    case: &str,
    control: &str,
    target: &str,
    min_wt: usize,
) -> Result<()> {
    let _ = min_wt; // Legacy CLI parameter; no minimum for observed WT cells.
    let clusters = require_annotation(data, cluster_col, "population")?;
    let conditions = require_annotation(data, condition_col, "condition")?;
    let cases = selected(case);
    let controls = selected(control);
    if cases.is_empty() || controls.is_empty() || !cases.is_disjoint(&controls) {
        bail!("--case and --control must be nonempty and disjoint");
    }
    let mut groups: BTreeMap<String, (Vec<usize>, Vec<usize>)> = BTreeMap::new();
    for (i, (cluster, condition)) in clusters.iter().zip(conditions.iter()).enumerate() {
        if cases.contains(condition.as_str()) {
            groups.entry(cluster.clone()).or_default().0.push(i);
        } else if controls.contains(condition.as_str()) {
            groups.entry(cluster.clone()).or_default().1.push(i);
        }
    }
    let (target_case, target_control) = groups
        .get(target)
        .with_context(|| format!("target cluster {target} absent from selected samples"))?;
    if target_case.is_empty() {
        bail!("no case cells in target cluster {target}");
    }
    let case_profile = profile(data, target_case);
    let mut candidates = Vec::new();
    for (cluster, (case_cells, control_cells)) in &groups {
        let mut pooled = case_cells.clone();
        pooled.extend(control_cells);
        if pooled.is_empty() {
            continue;
        }
        let similarity = correlation(&case_profile, &profile(data, &pooled));
        candidates.push((
            cluster.clone(),
            case_cells.len(),
            control_cells.len(),
            similarity,
        ));
    }
    candidates.sort_by(|a, b| {
        b.3.unwrap_or(f64::NEG_INFINITY)
            .total_cmp(&a.3.unwrap_or(f64::NEG_INFINITY))
            .then_with(|| a.0.cmp(&b.0))
    });
    // Only observed WT cells from the requested cluster may be used.
    if target_control.is_empty() {
        bail!("target cluster {target} has no WT cells; cannot compute a within-cluster effect");
    }
    let reference = target.to_owned();
    let (_, reference_wt) = &groups[&reference];
    let wt_profile = profile(data, reference_wt);
    fs::create_dir_all(out)?;
    let mut ranks = BufWriter::new(File::create(out.join("closest_wt_clusters.tsv"))?);
    writeln!(
        ranks,
        "target_cluster\tcandidate_cluster\tcase_cells\twt_cells\tpearson_log1p_mean\teligible\tselected"
    )?;
    for (cluster, case_n, wt_n, score) in &candidates {
        writeln!(
            ranks,
            "{target}\t{cluster}\t{case_n}\t{wt_n}\t{}\t{}\t{}",
            score.map(|x| x.to_string()).unwrap_or_else(|| "NA".into()),
            cluster == &reference,
            cluster == &reference
        )?;
    }
    let mut diff = BufWriter::new(File::create(out.join("dba_vs_closest_wt.tsv"))?);
    writeln!(
        diff,
        "gene\ttarget_cluster\twt_cluster\tcase_n\twt_n\tcase_mean\twt_mean\tlog2fc\tcase_detected\twt_detected"
    )?;
    let mut case_mask = vec![false; data.cells.len()];
    let mut wt_mask = vec![false; data.cells.len()];
    for &i in target_case {
        case_mask[i] = true;
    }
    for &i in reference_wt {
        wt_mask[i] = true;
    }
    for (g, row) in data.matrix.outer_iterator().enumerate() {
        let mut cd = 0usize;
        let mut wd = 0usize;
        for (i, value) in row.iter() {
            if *value != 0. {
                if case_mask[i] {
                    cd += 1;
                }
                if wt_mask[i] {
                    wd += 1;
                }
            }
        }
        let cm = case_profile[g];
        let wm = wt_profile[g];
        let fc = ((cm + 0.01) / (wm + 0.01)).log2();
        writeln!(
            diff,
            "{}\t{target}\t{reference}\t{}\t{}\t{cm}\t{wm}\t{fc}\t{}\t{}",
            data.features[g],
            target_case.len(),
            reference_wt.len(),
            cd as f64 / target_case.len() as f64,
            wd as f64 / reference_wt.len() as f64
        )?;
    }
    eprintln!(
        "[traconmech] DBA L{target}: {} cells; WT L{reference}: {} cells ({}); {} genes",
        target_case.len(),
        reference_wt.len(),
        "within-cluster",
        data.features.len()
    );
    Ok(())
}
