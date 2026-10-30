//! SVG trajectories of mean original log1p expression for frozen module genes.
//! Clusters are categorical unless the user supplied an explicit biological order.
use anyhow::Result;
use statrs::distribution::{ContinuousCDF, StudentsT};
use sc_analysis::SingleCellData;
use std::collections::BTreeMap;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
fn cell_stats(v: &[f64]) -> (f64, f64, f64, f64, f64) {
    if v.is_empty() { return (f64::NAN, f64::NAN, f64::NAN, f64::NAN, f64::NAN); }
    let n = v.len() as f64;
    let mean = v.iter().sum::<f64>() / n;
    if v.len() < 2 { return (mean, f64::NAN, f64::NAN, f64::NAN, f64::NAN); }
    let sd = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
    let se = sd / n.sqrt();
    (mean, sd, se, mean - 1.96 * se, mean + 1.96 * se)
}
/// Two-sided Welch test on independent cell-level module means.
/// Returns MUT-WT, t statistic, degrees of freedom, and unadjusted p.
fn welch_test(wt: &[f64], mutant: &[f64]) -> Option<(f64, f64, f64, f64)> {
    if wt.len() < 2 || mutant.len() < 2 { return None; }
    let (mw, sw, _, _, _) = cell_stats(wt);
    let (mm, sm, _, _, _) = cell_stats(mutant);
    let a = sw * sw / wt.len() as f64;
    let b = sm * sm / mutant.len() as f64;
    let variance = a + b;
    if !variance.is_finite() || variance <= 0.0 { return None; }
    let denominator = a * a / (wt.len() - 1) as f64 + b * b / (mutant.len() - 1) as f64;
    if denominator <= 0.0 { return None; }
    let df = variance * variance / denominator;
    if !df.is_finite() || df <= 0.0 { return None; }
    let t = (mm - mw) / variance.sqrt();
    let distribution = StudentsT::new(0.0, 1.0, df).ok()?;
    let p = (2.0 * distribution.sf(t.abs())).clamp(0.0, 1.0);
    Some((mm - mw, t, df, p))
}
fn p_label(p: Option<f64>) -> String {
    match p {
        Some(p) => {
            let stars = if p < 0.001 { "** " } else if p < 0.01 { "* " } else { "" };
            format!("{stars}p={p:.2e}")
        }
        None => "p=NA".to_owned(),
    }
}
fn fmt(v: f64) -> String { if v.is_finite() { v.to_string() } else { "NA".to_owned() } }
fn polyline(points: &[(f64, f64)]) -> String {
    points
        .iter()
        .map(|(x, y)| format!("{x:.1},{y:.1}"))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Plot arithmetic means of the supplied log1p expression, with no rescaling.
pub fn write_module_plots(
    out: &Path,
    data: &SingleCellData,
    modules: &BTreeMap<String, usize>,
    ordered: &[String],
    groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    actual_groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    reference_clusters: &BTreeMap<String, String>,
    case_label: &str,
    control_label: &str,
    label_angle: Option<f64>,
    count_label_layout: &str,
) -> Result<()> {
    write_module_plots_named(
        out,
        data,
        modules,
        ordered,
        groups,
        actual_groups,
        reference_clusters,
        case_label,
        control_label,
        "",
        label_angle,
        count_label_layout,
    )
}

pub fn write_module_plots_named(
    out: &Path,
    data: &SingleCellData,
    modules: &BTreeMap<String, usize>,
    ordered: &[String],
    groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    actual_groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    reference_clusters: &BTreeMap<String, String>,
    case_label: &str,
    control_label: &str,
    suffix: &str,
    label_angle: Option<f64>,
    count_label_layout: &str,
) -> Result<()> {
    let mut cell_types = BTreeMap::<String, String>::new();
    let annotation_path = out.join("cluster_ai_annotations.tsv");
    if annotation_path.is_file() {
        for line in fs::read_to_string(&annotation_path)?.lines().skip(1) {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() >= 2 {
                cell_types.insert(fields[0].to_owned(), fields[1].to_owned());
            }
        }
    }
    let mut by_module: BTreeMap<usize, Vec<&str>> = BTreeMap::new();
    for (gene, module) in modules {
        by_module.entry(*module).or_default().push(gene);
    }
    let mut tsv = BufWriter::new(File::create(
        out.join(format!("module_trajectories{suffix}.tsv")),
    )?);
    // One row per module and category: both plotted groups and the comparison
    // are kept together, avoiding duplicated test results in a long-form table.
    writeln!(tsv, "module\tcluster\tgenes\tcontrol_label\tcase_label\tcontrol_mean_log1p\tcontrol_sd\tcontrol_se\tcontrol_ci95_low\tcontrol_ci95_high\tcontrol_count\tcase_mean_log1p\tcase_sd\tcase_se\tcase_ci95_low\tcase_ci95_high\tcase_count\tactual_control_count\tactual_case_count\tdelta_case_minus_control\twelch_t\twelch_df\tp_raw\tp_bonferroni\tbonferroni_factor")?;
    let quantified = crate::module_score::score(data, modules)?;
    // One correction family covers every module x phase in this output run.
    // Include untestable combinations in the factor (conservative and stable).
    let n_modules = by_module.keys().filter(|&&id| id != 0).count();
    let bonferroni_factor = n_modules * ordered.len();
    for (module, genes) in by_module {
        if module == 0 {
            continue;
        }
        let module_index = quantified.ids.iter().position(|&id| id == module).expect("module id");
        let mut values = vec![vec![Vec::<f64>::new(); 2]; ordered.len()];
        for (i, cluster) in ordered.iter().enumerate() {
            let (cases, controls) = &groups[cluster];
            for (j, cells) in [controls, cases].iter().enumerate() {
                for &cell in *cells { values[i][j].push(quantified.module_means[cell][module_index]); }
            }
        }
        let mut means = vec![[f64::NAN; 2]; ordered.len()];
        let mut lows = means.clone();
        let mut highs = means.clone();
        let tests: Vec<_> = values.iter().map(|v| welch_test(&v[0], &v[1])).collect();
        for (i, cluster) in ordered.iter().enumerate() {
            let mut stats = Vec::with_capacity(2);
            for j in 0..2 {
                let v = &values[i][j];
                let (mean, sd, se, low, high) = cell_stats(v);
                means[i][j] = mean;
                lows[i][j] = low;
                highs[i][j] = high;
                stats.push([fmt(mean), fmt(sd), fmt(se), fmt(low), fmt(high), v.len().to_string()]);
            }
            let test = tests[i];
            let adjusted = test.map(|(_, _, _, p)| (p * bonferroni_factor as f64).min(1.0));
            let (actual_cases, actual_controls) = actual_groups.get(cluster)
                .map(|(cases, controls)| (cases.len(), controls.len()))
                .unwrap_or((0, 0));
            let mut fields = vec![module.to_string(), cluster.clone(), genes.len().to_string(), control_label.to_owned(), case_label.to_owned()];
            fields.extend(stats[0].iter().cloned());
            fields.extend(stats[1].iter().cloned());
            fields.extend([actual_controls.to_string(), actual_cases.to_string()]);
            for value in [test.map(|t| t.0), test.map(|t| t.1), test.map(|t| t.2), test.map(|t| t.3), adjusted] {
                fields.push(value.map_or("NA".to_owned(), fmt));
            }
            fields.push(bonferroni_factor.to_string());
            writeln!(tsv, "{}", fields.join("\t"))?;
        }
        let width = (ordered.len().max(2) as f64 * 115.0 + 190.0).max(660.0);
        // Long population names need vertical space below the plot. Rotate
        // only when necessary so short Leiden/Louvain IDs stay horizontal.
        let angle = label_angle.unwrap_or_else(|| {
            if ordered.iter().any(|name| name.chars().count() > 14) {
                55.0
            } else {
                0.0
            }
        });
        anyhow::ensure!(
            angle.is_finite() && (0.0..=85.0).contains(&angle),
            "--label-angle must be in 0..=85"
        );
        anyhow::ensure!(
            matches!(count_label_layout, "horizontal" | "rotated"),
            "invalid --count-label-layout"
        );
        let rotate_labels = angle > 0.0;
        let longest = ordered.iter().map(|s| s.chars().count()).max().unwrap_or(0) as f64;
        let label_space = if rotate_labels {
            (longest * 7.0 * angle.to_radians().sin() + 35.0).ceil()
        } else {
            55.0
        };
        let height = (315.0 + label_space + 150.0).max(507.0);
        let (left, right, top, bottom) = (85.0, width - 40.0, 75.0, 315.0);
        let x = |i: usize| left + (right - left) * (i as f64 + 0.5) / ordered.len() as f64;
        let ymax = highs.iter().chain(means.iter()).flat_map(|v| v.iter().copied()).filter(|v| v.is_finite()).fold(0.0, f64::max) * 1.05 + 0.01;
        let ymin = 0.0;
        let y = |v: f64| bottom - (v - ymin) / (ymax - ymin).max(1e-10) * (bottom - top);
        let mut svg = BufWriter::new(File::create(
            out.join(format!("module_{module}_trajectory{suffix}.svg")),
        )?);
        writeln!(
            svg,
            "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
        )?;
        writeln!(
            svg,
            "<text x=\"{left}\" y=\"30\" font-family=\"sans-serif\" font-size=\"19\" font-weight=\"bold\">Module {module} · {} genes ({} informative)</text>",
            genes.len(),
            genes.iter().filter(|g| data.feature_index(g).is_some()).count()
        )?;
        writeln!(
            svg,
            "<text x=\"{left}\" y=\"51\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#555\">Mean original log1p expression of module genes (including zeros)</text>"
        )?;
        let tick_step = (ymax - ymin) / 4.0;
        let precision = (0..=8).find(|&digits| tick_step * 10f64.powi(digits) >= 2.0).unwrap_or(8) as usize;
        for tick in 0..5 {
            let val = ymin + (ymax - ymin) * tick as f64 / 4.0;
            writeln!(
                svg,
                "<line x1=\"{left}\" y1=\"{}\" x2=\"{right}\" y2=\"{}\" stroke=\"#eee\"/><text x=\"{}\" y=\"{}\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\">{val:.precision$}</text>",
                y(val),
                y(val),
                left - 9.0,
                y(val) + 4.0
            )?;
        }
        // Actual within-cluster counts, not the possibly substituted WT
        // reference used for computing the expression trajectories.
        for (i, cluster) in ordered.iter().enumerate() {
            let (case_n, wt_n) = actual_groups
                .get(cluster)
                .map(|(cases, controls)| (cases.len(), controls.len()))
                .unwrap_or((0, 0));
            let source = reference_clusters.get(cluster);
            let substituted = source.is_some_and(|reference| reference != cluster);
            let marker = if substituted { "*" } else { "" };
            // Each count is anchored to its own cluster label, never to the SVG bottom.
            let label_y = bottom + 63.0;
            let count_y = bottom + 18.0;
            let p_y = bottom + 39.0;
            if count_label_layout == "rotated" && rotate_labels {
                writeln!(
                    svg,
                    "<g transform=\"translate({:.1},{:.1}) rotate(-{angle:.1})\" text-anchor=\"end\" font-family=\"sans-serif\">",
                    x(i),
                    label_y
                )?;
                writeln!(
                    svg,
                    "<text x=\"0\" y=\"-45\" font-size=\"10\"><tspan fill=\"#2563eb\">{wt_n}</tspan><tspan fill=\"#555\"> / </tspan><tspan fill=\"#dc2626\">{case_n}</tspan></text>"
                )?;
                writeln!(
                    svg,
                    "<text x=\"0\" y=\"0\" font-size=\"11\">{}{}</text></g>",
                    escape(cluster),
                    marker
                )?;
            } else {
                writeln!(
                    svg,
                    "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"10\"><tspan fill=\"#2563eb\">{wt_n}</tspan><tspan fill=\"#555\"> / </tspan><tspan fill=\"#dc2626\">{case_n}</tspan></text>",
                    x(i),
                    count_y
                )?;
                if rotate_labels {
                    writeln!(
                        svg,
                        "<text transform=\"translate({:.1},{:.1}) rotate(-{angle:.1})\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\">{}{}</text>",
                        x(i),
                        label_y,
                        escape(cluster),
                        marker
                    )?;
                } else {
                    writeln!(
                        svg,
                        "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\">{}{}</text>",
                        x(i),
                        label_y,
                        escape(cluster),
                        marker
                    )?;
                }
            }
            let adjusted = tests[i].map(|(_, _, _, p)| (p * bonferroni_factor as f64).min(1.0));
            writeln!(svg,
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#333\">{}</text>",
                x(i), p_y, p_label(adjusted))?;
            if let Some(cell_type) = cell_types.get(cluster) {
                // SVG title retains the full label even when the visible label
                // needs truncation to avoid colliding with adjacent clusters.
                let short: String = cell_type.chars().take(19).collect();
                let short = if cell_type.chars().count() > 19 {
                    format!("{short}…")
                } else {
                    short
                };
                writeln!(
                    svg,
                    "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"9\" fill=\"#444\"><title>{}</title>{}</text>",
                    x(i),
                    if rotate_labels {
                        bottom + label_space + 64.0
                    } else {
                        bottom + 82.0
                    },
                    escape(cell_type),
                    escape(&short)
                )?;
            }
        }
        for (j, color, label) in [(0, "#2563eb", control_label), (1, "#dc2626", case_label)] {
            // Draw contiguous observed segments only; missing populations remain gaps.
            let mut segment: Vec<(f64, f64)> = Vec::new();
            for (i, row) in means.iter().enumerate() {
                if row[j].is_finite() {
                    segment.push((x(i), y(row[j])));
                    if lows[i][j].is_finite() && highs[i][j].is_finite() {
                        let xx = x(i);
                        let yl = y(lows[i][j].max(0.0));
                        let yh = y(highs[i][j]);
                        writeln!(svg, "<path d=\"M{xx:.1},{yl:.1} V{yh:.1} M{:.1},{yl:.1} H{:.1} M{:.1},{yh:.1} H{:.1}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"1.5\"/>", xx-5.0, xx+5.0, xx-5.0, xx+5.0)?;
                    }
                    writeln!(svg, "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"4\" fill=\"{color}\"/>", x(i), y(row[j]))?;
                } else {
                    if segment.len() > 1 {
                        writeln!(svg, "<polyline points=\"{}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2.5\"/>", polyline(&segment))?;
                    }
                    segment.clear();
                }
            }
            if segment.len() > 1 {
                writeln!(svg, "<polyline points=\"{}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2.5\"/>", polyline(&segment))?;
            }
            let lx = left + j as f64 * 190.0;
            writeln!(svg, "<line x1=\"{lx}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{color}\" stroke-width=\"3\"/><text x=\"{}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"13\">{}</text>", height - 33.0, lx + 26.0, height - 33.0, lx + 34.0, height - 29.0, escape(label))?;
        }
        writeln!(
            svg,
            "<text x=\"{left}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#666\">Below clusters: observed <tspan fill=\"#2563eb\">{}</tspan> / <tspan fill=\"#dc2626\">{}</tspan> cell counts. Bonferroni correction: {bonferroni_factor} tests ({n_modules} modules × {} phases).</text>",
            height - 58.0,
            escape(control_label),
            escape(case_label),
            ordered.len()
        )?;
        writeln!(
            svg,
            "<text x=\"{left}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#666\">Clusters are categorical; connecting lines do not imply developmental order.</text></svg>",
            height - 8.0
        )?;
    }
    tsv.flush()?;
    Ok(())
}
