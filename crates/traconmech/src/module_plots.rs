//! Self-contained SVG summaries of gene-normalized module expression.
//! Clusters are categorical unless the user supplied an explicit biological order.
use anyhow::Result;
use sc_analysis::SingleCellData;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn median(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return 0.0; }
    let pos = p * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
        .replace('"', "&quot;").replace('\'', "&apos;")
}
fn polyline(points: &[(f64, f64)]) -> String {
    points.iter().map(|(x,y)| format!("{x:.1},{y:.1}"))
        .collect::<Vec<_>>().join(" ")
}

/// Each gene is log1p-transformed and z-scored *jointly* across both conditions
/// and all displayed clusters. Thus WT/EXP separation is not normalized away.
/// Bands show the interquartile range across genes, not uncertainty in cells.
pub fn write_module_plots(
    out: &Path, data: &SingleCellData,
    modules: &BTreeMap<String, usize>,
    ordered: &[String],
    groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    actual_groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    reference_clusters: &BTreeMap<String, String>,
    case_label: &str, control_label: &str,
) -> Result<()> {
    write_module_plots_named(out, data, modules, ordered, groups, actual_groups,
        reference_clusters, case_label, control_label, "")
}

pub fn write_module_plots_named(
    out: &Path, data: &SingleCellData,
    modules: &BTreeMap<String, usize>, ordered: &[String],
    groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    actual_groups: &BTreeMap<String, (Vec<usize>, Vec<usize>)>,
    reference_clusters: &BTreeMap<String, String>,
    case_label: &str, control_label: &str, suffix: &str,
) -> Result<()> {
    let mut by_module: BTreeMap<usize, Vec<&str>> = BTreeMap::new();
    for (gene, module) in modules { by_module.entry(*module).or_default().push(gene); }
    let mut tsv = BufWriter::new(File::create(out.join(format!("module_trajectories{suffix}.tsv")))?);
    writeln!(tsv, "module\tcluster\tcondition\tgenes\tmedian_z\tq25_z\tq75_z\tcase_n\tcontrol_n")?;
    for (module, genes) in by_module {
        if module == 0 { continue; }
        let mut values = vec![vec![Vec::<f64>::new(); 2]; ordered.len()];
        for gene in &genes {
            let Some(idx) = data.feature_index(gene) else { continue; };
            let Some(row) = data.matrix.outer_view(idx) else { continue; };
            let mut v = Vec::with_capacity(ordered.len() * 2);
            for cluster in ordered {
                let (cases, controls) = &groups[cluster];
                for cells in [controls, cases] {
                    let mean = if cells.is_empty() { 0.0 } else {
                        cells.iter().map(|&cell| row.get(cell).copied().unwrap_or(0.0) as f64)
                            .sum::<f64>() / cells.len() as f64
                    };
                    v.push(mean.ln_1p());
                }
            }
            let avg = v.iter().sum::<f64>() / v.len() as f64;
            let sd = (v.iter().map(|x| (x-avg).powi(2)).sum::<f64>() / v.len() as f64).sqrt();
            // A constant gene cannot contribute a trajectory; retain it in
            // module membership, but omit it from normalized visual summaries.
            if sd < 1e-12 { continue; }
            for (i, x) in v.into_iter().enumerate() {
                values[i / 2][i % 2].push((x - avg) / sd);
            }
        }
        let mut bands = vec![vec![(0.0, 0.0, 0.0); 2]; ordered.len()];
        for (i, cluster) in ordered.iter().enumerate() {
            let (cases, controls) = &groups[cluster];
            for j in 0..2 {
                let v = &mut values[i][j];
                v.sort_by(f64::total_cmp);
                let (q25, med, q75) = (median(v, 0.25), median(v, 0.5), median(v, 0.75));
                bands[i][j] = (q25, med, q75);
                writeln!(tsv, "{module}\t{cluster}\t{}\t{}\t{med}\t{q25}\t{q75}\t{}\t{}",
                    if j == 0 { control_label } else { case_label }, v.len(), cases.len(), controls.len())?;
            }
        }
        let width = (ordered.len().max(2) as f64 * 115.0 + 190.0).max(660.0);
        let height = 432.0;
        let (left, right, top, bottom) = (85.0, width - 40.0, 75.0, 315.0);
        let x = |i: usize| left + (right-left) * (i as f64 + 0.5) / ordered.len() as f64;
        let ymin = bands.iter().flat_map(|a| a.iter().map(|v| v.0)).fold(0.0, f64::min) - 0.15;
        let ymax = bands.iter().flat_map(|a| a.iter().map(|v| v.2)).fold(0.0, f64::max) + 0.15;
        let y = |v: f64| bottom - (v-ymin)/(ymax-ymin).max(1e-10)*(bottom-top);
        let mut svg = BufWriter::new(File::create(out.join(format!("module_{module}_trajectory{suffix}.svg")))?);
        writeln!(svg, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/>")?;
        writeln!(svg, "<text x=\"{left}\" y=\"30\" font-family=\"sans-serif\" font-size=\"19\" font-weight=\"bold\">Module {module} · {} genes ({} informative)</text>", genes.len(), values[0][0].len())?;
        writeln!(svg, "<text x=\"{left}\" y=\"51\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#555\">Gene-wise joint log1p z-score · median and interquartile range</text>")?;
        for tick in 0..5 {
            let val = ymin + (ymax-ymin) * tick as f64 / 4.0;
            writeln!(svg, "<line x1=\"{left}\" y1=\"{}\" x2=\"{right}\" y2=\"{}\" stroke=\"#eee\"/><text x=\"{}\" y=\"{}\" text-anchor=\"end\" font-family=\"sans-serif\" font-size=\"11\">{val:.1}</text>",y(val),y(val),left-9.0,y(val)+4.0)?;
        }
        // Actual within-cluster counts, not the possibly substituted WT
        // reference used for computing the expression trajectories.
        for (i, cluster) in ordered.iter().enumerate() {
            let (case_n, wt_n) = actual_groups.get(cluster)
                .map(|(cases, controls)| (cases.len(), controls.len()))
                .unwrap_or((0, 0));
            let source = reference_clusters.get(cluster);
            let substituted = source.is_some_and(|reference| reference != cluster);
            let marker = if substituted { "*" } else { "" };
            writeln!(svg, "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\">{}{}</text>",
                x(i), bottom+20.0, escape(cluster), marker)?;
            // Separate colored counts remain legible even for 4-digit values.
            writeln!(svg, "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"10\"><tspan fill=\"#2563eb\">{wt_n}</tspan><tspan fill=\"#555\"> / </tspan><tspan fill=\"#dc2626\">{case_n}</tspan></text>",
                x(i), bottom+35.0)?;
        }
        for (j, color, label) in [(0, "#2563eb", control_label), (1, "#dc2626", case_label)] {
            let upper: Vec<_> = bands.iter().enumerate().map(|(i,v)| (x(i),y(v[j].2))).collect();
            let lower: Vec<_> = bands.iter().enumerate().rev().map(|(i,v)| (x(i),y(v[j].0))).collect();
            let mut polygon = upper;
            polygon.extend(lower);
            writeln!(svg, "<polygon points=\"{}\" fill=\"{color}\" opacity=\"0.13\"/>",polyline(&polygon))?;
            let line: Vec<_> = bands.iter().enumerate().map(|(i,v)| (x(i),y(v[j].1))).collect();
            writeln!(svg, "<polyline points=\"{}\" fill=\"none\" stroke=\"{color}\" stroke-width=\"2.5\"/>",polyline(&line))?;
            for (px,py) in line { writeln!(svg, "<circle cx=\"{px:.1}\" cy=\"{py:.1}\" r=\"3\" fill=\"{color}\"/>")?; }
            let lx = left + j as f64 * 190.0;
            writeln!(svg, "<line x1=\"{lx}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{color}\" stroke-width=\"3\"/><text x=\"{}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"13\">{}</text>",height-33.0,lx+26.0,height-33.0,lx+34.0,height-29.0,escape(label))?;
        }
        writeln!(svg, "<text x=\"{left}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#666\">Below clusters: actual <tspan fill=\"#2563eb\">{}</tspan> / <tspan fill=\"#dc2626\">{}</tspan> cell counts. * WT reference substituted (see wt_reference_clusters.tsv).</text>", height-58.0, escape(control_label), escape(case_label))?;
        writeln!(svg, "<text x=\"{left}\" y=\"{}\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#666\">Clusters are categorical; connecting lines do not imply developmental order.</text></svg>",height-8.0)?;
    }
    tsv.flush()?;
    Ok(())
}
