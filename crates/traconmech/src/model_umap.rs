//! Frozen-module UMAPs: Seurat-style background-corrected module scores.
use anyhow::{Context, Result, ensure};
use hdf5::File as H5File;
use plotters::prelude::*;
use sc_analysis::SingleCellData;
use std::{collections::{BTreeMap, HashMap, HashSet}, fs::File, io::{BufRead, BufReader}, path::Path};

const PLOT_W: u32 = 1000;
const PLOT_H: u32 = 900;
const LEGEND_W: u32 = 300;
const ROWS_PER_LEGEND: usize = 28;

// Matplotlib viridis_r sampled at 16 evenly spaced stops (high scores = purple).
// Linear RGB interpolation keeps this dependency-free and deterministic.
fn viridis_r(t: f64) -> RGBColor {
    const STOPS: [(u8, u8, u8); 16] = [
        (253, 231, 37), (221, 227, 24), (181, 222, 43), (144, 215, 67),
        (110, 206, 88), (81, 196, 104), (57, 185, 119), (39, 173, 129),
        (31, 158, 137), (33, 142, 141), (40, 125, 142), (49, 104, 142),
        (59, 82, 139), (68, 59, 131), (72, 35, 116), (68, 1, 84),
    ];
    let pos = t.clamp(0.0, 1.0) * (STOPS.len() - 1) as f64;
    let i = (pos.floor() as usize).min(STOPS.len() - 2);
    let f = pos - i as f64;
    let a = STOPS[i];
    let b = STOPS[i + 1];
    let mix = |x: u8, y: u8| (x as f64 + f * (y as f64 - x as f64)).round() as u8;
    RGBColor(mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

fn safe_name(name: &str) -> String {
    name.chars().map(|c| if c.is_ascii_alphanumeric() || c == '-' || c == '_' { c } else { '_' }).collect()
}

fn coords(path: &Path, n: usize) -> Result<Vec<(f64, f64)>> {
    let file = H5File::open(path).with_context(|| format!("opening H5AD {}", path.display()))?;
    let ds = file.dataset("obsm/X_umap").context("H5AD missing /obsm/X_umap")?;
    let shape = ds.shape();
    ensure!(shape.len() == 2 && shape[0] == n && shape[1] >= 2, "UMAP dimensions {:?} do not match {} cells", shape, n);
    let raw: Vec<f64> = if ds.dtype()?.is::<f32>() {
        ds.read_raw::<f32>()?.into_iter().map(f64::from).collect()
    } else {
        ds.read_raw::<f64>()?
    };
    let xy: Vec<_> = raw.chunks_exact(shape[1]).map(|r| (r[0], r[1])).collect();
    ensure!(xy.iter().all(|(x,y)| x.is_finite() && y.is_finite()), "non-finite UMAP coordinates");
    Ok(xy)
}

fn bounds(xy: &[(f64, f64)]) -> ((f64,f64),(f64,f64)) {
    let (mut x0,mut x1,mut y0,mut y1) = (f64::INFINITY,f64::NEG_INFINITY,f64::INFINITY,f64::NEG_INFINITY);
    for &(x,y) in xy { x0=x0.min(x); x1=x1.max(x); y0=y0.min(y); y1=y1.max(y); }
    let dx = (x1-x0).max(1e-6)*0.05;
    let dy = (y1-y0).max(1e-6)*0.05;
    ((x0-dx,x1+dx),(y0-dy,y1+dy))
}

fn draw_base<'a, 'b>(
    area: &'a DrawingArea<SVGBackend<'b>, plotters::coord::Shift>,
    title: &str, limits: ((f64,f64),(f64,f64)),
) -> Result<ChartContext<'a, SVGBackend<'b>, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>> {
    let ((x0,x1),(y0,y1)) = limits;
    let mut chart = ChartBuilder::on(area).caption(title, ("sans-serif", 27))
        .margin(28).x_label_area_size(35).y_label_area_size(40)
        .build_cartesian_2d(x0..x1,y0..y1)?;
    chart.configure_mesh().x_desc("UMAP 1").y_desc("UMAP 2").disable_mesh().draw()?;
    Ok(chart)
}

pub fn write(data: &SingleCellData, h5ad: &Path, out: &Path, grouping: &str) -> Result<()> {
    let xy = coords(h5ad, data.cells.len())?;
    let labels = data.annotations.column(grouping)
        .with_context(|| format!("missing UMAP grouping annotation {grouping}"))?;
    ensure!(labels.len() == xy.len(), "UMAP grouping length mismatch");
    let limits = bounds(&xy);
    let folder = out.join("umap");
    std::fs::create_dir_all(&folder)?;
    let mut counts = BTreeMap::<String, usize>::new();
    for label in labels { *counts.entry(label.clone()).or_default() += 1; }
    let names: Vec<String> = counts.keys().cloned().collect();
    let n_legend = names.len().div_ceil(ROWS_PER_LEGEND).max(1) as u32;
    let label_file = folder.join(format!("{}.svg", safe_name(grouping)));
    {
        let root = SVGBackend::new(&label_file, (PLOT_W + n_legend*LEGEND_W,PLOT_H)).into_drawing_area();
        root.fill(&WHITE)?;
        let (plot, legend) = root.split_horizontally(PLOT_W);
        let mut chart = draw_base(&plot, &format!("UMAP — {grouping}"), limits)?;
        let lookup: HashMap<&str, usize> = names.iter().enumerate().map(|(i,s)| (s.as_str(),i)).collect();
        chart.draw_series(xy.iter().zip(labels).map(|(&(x,y),label)| {
            let color = Palette99::pick(lookup[label.as_str()]);
            Circle::new((x,y),2,color.filled())
        }))?;
        for (i,name) in names.iter().enumerate() {
            let col = (i / ROWS_PER_LEGEND) as i32;
            let row = (i % ROWS_PER_LEGEND) as i32;
            let x = 20 + col*LEGEND_W as i32;
            let y = 35 + row*30;
            legend.draw(&Circle::new((x,y),6,Palette99::pick(i).filled()))?;
            legend.draw(&Text::new(format!("{} ({})",name,counts[name]),(x+15,y+5),("sans-serif",14).into_font()))?;
        }
        root.present()?;
    }
    let scores_file = out.join("gene_model_scores.tsv");
    let mut lines = BufReader::new(File::open(&scores_file).with_context(|| format!("opening {}",scores_file.display()))?).lines();
    let header = lines.next().context("empty model score file")??;
    let columns: Vec<_> = header.split('\t').collect();
    ensure!(columns.len() >= 4 && columns[0] == "cell" && columns[columns.len()-2..] == ["best_module","condition"], "unexpected gene_model_scores.tsv header");
    let model_names = &columns[1..columns.len()-2];
    let mut seen = HashSet::new();
    let mut scores: Vec<Vec<Option<f64>>> = vec![vec![None; xy.len()];model_names.len()];
    let cell_lookup: HashMap<&str,usize> = data.cells.iter().enumerate().map(|(i,s)| (s.as_str(),i)).collect();
    let mut matched = 0;
    for line in lines {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        ensure!(fields.len() == columns.len(), "invalid score row column count");
        let Some(&cell) = cell_lookup.get(fields[0]) else { continue; };
        ensure!(seen.insert(cell), "duplicate model score for cell {}",fields[0]);
        matched += 1;
        for (j,value) in fields[1..fields.len()-2].iter().enumerate() {
            if *value != "NA" && !value.is_empty() {
                let score: f64 = value.parse().with_context(|| format!("invalid score {value}"))?;
                ensure!(score.is_finite(), "invalid module score {score}");
                scores[j][cell] = Some(score);
            }
        }
    }
    ensure!(matched == xy.len(), "only {matched}/{} model scores match H5AD cell IDs",xy.len());
    for (j,name) in model_names.iter().enumerate() {
        let path = folder.join(format!("{}.svg",safe_name(name)));
        ensure!(path != label_file, "model name conflicts with grouping plot filename");
        let root = SVGBackend::new(&path,(PLOT_W+210,PLOT_H)).into_drawing_area();
        root.fill(&WHITE)?;
        let (plot,legend) = root.split_horizontally(PLOT_W);
        let mut chart = draw_base(&plot,&format!("UMAP — {name} module score"),limits)?;
        // Missing scores appear as light grey points; draw them first.
        chart.draw_series(xy.iter().zip(&scores[j]).filter(|(_,s)| s.is_none())
            .map(|(&(x,y),_)| Circle::new((x,y),2,RGBColor(205,205,205).filled())))?;
        let observed: Vec<f64> = scores[j].iter().flatten().copied().collect();
        // Display-only clipping: keep original scores unchanged in the TSV.
        // Use linear interpolation between the 2nd and 98th percentiles.
        let mut sorted = observed;
        sorted.sort_by(|a, b| a.total_cmp(b));
        let percentile = |p: f64| -> f64 {
            if sorted.is_empty() { return 0.0; }
            let pos = p * (sorted.len() - 1) as f64;
            let lo = pos.floor() as usize;
            let hi = pos.ceil() as usize;
            sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
        };
        let min = percentile(0.02);
        let max = percentile(0.98);
        let range = (max - min).max(f64::EPSILON);
        chart.draw_series(xy.iter().zip(&scores[j]).filter_map(|(&(x,y),s)| s.map(|v| (x,y,v)))
            .map(|(x,y,v)| Circle::new((x,y),2,viridis_r((v-min)/range).filled())))?;
        let (x0,x1,y0,y1)=(45,85,135,730);
        for i in 0..120 {
            let t = i as f64/119.0;
            let ya = y1 - ((i+1)*(y1-y0)/120);
            let yb = y1 - (i*(y1-y0)/120);
            legend.draw(&Rectangle::new([(x0,ya),(x1,yb)],viridis_r(t).filled()))?;
        }
        for (v,y) in [(max,y0),((min+max)/2.0,(y0+y1)/2),(min,y1)] {
            legend.draw(&Text::new(format!("{v:.2}"),(x1+12,y+5),("sans-serif",17).into_font()))?;
        }
        legend.draw(&Text::new("Module score (p2–p98)",(15,95),("sans-serif",18).into_font()))?;
        let missing = scores[j].iter().filter(|v| v.is_none()).count();
        if missing > 0 {
            legend.draw(&Circle::new((52,785),6,RGBColor(205,205,205).filled()))?;
            legend.draw(&Text::new(format!("NA ({missing})"),(67,790),("sans-serif",15).into_font()))?;
        }
        root.present()?;
    }
    println!("Wrote {} model UMAPs and {} grouping plot to {}",model_names.len(),grouping,folder.display());
    Ok(())
}
