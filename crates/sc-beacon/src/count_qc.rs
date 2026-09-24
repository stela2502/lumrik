use std::path::Path;

use anyhow::{Context, Result, bail};
use plotters::prelude::*;

use crate::KneeCountFit;

/// Write QC plots for knee-based cell calling from the complete positive barcode distribution.
pub fn write_count_qc<P: AsRef<Path>>(
    counts: &[u32],
    fit: &KneeCountFit,
    out_dir: P,
) -> Result<()> {
    if counts.is_empty() {
        bail!("cannot plot cell-count QC for an empty count distribution");
    }
    let out_dir = out_dir.as_ref();
    std::fs::create_dir_all(out_dir).with_context(|| format!("creating {}", out_dir.display()))?;
    write_log_umi_distribution(counts, fit, out_dir.join("cell_umi_distribution.svg"))?;
    write_rank_plot(counts, fit, out_dir.join("cell_umi_rank.svg"))?;
    Ok(())
}

fn write_log_umi_distribution<P: AsRef<Path>>(
    counts: &[u32],
    fit: &KneeCountFit,
    path: P,
) -> Result<()> {
    const WIDTH: u32 = 1200;
    const HEIGHT: u32 = 720;
    const BINS: usize = 100;

    let logs: Vec<f64> = counts
        .iter()
        .copied()
        .filter(|&n| n > 0)
        .map(|n| (n as f64).log10())
        .collect();
    let max_x = logs.iter().copied().fold(0.0_f64, f64::max).max(1.0);
    let bin_width = max_x / BINS as f64;
    let mut bins = vec![0usize; BINS];
    for &x in &logs {
        let idx = ((x / bin_width).floor() as usize).min(BINS - 1);
        bins[idx] += 1;
    }
    let max_y = bins.iter().copied().max().unwrap_or(1) as f64;
    let cutoff_x = (fit.umi_cutoff as f64).log10();

    let root = SVGBackend::new(path.as_ref(), (WIDTH, HEIGHT)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Unfiltered exonic UMI distribution", ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(55)
        .y_label_area_size(75)
        .build_cartesian_2d(0.0..max_x, 0.0..(max_y * 1.08))?;
    chart
        .configure_mesh()
        .x_desc("log10(exonic UMIs per barcode)")
        .y_desc("Barcodes")
        .draw()?;
    chart.draw_series(bins.iter().enumerate().map(|(i, &n)| {
        let x0 = i as f64 * bin_width;
        Rectangle::new(
            [(x0, 0.0), (x0 + bin_width, n as f64)],
            BLUE.mix(0.35).filled(),
        )
    }))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(cutoff_x, 0.0), (cutoff_x, max_y * 1.05)],
        RED.stroke_width(2),
    )))?;
    chart.draw_series(std::iter::once(Text::new(
        format!("knee cutoff: {} UMIs", fit.umi_cutoff),
        (cutoff_x, max_y * 1.02),
        ("sans-serif", 18).into_font(),
    )))?;
    root.present()
        .with_context(|| format!("writing {}", path.as_ref().display()))?;
    Ok(())
}

fn write_rank_plot<P: AsRef<Path>>(counts: &[u32], fit: &KneeCountFit, path: P) -> Result<()> {
    const WIDTH: u32 = 1200;
    const HEIGHT: u32 = 720;
    let mut sorted: Vec<u32> = counts.iter().copied().filter(|&n| n > 0).collect();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let max_rank = (sorted.len() as f64).log10().max(1.0);
    let max_count = (*sorted.first().unwrap_or(&1) as f64).log10().max(1.0);
    let knee_x = (fit.knee_rank as f64).log10();
    let knee_y = (fit.umi_cutoff as f64).log10();

    let root = SVGBackend::new(path.as_ref(), (WIDTH, HEIGHT)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Unfiltered barcode rank", ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(55)
        .y_label_area_size(75)
        .build_cartesian_2d(0.0..max_rank, 0.0..(max_count * 1.05))?;
    chart
        .configure_mesh()
        .x_desc("log10(barcode rank)")
        .y_desc("log10(exonic UMIs)")
        .draw()?;
    chart.draw_series(LineSeries::new(
        sorted
            .iter()
            .enumerate()
            .map(|(i, &n)| (((i + 1) as f64).log10(), (n as f64).log10())),
        &BLUE,
    ))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(knee_x, 0.0), (knee_x, max_count * 1.03)],
        RED.stroke_width(2),
    )))?;
    chart.draw_series(std::iter::once(PathElement::new(
        vec![(0.0, knee_y), (max_rank, knee_y)],
        RED.mix(0.55).stroke_width(1),
    )))?;
    chart.draw_series(std::iter::once(Circle::new(
        (knee_x, knee_y),
        5,
        RED.filled(),
    )))?;
    chart.draw_series(std::iter::once(Text::new(
        format!("knee: rank {}, {} UMIs", fit.knee_rank, fit.umi_cutoff),
        (knee_x, (knee_y + 0.08).min(max_count)),
        ("sans-serif", 18).into_font(),
    )))?;
    root.present()
        .with_context(|| format!("writing {}", path.as_ref().display()))?;
    Ok(())
}
