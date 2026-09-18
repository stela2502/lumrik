use crate::CloneData;
use ndarray::ArrayView1;
#[allow(dead_code, unused)] // creates a warning otherwise
use plotters::prelude::*;
use std::collections::VecDeque;
use std::path::Path;

use crate::PcaModel;

pub struct MstTree {
    pub edges: Vec<(usize, usize, f32)>,
}


struct RootedLayout {
    parent: Vec<usize>,
    depth: Vec<usize>,
    path_distance: Vec<f32>,
    order: Vec<usize>,
    children: Vec<Vec<usize>>,
}

fn categorical_palette() -> [RGBColor; 6] {
    // Keep the actual colours brutally distinct.  Extra categorical capacity
    // comes from rotating the split axis, not from adding hard-to-distinguish
    // shades.
    [
        RGBColor(0, 0, 0),       // black
        RGBColor(255, 255, 255), // white
        RGBColor(230, 159, 0),   // orange
        RGBColor(0, 114, 178),   // blue
        RGBColor(255, 0, 0),     // bright red
        RGBColor(190, 120, 255), // light purple
    ]
}

fn categorical_code(i: usize, palette: &[RGBColor]) -> (RGBColor, RGBColor, usize) {
    // Use ordered colour pairs and four visibly different split axes:
    // |, /, -, \.  Swapping the two colours is itself a visible code, so
    // six strong colours provide 6 * 5 * 4 = 120 categorical symbols.
    let mut pairs = Vec::new();
    for a in 0..palette.len() {
        for b in 0..palette.len() {
            if a != b {
                pairs.push((palette[a], palette[b]));
            }
        }
    }
    let pair = pairs[(i / 4) % pairs.len()];
    (pair.0, pair.1, i % 4)
}

fn split_tilt(axis: usize) -> f64 {
    // Angle of the diameter separating the two coloured halves.
    match axis % 4 {
        0 => 90.0f64.to_radians(), // |
        1 => 45.0f64.to_radians(), // /
        2 => 0.0f64,               // -
        _ => -45.0f64.to_radians(),// \
    }
}

fn abundance_radius(n: usize, max_n: usize) -> i32 {
    (5.0 + 10.0 * ((n.max(1) as f64 / max_n.max(1) as f64).sqrt())).round() as i32
}

fn blue_yellow_red(value: f32, min_value: f32, max_value: f32) -> RGBColor {
    let t = if max_value > min_value {
        ((value - min_value) / (max_value - min_value)).clamp(0.0, 1.0)
    } else {
        0.5
    };
    let lerp = |a: u8, b: u8, u: f32| -> u8 {
        (a as f32 + (b as f32 - a as f32) * u).round() as u8
    };
    if t <= 0.5 {
        let u = t * 2.0;
        RGBColor(lerp(0, 255, u), lerp(114, 220, u), lerp(178, 0, u))
    } else {
        let u = (t - 0.5) * 2.0;
        RGBColor(255, lerp(220, 0, u), 0)
    }
}

fn rgb_hex(c: RGBColor) -> String {
    format!("#{:02X}{:02X}{:02X}", c.0, c.1, c.2)
}

/// Return the exact categorical colour codes used by the rooted SVG renderer.
/// Split nodes are represented as `#RRGGBB/#RRGGBB`; the split orientation is
/// deliberately not folded into the colour string.
pub fn rooted_categorical_hex(
    categories: &[String],
    category_order: Option<&[String]>,
    root_node: usize,
) -> Vec<String> {
    use std::collections::{BTreeMap, BTreeSet};
    let cats: BTreeSet<String> = categories.iter().enumerate()
        .filter(|(i, _)| *i != root_node)
        .map(|(_, c)| c.clone())
        .collect();
    let mut ordered = Vec::with_capacity(cats.len());
    if let Some(preferred) = category_order {
        for c in preferred {
            if cats.contains(c) && !ordered.contains(c) { ordered.push(c.clone()); }
        }
    }
    for c in &cats {
        if !ordered.contains(c) { ordered.push(c.clone()); }
    }
    let palette = categorical_palette();
    let cmap: BTreeMap<String, (RGBColor, RGBColor, usize)> = ordered.into_iter().enumerate()
        .map(|(i, c)| {
            if c == "unpaired" || c == "HC NAIVE" {
                (c, (RGBColor(170,170,170), RGBColor(170,170,170), 0))
            } else {
                (c, categorical_code(i, &palette))
            }
        }).collect();
    categories.iter().enumerate().map(|(i, c)| {
        if i == root_node { return "#AAAAAA".to_string(); }
        let (a,b,_) = cmap.get(c).copied().unwrap_or((BLACK, WHITE, 0));
        format!("{}/{}", rgb_hex(a), rgb_hex(b))
    }).collect()
}

/// Return the exact continuous colour used by the rooted SVG renderer.
pub fn rooted_continuous_hex(values: &[Option<f32>], root_node: usize) -> Vec<String> {
    let vmax = values.iter().enumerate().filter(|(i,_)| *i != root_node)
        .filter_map(|(_,v)| *v).fold(0.0f32, f32::max).max(1.0);
    values.iter().enumerate().map(|(i,v)| {
        if i == root_node { "#AAAAAA".to_string() }
        else { v.map(|x| rgb_hex(blue_yellow_red(x, 0.0, vmax))).unwrap_or_else(|| "#AAAAAA".to_string()) }
    }).collect()
}

fn radial_positions(layout: &RootedLayout, root_node: usize) -> Vec<(f32, f32)> {
    fn assign_angles(v: usize, children: &[Vec<usize>], angles: &mut [f32], next: &mut f32) {
        if children[v].is_empty() {
            angles[v] = *next;
            *next += 1.0;
        } else {
            for &c in &children[v] {
                assign_angles(c, children, angles, next);
            }
            angles[v] = children[v].iter().map(|&c| angles[c]).sum::<f32>() / children[v].len() as f32;
        }
    }
    let n = layout.parent.len();
    let mut angles = vec![0.0f32; n];
    let mut leaves = 0.0;
    assign_angles(root_node, &layout.children, &mut angles, &mut leaves);
    let denom = leaves.max(1.0);
    let max_depth = layout.depth.iter().copied().max().unwrap_or(0).max(1) as f32;
    let mut out = vec![(0.0, 0.0); n];
    for &v in &layout.order {
        let theta = std::f32::consts::TAU * ((angles[v] + 0.5) / denom) - std::f32::consts::FRAC_PI_2;
        // The observed root starts away from the virtual NAIVE centre; each
        // additional MST depth occupies another ring.
        let radius = 0.18 + 0.78 * (layout.depth[v] as f32 / max_depth);
        out[v] = (radius * theta.cos(), radius * theta.sin());
    }
    out
}

fn layered_positions(
    layout: &RootedLayout,
    root_node: usize,
    abundance: &[usize],
    max_abundance: usize,
) -> (Vec<(f32, f32)>, u32) {
    const NODE_GAP_PX: f32 = 10.0; // deliberately generous test clearance

    let n = layout.parent.len();
    let radii: Vec<f32> = (0..n)
        .map(|v| {
            if v == root_node {
                9.0
            } else {
                abundance_radius(abundance[v], max_abundance) as f32
            }
        })
        .collect();

    // Horizontal columns are MST depths.  Give every depth the radius of its
    // largest node, then separate adjacent depth centres by R[a] + R[b] + gap.
    // This keeps a whole depth column aligned while preventing a large node in
    // one column from colliding with the next column.
    let max_depth = layout.depth.iter().copied().max().unwrap_or(0);
    let mut depth_radius = vec![0.0f32; max_depth + 1];
    for &v in &layout.order {
        depth_radius[layout.depth[v]] = depth_radius[layout.depth[v]].max(radii[v]);
    }
    let mut depth_x = vec![0.0f32; max_depth + 1];
    for d in 1..=max_depth {
        depth_x[d] = depth_x[d - 1] + depth_radius[d - 1] + depth_radius[d] + NODE_GAP_PX;
    }
    let x_total = depth_x.last().copied().unwrap_or(0.0).max(1.0);

    // Vertically, reserve a non-overlapping envelope for every branch.  A
    // branch is at least as tall as its largest node and otherwise as tall as
    // its child envelopes plus the requested circumference-to-circumference
    // gaps.  This is the radius-aware equivalent of the old unit leaf spacing.
    fn envelope(
        v: usize,
        children: &[Vec<usize>],
        radii: &[f32],
        half_height: &mut [f32],
    ) -> f32 {
        if children[v].is_empty() {
            half_height[v] = radii[v];
            return half_height[v];
        }
        let child_height: f32 = children[v]
            .iter()
            .map(|&c| 2.0 * envelope(c, children, radii, half_height))
            .sum::<f32>()
            + NODE_GAP_PX * (children[v].len().saturating_sub(1) as f32);
        half_height[v] = radii[v].max(child_height / 2.0);
        half_height[v]
    }

    fn assign_y(
        v: usize,
        centre: f32,
        children: &[Vec<usize>],
        half_height: &[f32],
        y: &mut [f32],
    ) {
        y[v] = centre;
        if children[v].is_empty() {
            return;
        }
        let total: f32 = children[v]
            .iter()
            .map(|&c| 2.0 * half_height[c])
            .sum::<f32>()
            + NODE_GAP_PX * (children[v].len().saturating_sub(1) as f32);
        let mut cursor = centre - total / 2.0;
        for &c in &children[v] {
            let child_centre = cursor + half_height[c];
            assign_y(c, child_centre, children, half_height, y);
            cursor += 2.0 * half_height[c] + NODE_GAP_PX;
        }
    }

    let mut half_height = vec![0.0f32; n];
    envelope(root_node, &layout.children, &radii, &mut half_height);
    let mut y = vec![0.0f32; n];
    assign_y(root_node, 0.0, &layout.children, &half_height, &mut y);
    let y_extent = half_height[root_node].max(1.0);

    let mut out = vec![(0.0f32, 0.0f32); n];
    for &v in &layout.order {
        let x = -0.92 + 1.84 * depth_x[layout.depth[v]] / x_total;
        let yy = 0.98 * y[v] / y_extent;
        out[v] = (x, yy);
    }
    // The branch envelopes above are in rendered-pixel radii.  A fixed-height
    // canvas would squeeze those coordinates back together during the chart
    // transform and undo the collision calculation. Grow the layered canvas
    // with the required vertical envelope instead.
    let canvas_height = (2.5 * y_extent + 160.0).ceil().max(1000.0) as u32;
    (out, canvas_height)
}

fn draw_split_legend_node<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    centre: (i32, i32),
    r: i32,
    c1: RGBColor,
    c2: RGBColor,
    axis: usize,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    let tilt = split_tilt(axis);
    let half = |offset: f64| {
        let mut pts = vec![centre];
        for step in 0..=18 {
            let a = tilt + offset + std::f64::consts::PI * step as f64 / 18.0;
            pts.push((centre.0 + (r as f64 * a.cos()).round() as i32, centre.1 + (r as f64 * a.sin()).round() as i32));
        }
        pts
    };
    area.draw(&Polygon::new(half(0.0), c1.filled()))?;
    area.draw(&Polygon::new(half(std::f64::consts::PI), c2.filled()))?;
    area.draw(&Circle::new(centre, r, BLACK.stroke_width(1)))?;
    Ok(())
}

fn draw_abundance_legend<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    mut y: i32,
    max_abundance: usize,
) -> Result<i32, DrawingAreaErrorKind<DB::ErrorType>> {
    area.draw(&Text::new("Clone/state size (cells)", (10, y), ("sans-serif", 15).into_font()))?;
    y += 28;
    // Human-readable reference counts rather than an arbitrary sqrt(max)
    // midpoint. Every circle still uses exactly the same radius transform as
    // the plotted nodes.
    let mut examples = vec![1usize];
    for n in [10usize, 50usize] {
        if n <= max_abundance {
            examples.push(n);
        }
    }
    if max_abundance > 1 && !examples.contains(&max_abundance) {
        examples.push(max_abundance);
    }
    examples.sort_unstable();
    examples.dedup();
    for n in examples {
        let r = abundance_radius(n, max_abundance);
        area.draw(&Circle::new((22, y), r, RGBColor(180, 180, 180).filled()))?;
        area.draw(&Circle::new((22, y), r, BLACK.stroke_width(1)))?;
        use plotters::style::text_anchor::{HPos, Pos, VPos};
        let text_style = ("sans-serif", 15)
            .into_text_style(area)
            .pos(Pos::new(HPos::Left, VPos::Center));
        area.draw(&Text::new(format!("{n} cell{}", if n == 1 { "" } else { "s" }), (44, y), text_style))?;
        y += (2 * r + 12).max(26);
    }
    Ok(y + 6)
}

impl MstTree {
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// Write MST edges as TSV: parent<TAB>child<TAB>distance
    pub fn to_tsv<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        self.to_delimited(path, '\t')
    }

    /// Write MST edges using a custom delimiter.
    pub fn to_delimited<P: AsRef<Path>>(&self, path: P, sep: char) -> std::io::Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};

        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        for (p, c, d) in &self.edges {
            writeln!(w, "{}{}{}{}.{:.6}", p, sep, c, sep, d)?;
        }

        Ok(())
    }

    pub fn total_length(&self) -> f32 {
        self.edges.iter().map(|(_, _, d)| d).sum()
    }

    /// Plot an already-built MST as a rooted tree with a categorical overlay.
    ///
    /// The default layered layout puts `HC NAIVE` at the left and successive MST depths
    /// to the right. A radial layout remains available as an explicit option. Large categorical
    /// repertoires use deterministic two-colour split nodes so the number of
    /// distinguishable identities is not limited to the base palette size.
    pub fn plot_rooted_annotated_cached(
        &self,
        n_nodes: usize,
        root_node: usize,
        categories: &[String],
        category_order: Option<&[String]>,
        mixed: &[bool],
        abundance: &[usize],
        annotation_name: &str,
        title: &str,
        outfile: &str,
        radial_layout: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        use plotters::prelude::*;
        use std::collections::{BTreeMap, BTreeSet};
        if n_nodes == 0 {
            return Ok(());
        }
        if categories.len() != n_nodes || mixed.len() != n_nodes || abundance.len() != n_nodes {
            return Err("rooted annotation vectors are not row-aligned with ClonoMap states".into());
        }

        let rooted = self.rooted_layout(n_nodes, root_node);
        let max_abundance = abundance.iter().copied().max().unwrap_or(1).max(1);
        let (positions, canvas_height) = if radial_layout {
            (radial_positions(&rooted, root_node), 1000u32)
        } else {
            layered_positions(&rooted, root_node, abundance, max_abundance)
        };
        let cats: BTreeSet<String> = categories
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != root_node)
            .map(|(_, c)| c.clone())
            .collect();
        let palette = categorical_palette();
        let category_abundance: BTreeMap<String, usize> = categories
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != root_node)
            .fold(BTreeMap::new(), |mut acc, (_, c)| {
                *acc.entry(c.clone()).or_default() += 1;
                acc
            });
        // Categories can optionally arrive in a biologically meaningful order.
        // Valkyrn uses this for LC: neighbouring visual codes then represent
        // neighbouring LC CDR3 consensus sequences rather than lexical IDs.
        let mut ordered_cats = Vec::with_capacity(cats.len());
        if let Some(preferred) = category_order {
            for c in preferred {
                if cats.contains(c) && !ordered_cats.contains(c) {
                    ordered_cats.push(c.clone());
                }
            }
        }
        for c in &cats {
            if !ordered_cats.contains(c) {
                ordered_cats.push(c.clone());
            }
        }
        let cmap: BTreeMap<String, (RGBColor, RGBColor, usize)> = ordered_cats
            .into_iter()
            .enumerate()
            .map(|(i, c)| {
                if c == "unpaired" || c == "HC NAIVE" {
                    (c, (RGBColor(170, 170, 170), RGBColor(170, 170, 170), 0))
                } else {
                    (c, categorical_code(i, &palette))
                }
            })
            .collect();
        let category_count = cmap.len();
        let legend_columns = category_count.div_ceil(25).max(1);
        let legend_width = (legend_columns as u32 * 220).max(280);
        let root = SVGBackend::new(outfile, (1220 + legend_width, canvas_height)).into_drawing_area();
        root.fill(&WHITE)?;
        let (plot, legend) = root.split_horizontally(1220);
        let mut chart = ChartBuilder::on(&plot)
            .caption(title, ("sans-serif", 25))
            .margin(35)
            .build_cartesian_2d(-1.12f32..1.12f32, -1.12f32..1.12f32)?;
        chart.configure_mesh().disable_mesh().x_labels(0).y_labels(0).draw()?;

        for &v in &rooted.order {
            if v != root_node {
                let p = rooted.parent[v];
                chart.draw_series([PathElement::new(vec![positions[p], positions[v]], BLACK.mix(0.38))])?;
            }
        }
        for &v in &rooted.order {
            let r = abundance_radius(abundance[v], max_abundance);
            if v == root_node {
                chart.draw_series([Circle::new(positions[v], 9, RGBColor(170, 170, 170).filled())])?;
                chart.draw_series([Circle::new(positions[v], 9, BLACK.stroke_width(2))])?;
                chart.draw_series([Text::new(
                    "HC NAIVE",
                    (positions[v].0 + 0.025, positions[v].1 + 0.035),
                    ("sans-serif", 18).into_font(),
                )])?;
                continue;
            }
            let (c1, c2, axis) = cmap.get(&categories[v]).copied().unwrap_or((BLACK, WHITE, 0));
            chart.draw_series(PointSeries::of_element(
                [positions[v]], r, &c1, &move |coord, size, _| {
                    let mut a = Vec::with_capacity(20);
                    let mut b = Vec::with_capacity(20);
                    // Diameter tilted 40 degrees. Each polygon includes the centre
                    // and one half of the circumference.
                    let tilt = split_tilt(axis);
                    a.push((0, 0));
                    b.push((0, 0));
                    for step in 0..=18 {
                        let ang = tilt + std::f64::consts::PI * step as f64 / 18.0;
                        a.push(((size as f64 * ang.cos()).round() as i32, (size as f64 * ang.sin()).round() as i32));
                    }
                    for step in 0..=18 {
                        let ang = tilt + std::f64::consts::PI + std::f64::consts::PI * step as f64 / 18.0;
                        b.push(((size as f64 * ang.cos()).round() as i32, (size as f64 * ang.sin()).round() as i32));
                    }
                    EmptyElement::at(coord)
                        + Polygon::new(a, c1.filled())
                        + Polygon::new(b, c2.filled())
                        + Circle::new((0, 0), size, BLACK.stroke_width(2))
                },
            ))?;
            if mixed[v] {
                chart.draw_series(PointSeries::of_element(
                    [positions[v]], r + 2, &BLACK, &|coord, size, style| {
                        EmptyElement::at(coord) + Circle::new((0, 0), size, style.stroke_width(1))
                    },
                ))?;
            }
        }

        legend.draw(&Text::new(annotation_name, (10, 35), ("sans-serif", 20).into_font()))?;
        let rows_per_column = category_count.div_ceil(legend_columns).max(1);
        for (i, (cat, (c1, c2, axis))) in cmap.iter().enumerate() {
            let col = i / rows_per_column;
            let row = i % rows_per_column;
            let x = 10 + col as i32 * 220;
            let y = 67 + row as i32 * 24;
            draw_split_legend_node(&legend, (x + 8, y), 7, *c1, *c2, *axis)?;
            use plotters::style::text_anchor::{HPos, Pos, VPos};
            let text_style = ("sans-serif", 15)
                .into_text_style(&legend)
                .pos(Pos::new(HPos::Left, VPos::Center));
            let n_spheres = category_abundance.get(cat).copied().unwrap_or(0);
            legend.draw(&Text::new(
                format!("{}  [{}]", cat, n_spheres),
                (x + 24, y),
                text_style,
            ))?;
        }
        let mut ly = 67 + rows_per_column as i32 * 24 + 18;
        ly = draw_abundance_legend(&legend, ly, max_abundance)?;
        legend.draw(&Text::new("Black ring = mixed state", (10, ly), ("sans-serif", 15).into_font()))?;
        ly += 24;
        legend.draw(&Text::new("Grey HC NAIVE = inferred unmutated HC; LC empty", (10, ly), ("sans-serif", 15).into_font()))?;
        ly += 24;
        let max_path = rooted.path_distance.iter().copied().fold(0.0f32, f32::max);
        legend.draw(&Text::new(format!("Max MST path: {:.1}", max_path), (10, ly), ("sans-serif", 15).into_font()))?;
        root.present()?;
        Ok(())
    }

    /// Plot the cached rooted MST with a continuous per-state annotation.
    /// The topology is identical to the categorical view and HC NAIVE remains the
    /// root anchor.
    pub fn plot_rooted_continuous_cached(
        &self,
        n_nodes: usize,
        root_node: usize,
        values: &[Option<f32>],
        abundance: &[usize],
        annotation_name: &str,
        title: &str,
        outfile: &str,
        radial_layout: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        use plotters::prelude::*;
        if n_nodes == 0 {
            return Ok(());
        }
        if values.len() != n_nodes || abundance.len() != n_nodes {
            return Err("rooted continuous annotation vectors are not row-aligned with ClonoMap states".into());
        }
        let rooted = self.rooted_layout(n_nodes, root_node);
        let max_abundance = abundance.iter().copied().max().unwrap_or(1).max(1);
        let (positions, canvas_height) = if radial_layout {
            (radial_positions(&rooted, root_node), 1000u32)
        } else {
            layered_positions(&rooted, root_node, abundance, max_abundance)
        };
        let observed_values: Vec<f32> = values
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != root_node)
            .filter_map(|(_, value)| *value)
            .collect();
        // Mutation/depth heat maps have a meaningful zero. Keep the colour
        // scale absolute (0..observed maximum), rather than stretching the
        // observed minimum to blue. Thus distance 2 always means the same
        // colour within a plot whose maximum is 3, even when no observed node
        // happens to have distance 0.
        let vmin = 0.0f32;
        let vmax = observed_values.iter().copied().fold(0.0f32, f32::max).max(1.0);

        let root = SVGBackend::new(outfile, (1500, canvas_height)).into_drawing_area();
        root.fill(&WHITE)?;
        let (plot, legend) = root.split_horizontally(1220);
        let mut chart = ChartBuilder::on(&plot)
            .caption(title, ("sans-serif", 25))
            .margin(35)
            .build_cartesian_2d(-1.12f32..1.12f32, -1.12f32..1.12f32)?;
        chart.configure_mesh().disable_mesh().x_labels(0).y_labels(0).draw()?;
        for &v in &rooted.order {
            if v != root_node {
                let p = rooted.parent[v];
                chart.draw_series([PathElement::new(vec![positions[p], positions[v]], BLACK.mix(0.38))])?;
            }
        }
        for &v in &rooted.order {
            let r = abundance_radius(abundance[v], max_abundance);
            if v == root_node {
                chart.draw_series([Circle::new(positions[v], 9, RGBColor(170, 170, 170).filled())])?;
                chart.draw_series([Circle::new(positions[v], 9, BLACK.stroke_width(2))])?;
                chart.draw_series([Text::new(
                    "HC NAIVE",
                    (positions[v].0 + 0.025, positions[v].1 + 0.035),
                    ("sans-serif", 18).into_font(),
                )])?;
                continue;
            }
            let style = match values[v] {
                Some(x) => blue_yellow_red(x, vmin, vmax).filled(),
                None => RGBColor(170, 170, 170).filled(),
            };
            chart.draw_series([Circle::new(positions[v], r, style)])?;
            chart.draw_series([Circle::new(positions[v], r, BLACK.stroke_width(2))])?;
        }
        let mut ly = 35i32;
        legend.draw(&Text::new(annotation_name, (10, ly), ("sans-serif", 20).into_font()))?;
        ly += 34;

        // Standard continuous colour bar: every numeric value maps onto the
        // same blue -> yellow -> red interpolation used by the nodes.
        const BAR_X: i32 = 10;
        const BAR_W: i32 = 220;
        const BAR_H: i32 = 16;
        for px in 0..BAR_W {
            let t = px as f32 / (BAR_W - 1) as f32;
            let value = vmin + t * (vmax - vmin);
            legend.draw(&Rectangle::new(
                [(BAR_X + px, ly), (BAR_X + px + 1, ly + BAR_H)],
                blue_yellow_red(value, vmin, vmax).filled(),
            ))?;
        }
        legend.draw(&Rectangle::new(
            [(BAR_X, ly), (BAR_X + BAR_W, ly + BAR_H)],
            BLACK.stroke_width(1),
        ))?;
        ly += BAR_H + 5;
        use plotters::style::text_anchor::{HPos, Pos, VPos};
        let left_style = ("sans-serif", 14)
            .into_text_style(&legend)
            .pos(Pos::new(HPos::Left, VPos::Top));
        let right_style = ("sans-serif", 14)
            .into_text_style(&legend)
            .pos(Pos::new(HPos::Right, VPos::Top));
        legend.draw(&Text::new(format!("{vmin:.0}"), (BAR_X, ly), left_style))?;
        legend.draw(&Text::new(format!("{vmax:.0}"), (BAR_X + BAR_W, ly), right_style))?;
        ly += 28;
        ly = draw_abundance_legend(&legend, ly, max_abundance)?;
        legend.draw(&Text::new("Grey HC NAIVE = inferred unmutated HC; LC empty", (10, ly), ("sans-serif", 15).into_font()))?;
        ly += 24;
        let max_path = rooted.path_distance.iter().copied().fold(0.0f32, f32::max);
        legend.draw(&Text::new(format!("Max MST path: {:.1}", max_path), (10, ly), ("sans-serif", 15).into_font()))?;
        root.present()?;
        Ok(())
    }

    fn rooted_layout(&self, n_nodes: usize, root_node: usize) -> RootedLayout {
        let mut adj = vec![Vec::<(usize, f32)>::new(); n_nodes];
        for &(a, b, d) in &self.edges {
            if a < n_nodes && b < n_nodes {
                adj[a].push((b, d));
                adj[b].push((a, d));
            }
        }
        let mut parent = vec![usize::MAX; n_nodes];
        let mut depth = vec![0usize; n_nodes];
        let mut path_distance = vec![0.0f32; n_nodes];
        let mut order = Vec::with_capacity(n_nodes);
        let mut q = VecDeque::new();
        parent[root_node] = root_node;
        q.push_back(root_node);
        while let Some(v) = q.pop_front() {
            order.push(v);
            for &(u, d) in &adj[v] {
                if parent[u] == usize::MAX {
                    parent[u] = v;
                    depth[u] = depth[v] + 1;
                    path_distance[u] = path_distance[v] + d;
                    q.push_back(u);
                }
            }
        }
        let mut children = vec![Vec::<usize>::new(); n_nodes];
        for &v in &order {
            if v != root_node {
                children[parent[v]].push(v);
            }
        }
        RootedLayout { parent, depth, path_distance, order, children }
    }

    /// Build an MST directly from selected rows of a caller-supplied feature matrix.
    ///
    /// Row indices in the returned edges always refer to the original matrix.
    /// This deliberately bypasses PCA and is intended for biologically coherent
    /// subspaces where dimensionality reduction must not alter topology.
    pub fn build_feature_rows(features: &ndarray::Array2<f32>, rows: &[usize]) -> Self {
        if rows.len() < 2 {
            return Self { edges: Vec::new() };
        }

        let n = rows.len();
        let mut in_tree = vec![false; n];
        let mut dist = vec![f32::INFINITY; n];
        let mut parent = vec![None; n];
        in_tree[0] = true;

        for i in 1..n {
            dist[i] = Self::euclidean(features.row(rows[0]), features.row(rows[i]));
            parent[i] = Some(0);
        }

        for _ in 1..n {
            let mut best = None;
            let mut best_d = f32::INFINITY;
            for i in 0..n {
                if !in_tree[i] && dist[i] < best_d {
                    best = Some(i);
                    best_d = dist[i];
                }
            }
            let Some(v) = best else { break };
            in_tree[v] = true;
            for u in 0..n {
                if in_tree[u] {
                    continue;
                }
                let d = Self::euclidean(features.row(rows[v]), features.row(rows[u]));
                if d < dist[u] {
                    dist[u] = d;
                    parent[u] = Some(v);
                }
            }
        }

        let mut edges = Vec::with_capacity(n - 1);
        for i in 1..n {
            if let Some(p) = parent[i] {
                edges.push((rows[p], rows[i], dist[i]));
            }
        }
        Self { edges }
    }

    pub fn build(model: &PcaModel) -> Self {
        use rayon::prelude::*;
        let coords = &model.coords;

        let n = coords.nrows();
        if n == 0 {
            return Self { edges: Vec::new() };
        }

        // MST state
        let mut in_tree = vec![false; n];
        let mut dist = vec![f32::INFINITY; n];
        let mut parent = vec![None; n];

        in_tree[0] = true;

        // Initialize distances relative to node 0
        for i in 1..n {
            dist[i] = Self::euclidean(coords.row(0), coords.row(i));
            parent[i] = Some(0);
        }

        // --- main MST loop ---
        for _ in 1..(n - 1) {
            // Find best vertex outside MST (sequential O(n))
            let mut best = None;
            let mut best_d = f32::INFINITY;

            for i in 0..n {
                if !in_tree[i] && dist[i] < best_d {
                    best = Some(i);
                    best_d = dist[i];
                }
            }

            let v = best.expect("MST cannot proceed – no available node");
            in_tree[v] = true;

            let v_row = coords.row(v);

            // --- PARALLEL relax edges ---
            // Instead of updating dist[u] directly (data race)
            // we compute all improvements first, then apply them.
            let updates: Vec<(usize, f32, usize)> = (0..n)
                .into_par_iter()
                .filter_map(|u| {
                    if in_tree[u] {
                        return None;
                    }

                    let d = Self::euclidean(v_row, coords.row(u));
                    if d < dist[u] {
                        Some((u, d, v)) // (node, new-dist, parent)
                    } else {
                        None
                    }
                })
                .collect();

            // Apply computed updates (sequential, very cheap)
            for (u, newd, pv) in updates {
                dist[u] = newd;
                parent[u] = Some(pv);
            }
        }

        // Collect MST edges
        let mut edges = Vec::with_capacity(n - 1);
        for i in 1..n {
            edges.push((parent[i].unwrap(), i, dist[i]));
        }

        Self { edges }
    }

    /// Find sparse root: the node located in the lowest-density PCA region.
    /// Uses the distance to the k-th nearest neighbor.
    pub fn find_sparse_root(model: &PcaModel, k: usize) -> usize {
        let coords = &model.coords;
        let n = coords.nrows();
        let mut scores = vec![0.0; n];

        for i in 0..n {
            // Distance to all other nodes
            let mut dists: Vec<f32> = (0..n)
                .filter(|&j| j != i)
                .map(|j| {
                    coords
                        .row(i)
                        .iter()
                        .zip(coords.row(j).iter())
                        .map(|(x, y)| (x - y).powi(2))
                        .sum::<f32>()
                        .sqrt()
                })
                .collect();

            // Sort and store k-th smallest distance
            dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
            scores[i] = dists[k.min(dists.len() - 1)];
        }

        // Node with the largest kNN radius = sparsest region = predicted root
        scores
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0
    }

    pub fn clusters_elbow(&self, n_nodes: usize) -> Vec<Vec<usize>> {
        let Some(threshold) = self.elbow_threshold() else {
            return vec![];
        };

        self.clusters_with_cut(n_nodes, threshold)
    }

    pub fn clusters_robust(&self, n_nodes: usize) -> Vec<Vec<usize>> {
        let Some(threshold) = self.robust_threshold_auto() else {
            return vec![];
        };

        self.clusters_with_cut(n_nodes, threshold)
    }

    pub fn clusters_with_cut(&self, n_nodes: usize, max_len: f32) -> Vec<Vec<usize>> {
        let mut adj = vec![Vec::new(); n_nodes];

        for (a, b, d) in &self.edges {
            if *d <= max_len {
                adj[*a].push(*b);
                adj[*b].push(*a);
            }
        }

        let mut visited = vec![false; n_nodes];
        let mut out = Vec::new();

        for i in 0..n_nodes {
            if visited[i] {
                continue;
            }

            let mut stack = VecDeque::new();
            let mut comp = Vec::new();

            stack.push_back(i);
            visited[i] = true;

            while let Some(u) = stack.pop_front() {
                comp.push(u);
                for &v in &adj[u] {
                    if !visited[v] {
                        visited[v] = true;
                        stack.push_back(v);
                    }
                }
            }

            out.push(comp);
        }

        out
    }
    /// Automatically chooses clustering threshold using elbow detection.
    pub fn elbow_threshold(&self) -> Option<f32> {
        if self.edges.len() < 2 {
            return None;
        }

        let mut lens: Vec<f32> = self.edges.iter().map(|(_, _, d)| *d).collect();
        lens.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut best_i = 0;
        let mut best_gap = 0.0;

        for i in 0..lens.len() - 1 {
            let gap = lens[i + 1] - lens[i];
            if gap > best_gap {
                best_gap = gap;
                best_i = i;
            }
        }

        Some(lens[best_i])
    }

    pub fn robust_threshold(&self, k: f32) -> Option<f32> {
        if self.edges.len() < 2 {
            return None;
        }

        let mut x: Vec<f32> = self.edges.iter().map(|(_, _, d)| *d).collect();
        x.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let med = median(&x);

        let dev: Vec<f32> = x.iter().map(|v| (v - med).abs()).collect();
        let mad = median(&dev);

        Some(med + k * mad)
    }

    pub fn robust_threshold_auto(&self) -> Option<f32> {
        if self.edges.len() < 4 {
            return None;
        }

        let mut x: Vec<f32> = self.edges.iter().map(|(_, _, d)| *d).collect();
        x.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // --- median ---
        let med = median(&x);

        // --- MAD ---
        let dev: Vec<f32> = x.iter().map(|v| (v - med).abs()).collect();
        let mad = median(&dev).max(1e-9);

        // --- normalized tail weights ---
        // z-score-like: (x - median) / MAD
        let z: Vec<f32> = x.iter().map(|v| (v - med) / mad).collect();

        // --- detect first big tail rise ---
        // find first value beyond a natural outlier region
        let mut cut = None;

        for i in 0..z.len() {
            // "unlikely under normal" threshold
            if z[i] > 3.5 && i > x.len() / 2 {
                cut = Some(x[i]);
                break;
            }
        }

        // --- fallback: percentile based ---
        if cut.is_none() {
            let idx = ((x.len() as f32) * 0.85) as usize;
            cut = Some(x[idx.min(x.len() - 1)]);
        }

        cut
    }

    pub fn cut(&self, max_len: f32) -> Vec<(usize, usize)> {
        self.edges
            .iter()
            .filter(|(_, _, d)| *d <= max_len)
            .map(|(a, b, _)| (*a, *b))
            .collect()
    }

        pub fn plot_2d(
        &self,
        coords: &ndarray::Array2<f32>,
        outfile: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        use plotters::prelude::*;

        let root = SVGBackend::new(outfile, (900, 900)).into_drawing_area();
        root.fill(&WHITE)?;

        let x = coords.column(0);
        let y = coords.column(1);

        let xmin = x.iter().cloned().fold(f32::INFINITY, f32::min);
        let xmax = x.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let ymin = y.iter().cloned().fold(f32::INFINITY, f32::min);
        let ymax = y.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        let mut chart = ChartBuilder::on(&root)
            .caption("PCA Tree", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

        chart.configure_mesh().draw()?;

        // Draw edges (lines)
        for &(a, b, _) in &self.edges {
            let pa = (coords[(a, 0)], coords[(a, 1)]);
            let pb = (coords[(b, 0)], coords[(b, 1)]);
            chart.draw_series([PathElement::new(vec![pa, pb], &BLACK)])?;
        }

        // Draw nodes
        chart.draw_series(
            x.iter()
                .zip(y.iter())
                .map(|(&x, &y)| Circle::new((x, y), 3, RED.filled())),
        )?;

        Ok(())
    }

    /// Re-root the MST at the given node index
    pub fn reroot(&self, n: usize, root: usize) -> MstTree {
        let mut adj = vec![Vec::new(); n];

        // Build adjacency
        for &(p, c, d) in &self.edges {
            adj[p].push((c, d));
            adj[c].push((p, d));
        }

        let mut parent = vec![None; n];
        let mut dist_to_parent = vec![0.0f32; n];

        let mut queue = VecDeque::new();
        queue.push_back(root);
        parent[root] = Some(root); // mark root

        // BFS to orient edges
        while let Some(v) = queue.pop_front() {
            for &(nbr, d) in &adj[v] {
                if parent[nbr].is_none() {
                    parent[nbr] = Some(v);
                    dist_to_parent[nbr] = d;
                    queue.push_back(nbr);
                }
            }
        }

        // Rebuild edges (skip the root)
        let mut new_edges = Vec::new();
        for i in 0..n {
            if i == root {
                continue;
            }
            let p = parent[i].unwrap();
            new_edges.push((p, i, dist_to_parent[i]));
        }

        MstTree { edges: new_edges }
    }

    pub fn to_newick(&self, n: usize, root: usize, labels: &CloneData) -> String {
        let mut children = vec![Vec::new(); n];
        for &(p, c, d) in &self.edges {
            children[p].push((c, d));
        }

        fn build(idx: usize, children: &Vec<Vec<(usize, f32)>>, aa_labels: &[String]) -> String {
            if children[idx].is_empty() {
                format!("{}", aa_labels[idx])
            } else {
                let inner: Vec<String> = children[idx]
                    .iter()
                    .map(|(c, d)| format!("{}:{:.4}", build(*c, children, aa_labels), d))
                    .collect();

                format!("({})", inner.join(","))
            }
        }
        let keys = labels.aa_with_count_labels();

        format!("{};", build(root, &children, &keys))
    }

    fn euclidean(a: ArrayView1<f32>, b: ArrayView1<f32>) -> f32 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f32>()
            .sqrt()
    }
}

fn median(v: &[f32]) -> f32 {
    let m = v.len() / 2;
    if v.len() % 2 == 0 {
        (v[m - 1] + v[m]) / 2.0
    } else {
        v[m]
    }
}


#[cfg(test)]
mod plot_encoding_tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn split_palette_supports_one_hundred_twenty_distinct_axis_codes() {
        let palette = categorical_palette();
        let codes: BTreeSet<_> = (0..120)
            .map(|i| {
                let (a, b, axis) = categorical_code(i, &palette);
                (a.0, a.1, a.2, b.0, b.1, b.2, axis)
            })
            .collect();
        assert_eq!(codes.len(), 120);
        assert!((0..120).all(|i| {
            let (a, b, _) = categorical_code(i, &palette);
            a != b
        }));
    }

    #[test]
    fn abundance_legend_uses_exact_node_radius_transform() {
        assert_eq!(abundance_radius(1, 100), 6);
        assert_eq!(abundance_radius(100, 100), 15);
        assert!(abundance_radius(25, 100) > abundance_radius(1, 100));
    }

    #[test]
    fn radial_layout_keeps_virtual_naive_centre_free() {
        let layout = RootedLayout {
            parent: vec![0, 0, 0],
            depth: vec![0, 1, 1],
            path_distance: vec![0.0, 1.0, 1.0],
            order: vec![0, 1, 2],
            children: vec![vec![1, 2], vec![], vec![]],
        };
        for (x, y) in radial_positions(&layout, 0) {
            assert!((x * x + y * y).sqrt() >= 0.18 - f32::EPSILON);
        }
    }
}
