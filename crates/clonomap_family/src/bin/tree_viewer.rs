//tree_viewer.rs
use plotters::prelude::*;
use std::env;
use std::error::Error;
use std::fs;
use std::iter::Peekable;
use std::str::Chars;

/// Minimal Newick node structure
#[derive(Debug, Clone)]
struct SimpleNode {
    name: String,
    length: f64,
    children: Vec<SimpleNode>,
}

/// Entry point
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: tree_viewer <tree.newick> [output.svg]");
        std::process::exit(1);
    }

    let in_path = &args[1];
    let out_path = if args.len() >= 3 {
        args[2].clone()
    } else {
        format!("{in_path}.svg")
    };

    let newick_text = fs::read_to_string(in_path)?;
    let root = parse_newick(&newick_text);

    // Layout tree into 2D coordinates
    let mut points: Vec<(f64, f64)> = Vec::new();
    let mut labels: Vec<String> = Vec::new();
    let mut edges: Vec<(usize, usize)> = Vec::new();

    let mut next_y = 0.0;
    layout_tree(
        &root,
        0.0,  // parent_x
        None, // parent_idx
        &mut next_y,
        &mut points,
        &mut labels,
        &mut edges,
    );

    // Compute bounds for plotting
    let (min_x, max_x, min_y, max_y) = bounds(&points);

    // Draw to SVG
    draw_svg(
        &out_path, &points, &labels, &edges, min_x, max_x, min_y, max_y,
    )?;

    println!("✅ Tree written to {}", out_path);
    Ok(())
}

/// Compute min/max over points
fn bounds(points: &[(f64, f64)]) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;

    for (x, y) in points {
        if *x < min_x {
            min_x = *x;
        }
        if *x > max_x {
            max_x = *x;
        }
        if *y < min_y {
            min_y = *y;
        }
        if *y > max_y {
            max_y = *y;
        }
    }

    // Add small padding so nodes are not at the edge
    let pad_x = (max_x - min_x).abs() * 0.05 + 1.0;
    let pad_y = (max_y - min_y).abs() * 0.05 + 1.0;

    (min_x - pad_x, max_x + pad_x, min_y - pad_y, max_y + pad_y)
}

/// Recursively layout the tree
///
/// - x = parent_x + branch_length
/// - leaves get consecutive y values
/// - internal nodes get mean y of their children
fn layout_tree(
    node: &SimpleNode,
    parent_x: f64,
    parent_idx: Option<usize>,
    next_y: &mut f64,
    points: &mut Vec<(f64, f64)>,
    labels: &mut Vec<String>,
    edges: &mut Vec<(usize, usize)>,
) -> usize {
    let x = parent_x + node.length;
    let idx = points.len();

    // temporary y; will be updated later
    points.push((x, 0.0));
    labels.push(node.name.clone());

    if let Some(p) = parent_idx {
        edges.push((p, idx));
    }

    if node.children.is_empty() {
        // leaf
        let y = *next_y;
        *next_y += 1.0;
        points[idx].1 = y;
    } else {
        let mut child_ys = Vec::new();
        for child in &node.children {
            let cidx = layout_tree(child, x, Some(idx), next_y, points, labels, edges);
            child_ys.push(points[cidx].1);
        }
        let y = child_ys.iter().sum::<f64>() / child_ys.len() as f64;
        points[idx].1 = y;
    }

    idx
}

/// Draw the tree to SVG using plotters
fn draw_svg(
    path: &str,
    points: &[(f64, f64)],
    labels: &[String],
    edges: &[(usize, usize)],
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
) -> Result<(), Box<dyn Error>> {
    let root = SVGBackend::new(path, (2000, 1200)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Phylogenetic Tree", ("sans-serif", 30))
        .margin(20)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart.configure_mesh().draw()?;

    // Draw edges
    chart.draw_series(edges.iter().map(|&(p, c)| {
        let (x1, y1) = points[p];
        let (x2, y2) = points[c];
        PathElement::new(vec![(x1, y1), (x2, y2)], &BLACK)
    }))?;

    // Mark nodes, and label leaves
    // compute which nodes are leaves (no outgoing edges as parent)
    let mut has_children = vec![false; points.len()];
    for &(_, c) in edges {
        has_children[c] = true;
    }

    chart.draw_series(points.iter().enumerate().map(|(i, (x, y))| {
        if !has_children[i] && !labels[i].is_empty() {
            EmptyElement::at((*x, *y))
                + Circle::new((0, 0), 3, BLUE.filled())
                + Text::new(
                    labels[i].clone(),
                    (6, 0),
                    ("sans-serif", 15).into_font().style(FontStyle::Normal),
                )
        } else {
            EmptyElement::at((*x, *y))
                + Circle::new((0, 0), 3, BLUE.filled())
                + Text::new(
                    "NA".to_string(),
                    (6, 0),
                    ("sans-serif", 15).into_font().style(FontStyle::Normal),
                )
        }
    }))?;

    root.present()?;
    Ok(())
}

//
// ----------- Minimal Newick parser ----------------
//

fn parse_newick(input: &str) -> SimpleNode {
    let mut chars = input.trim().chars().peekable();
    let node = parse_node(&mut chars);

    // Consume trailing semicolon if present
    if let Some(';') = chars.peek() {
        chars.next();
    }

    node
}

fn skip_ws(chars: &mut Peekable<Chars<'_>>) {
    while let Some(&c) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
        } else {
            break;
        }
    }
}

fn parse_node(chars: &mut Peekable<Chars<'_>>) -> SimpleNode {
    skip_ws(chars);

    let mut children = Vec::new();

    // Children block: '(' node (,node)* ')'
    if matches!(chars.peek(), Some('(')) {
        chars.next(); // '('
        loop {
            let child = parse_node(chars);
            children.push(child);

            skip_ws(chars);
            match chars.peek() {
                Some(',') => {
                    chars.next();
                    continue;
                }
                Some(')') => {
                    chars.next();
                    break;
                }
                other => {
                    panic!("Invalid Newick: expected ',' or ')', got {:?}", other);
                }
            }
        }
    }

    skip_ws(chars);

    // Label
    let mut name = String::new();
    while let Some(&c) = chars.peek() {
        if c == ':' || c == ',' || c == ')' || c == ';' {
            break;
        }
        if !c.is_whitespace() {
            name.push(c);
        }
        chars.next();
    }

    skip_ws(chars);

    // Branch length
    let mut length = 0.0_f64;
    if matches!(chars.peek(), Some(':')) {
        chars.next(); // ':'
        skip_ws(chars);
        let mut num = String::new();
        while let Some(&c) = chars.peek() {
            if c == ',' || c == ')' || c == ';' {
                break;
            }
            num.push(c);
            chars.next();
        }
        length = num.parse::<f64>().unwrap_or(0.0);
    }

    SimpleNode {
        name,
        length,
        children,
    }
}
