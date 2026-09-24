use anyhow::{Context, Result};
use csv::ReaderBuilder;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub(crate) struct BeaconBlock {
    pub name: String,
    pub assignment: Vec<String>,
    pub called_features: Vec<String>,
    pub label: Vec<String>,
    pub n_called: Vec<String>,
    pub best_feature: Vec<String>,
    pub best_log_odds: Vec<String>,
}

fn has_mex(path: &Path) -> bool {
    ["matrix.mtx.gz", "matrix.mtx"]
        .iter()
        .any(|x| path.join(x).is_file())
        && ["features.tsv.gz", "features.tsv"]
            .iter()
            .any(|x| path.join(x).is_file())
        && ["barcodes.tsv.gz", "barcodes.tsv"]
            .iter()
            .any(|x| path.join(x).is_file())
}

/// Accept either the exonic MEX itself, a Nelrune output directory, or the
/// Norn sample directory that contains `nelrune/nelrune_out/exonic`.
pub(crate) fn resolve_exonic(input: &Path) -> Result<PathBuf> {
    let candidates = [
        input.to_path_buf(),
        input.join("exonic"),
        input.join("nelrune_out/exonic"),
        input.join("nelrune/nelrune_out/exonic"),
    ];
    for candidate in candidates {
        if has_mex(&candidate) {
            return Ok(candidate);
        }
    }
    anyhow::bail!(
        "could not find an exonic Matrix Market dataset below {}; expected the MEX itself or one of exonic/, nelrune_out/exonic/, nelrune/nelrune_out/exonic/",
        input.display()
    )
}

fn clean_name(s: &str) -> String {
    let x = s
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>();
    x.trim_matches('_').to_ascii_lowercase()
}

pub(crate) fn load_beacon_blocks(exonic: &Path, cells: &[String]) -> Result<Vec<BeaconBlock>> {
    let Some(base) = exonic.parent() else {
        return Ok(Vec::new());
    };
    let cell_index = cells
        .iter()
        .enumerate()
        .map(|(i, x)| (x.as_str(), i))
        .collect::<HashMap<_, _>>();
    let mut paths = Vec::new();
    for entry in fs::read_dir(base).with_context(|| format!("reading {}", base.display()))? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            let p = entry.path().join("beacon/cell_guide_assignments.tsv");
            if p.is_file() {
                paths.push((entry.file_name().to_string_lossy().to_string(), p));
            }
        }
    }
    paths.sort_by(|a, b| a.0.cmp(&b.0));
    let mut blocks = Vec::new();
    for (feature_type, path) in paths {
        let mut assignment = vec!["unassigned".to_string(); cells.len()];
        let mut n_called = vec!["0".to_string(); cells.len()];
        let mut called_features = vec![String::new(); cells.len()];
        let mut label = vec!["none".to_string(); cells.len()];
        let mut best_feature = vec![String::new(); cells.len()];
        let mut best_log_odds = vec![String::new(); cells.len()];
        let mut rdr = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&path)
            .with_context(|| format!("reading Beacon assignments {}", path.display()))?;
        let headers = rdr.headers()?.clone();
        let col = |name: &str| {
            headers
                .iter()
                .position(|x| x == name)
                .with_context(|| format!("Beacon file {} lacks column {name}", path.display()))
        };
        let barcode_i = col("barcode")?;
        let n_i = col("n_called_guides")?;
        let assignment_i = col("assignment")?;
        let called_i = col("called_guides")?;
        let best_i = col("best_guide")?;
        let odds_i = col("best_log_odds")?;
        for rec in rdr.records() {
            let rec = rec?;
            let Some(&i) = rec.get(barcode_i).and_then(|x| cell_index.get(x)) else {
                continue;
            };
            assignment[i] = rec.get(assignment_i).unwrap_or("unassigned").to_string();
            n_called[i] = rec.get(n_i).unwrap_or("0").to_string();
            called_features[i] = rec.get(called_i).unwrap_or("").to_string();
            label[i] = match assignment[i].as_str() {
                "single" => called_features[i].clone(),
                "multi" => format!("multi:{}", called_features[i]),
                _ => "none".to_string(),
            };
            best_feature[i] = rec.get(best_i).unwrap_or("").to_string();
            best_log_odds[i] = rec.get(odds_i).unwrap_or("").to_string();
        }
        blocks.push(BeaconBlock {
            name: clean_name(&feature_type),
            assignment,
            called_features,
            label,
            n_called,
            best_feature,
            best_log_odds,
        });
    }
    Ok(blocks)
}
