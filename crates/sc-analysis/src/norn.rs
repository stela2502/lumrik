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

/// Accept either the canonical expression MEX itself, a Nelrune output directory,
/// or a Norn sample directory. Prefer the current `exprs` contract while retaining
/// the historical `exonic` layouts for older runs.
pub(crate) fn resolve_exonic(input: &Path) -> Result<PathBuf> {
    let candidates = [
        input.to_path_buf(),
        input.join("exprs"),
        input.join("filtered/exprs"),
        input.join("nelrune_out/filtered/exprs"),
        input.join("nelrune/nelrune_out/filtered/exprs"),
        input.join("exonic"),
        input.join("filtered/exonic"),
        input.join("nelrune_out/exonic"),
        input.join("nelrune/nelrune_out/exonic"),
    ];
    for candidate in candidates {
        if has_mex(&candidate) {
            return Ok(candidate);
        }
    }
    anyhow::bail!(
        "could not find a canonical expression Matrix Market dataset below {}; expected the MEX itself, a current exprs/ layout, or a legacy exonic/ layout",
        input.display()
    )
}

#[derive(Debug, Clone)]
pub(crate) struct ClonomapBlock {
    pub headers: Vec<String>,
    pub columns: Vec<Vec<String>>,
}

/// Find ClonoMap's authoritative per-cell result from a sample/Norn output root.
/// ClonoMap is optional: samples without receptor analysis simply return None.
pub(crate) fn resolve_clonomap_cells(input: &Path) -> Option<PathBuf> {
    let candidates = [
        input.join("clonomap_out/cells.tsv"),
        input.join("clonomap/cells.tsv"),
        input.join("nelrune/clonomap_out/cells.tsv"),
    ];
    candidates.into_iter().find(|p| p.is_file())
}

pub(crate) fn load_clonomap(input: &Path, cells: &[String]) -> Result<Option<ClonomapBlock>> {
    let Some(path) = resolve_clonomap_cells(input) else {
        return Ok(None);
    };
    let cell_index = cells
        .iter()
        .enumerate()
        .map(|(i, x)| (x.as_str(), i))
        .collect::<HashMap<_, _>>();
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(&path)
        .with_context(|| format!("reading ClonoMap cells {}", path.display()))?;
    let raw_headers = rdr.headers()?.clone();
    let cell_col = raw_headers
        .iter()
        .position(|x| x == "cell")
        .with_context(|| format!("ClonoMap file {} lacks column cell", path.display()))?;
    let keep = raw_headers
        .iter()
        .enumerate()
        .filter(|(_, h)| *h != "cell")
        .map(|(i, h)| (i, format!("clonomap_{}", clean_name(h))))
        .collect::<Vec<_>>();
    let headers = keep.iter().map(|(_, h)| h.clone()).collect::<Vec<_>>();
    let mut columns = vec![vec![String::new(); cells.len()]; keep.len()];
    for rec in rdr.records() {
        let rec = rec?;
        let Some(&cell_i) = rec.get(cell_col).and_then(|x| cell_index.get(x)) else {
            continue;
        };
        for (out_i, (src_i, _)) in keep.iter().enumerate() {
            columns[out_i][cell_i] = rec.get(*src_i).unwrap_or("").to_string();
        }
    }
    Ok(Some(ClonomapBlock { headers, columns }))
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

#[cfg(test)]
mod tests {
    use super::{load_clonomap, resolve_clonomap_cells, resolve_exonic};
    use std::fs;

    #[test]
    fn resolves_current_exprs_and_clonomap_from_sample_root() {
        let root = std::env::temp_dir().join(format!("sc_analysis_norn_contract_{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        let exprs = root.join("filtered/exprs");
        fs::create_dir_all(&exprs).unwrap();
        for name in ["matrix.mtx", "features.tsv", "barcodes.tsv"] {
            fs::write(exprs.join(name), "").unwrap();
        }
        let clonomap = root.join("clonomap_out");
        fs::create_dir_all(&clonomap).unwrap();
        fs::write(
            clonomap.join("cells.tsv"),
            "source\tcell\tfamily\tlc_clone\thc_mutation_count\nS1\tCELL_A\tHC_1\tLC_2\t7\n",
        ).unwrap();

        assert_eq!(resolve_exonic(&root).unwrap(), exprs);
        assert_eq!(resolve_clonomap_cells(&root).unwrap(), clonomap.join("cells.tsv"));
        let block = load_clonomap(&root, &["CELL_A".into(), "CELL_B".into()]).unwrap().unwrap();
        let family = block.headers.iter().position(|x| x == "clonomap_family").unwrap();
        assert_eq!(block.columns[family], ["HC_1", ""]);
        let lc = block.headers.iter().position(|x| x == "clonomap_lc_clone").unwrap();
        assert_eq!(block.columns[lc], ["LC_2", ""]);
        fs::remove_dir_all(&root).unwrap();
    }
}
