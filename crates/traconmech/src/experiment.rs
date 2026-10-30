//! Resolve experimental meaning once; keep scientist decisions separate from
//! deterministic per-population reference matching.
use anyhow::{Context, Result, bail};
use sc_analysis::SingleCellData;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{self, IsTerminal, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Design {
    pub cases: Vec<String>,
    pub controls: Vec<String>,
}
#[derive(Debug, Serialize, Deserialize)]
struct Experiment {
    version: u32,
    population_column: String,
    condition_column: String,
    observed_conditions: BTreeMap<String, usize>,
    reference_conditions: Vec<String>,
    comparison_conditions: Vec<String>,
    #[serde(default)]
    decision_source: String,
    min_reference_cells: usize,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    population_references: Vec<Reference>,
}
#[derive(Debug, Serialize, Deserialize)]
struct Reference {
    population: String,
    case_cells: usize,
    own_reference_cells: usize,
    selected_population: Option<String>,
    selected_reference_cells: usize,
    method: String,
    pearson_log1p_mean: Option<f64>,
}
fn labels(raw: &str) -> Vec<String> {
    let mut v: Vec<String> = raw
        .split(',')
        .map(str::trim)
        .filter(|x| !x.is_empty())
        .map(str::to_owned)
        .collect();
    v.sort();
    v.dedup();
    v
}
fn reference_name(s: &str) -> bool {
    let lower = s.to_ascii_lowercase();
    [
        "wt",
        "wildtype",
        "wild-type",
        "control",
        "ctrl",
        "untreated",
        "vehicle",
        "mock",
    ]
    .iter()
    .any(|prefix| {
        lower == *prefix
            || lower.starts_with(&format!("{prefix}_"))
            || lower.starts_with(&format!("{prefix}-"))
            || (lower.starts_with(prefix)
                && lower[prefix.len()..]
                    .chars()
                    .next()
                    .is_some_and(|c| c.is_ascii_digit()))
    })
}
fn ask(prompt: &str) -> Result<Vec<String>> {
    print!("{prompt}: ");
    io::stdout().flush()?;
    let mut answer = String::new();
    anyhow::ensure!(io::stdin().read_line(&mut answer)? != 0, "input closed while selecting labels");
    Ok(labels(&answer))
}
// Unknown responses are not decisions: ask again, without losing prior work.
fn confirm(prompt: &str, default_yes: bool) -> Result<bool> {
    loop {
        print!("{prompt} ");
        io::stdout().flush()?;
        let mut answer = String::new();
        anyhow::ensure!(io::stdin().read_line(&mut answer)? != 0, "input closed during confirmation");
        match parse_confirmation(&answer, default_yes) {
            Some(accepted) => return Ok(accepted),
            None => eprintln!("Please answer Y or N."),
        }
    }
}

fn parse_confirmation(answer: &str, default_yes: bool) -> Option<bool> {
    match answer.trim().to_ascii_lowercase().as_str() {
        "y" | "yes" => Some(true),
        "n" | "no" => Some(false),
        "" => Some(default_yes),
        _ => None,
    }
}

/// Report every loaded /obs header when a requested annotation is unavailable.
/// This also exposes undecoded H5AD columns, which may not be selectable.
pub(crate) fn require_annotation<'a>(
    data: &'a SingleCellData,
    name: &str,
    role: &str,
) -> Result<&'a [String]> {
    data.annotations.column(name).with_context(|| {
        let headers = data.annotations.headers();
        let available = if headers.is_empty() {
            "  (none)".to_owned()
        } else {
            headers
                .iter()
                .map(|header| format!("  {header}"))
                .collect::<Vec<_>>()
                .join("\n")
        };
        format!("missing {role} annotation {name}\nAvailable H5AD /obs columns:\n{available}")
    })
}

pub fn resolve(
    data: &SingleCellData,
    out: &Path,
    population: &str,
    condition: &str,
    cases: Option<&str>,
    controls: Option<&str>,
    min_reference_cells: usize,
    non_interactive: bool,
) -> Result<Design> {
    anyhow::ensure!(min_reference_cells > 0, "--min-wt-cells must be positive");
    require_annotation(data, population, "population")?;
    let column = require_annotation(data, condition, "condition")?;
    let mut observed = BTreeMap::new();
    for value in column {
        *observed.entry(value.clone()).or_insert(0usize) += 1;
    }
    let file = out.join("experiment.yaml");
    let old: Option<Experiment> = if file.exists() {
        Some(
            serde_yaml::from_str(&fs::read_to_string(&file)?)
                .with_context(|| format!("invalid {}", file.display()))?,
        )
    } else {
        None
    };
    if let Some(previous) = &old {
        anyhow::ensure!(
            previous.version == 1,
            "unsupported experiment.yaml version in {} (found {}, expected 1)",
            file.display(),
            previous.version
        );
        anyhow::ensure!(
            previous.condition_column == condition && previous.population_column == population,
            "experiment metadata columns changed in {}: saved population={:?}, condition={:?}; current population={:?}, condition={:?}. Review this file before rerunning",
            file.display(),
            previous.population_column,
            previous.condition_column,
            population,
            condition
        );
        anyhow::ensure!(
            previous.observed_conditions == observed,
            "experiment sample labels or cell counts changed in {}\nSaved: {:?}\nCurrent: {:?}\nReview this file before rerunning",
            file.display(),
            previous.observed_conditions,
            observed
        );
    }
    // Explicit CLI choices override saved selections. Otherwise reuse the
    // scientist-approved experiment.yaml; prompt only for a new experiment.
    let interactive = !non_interactive && io::stdin().is_terminal();
    let mut correction = false;
    let (references, comparisons, source) = loop {
        let (references, comparisons, source) = if correction {
            eprintln!("Observed condition labels (cell counts): {observed:?}");
            (
                ask("Reference labels (comma-separated, exact)")?,
                ask("Comparison labels (comma-separated, exact)")?,
                "scientist".to_owned(),
            )
        } else if cases.is_some() || controls.is_some() {
            // An omitted side is intentionally empty, not a malformed CLI.
            (
                labels(controls.unwrap_or("")),
                labels(cases.unwrap_or("")),
                "cli".to_owned(),
            )
        } else if let Some(previous) = &old {
            (
                previous.reference_conditions.clone(),
                previous.comparison_conditions.clone(),
                previous.decision_source.clone(),
            )
        } else if interactive {
            eprintln!("Observed condition labels (cell counts): {observed:?}");
            (
                ask("Reference labels (comma-separated, exact)")?,
                ask("Comparison labels (comma-separated, exact)")?,
                "scientist".to_owned(),
            )
        } else {
            eprintln!("[traconmech] WARNING: no experimental groups supplied; proceeding with empty groups");
            (Vec::new(), Vec::new(), "batch-unspecified".to_owned())
        };
        let ref_set: BTreeSet<_> = references.iter().collect();
        let case_set: BTreeSet<_> = comparisons.iter().collect();
        let overlap: Vec<&str> = ref_set.intersection(&case_set).map(|s| s.as_str()).collect();
        if !overlap.is_empty() {
            eprintln!("[traconmech] ERROR: reference and comparison labels overlap: {}", overlap.join(", "));
            if interactive {
                correction = true;
                continue;
            }
            bail!("reference and comparison labels must be disjoint");
        }
        let unknown: Vec<&str> = ref_set.union(&case_set)
            .filter(|name| !observed.contains_key(**name))
            .map(|s| s.as_str()).collect();
        if !unknown.is_empty() {
            eprintln!("[traconmech] WARNING: supplied labels absent from dataset (ignored): {}", unknown.join(", "));
        }
        let unassigned: Vec<&str> = observed.keys()
            .filter(|name| !ref_set.contains(name) && !case_set.contains(name))
            .map(String::as_str).collect();
        if !unassigned.is_empty() {
            eprintln!("[traconmech] WARNING: dataset labels not assigned to either group: {}", unassigned.join(", "));
        }
        let effective_references: Vec<String> = references.into_iter()
            .filter(|name| observed.contains_key(name)).collect();
        let effective_comparisons: Vec<String> = comparisons.into_iter()
            .filter(|name| observed.contains_key(name)).collect();
        if effective_references.is_empty() || effective_comparisons.is_empty() {
            if interactive {
                eprintln!("[traconmech] No usable {} group; please select both groups.",
                    if effective_references.is_empty() { "reference" } else { "comparison" });
                correction = true;
                continue;
            }
            eprintln!("[traconmech] WARNING: one or both effective groups are empty; continuing in batch mode");
        }
        // Preserve the full supplied design in YAML; analyze only observed labels.
        break (effective_references, effective_comparisons, source);
    };
    // Keep the scientist's design but recalculate reference matching on each run.
    let experiment = Experiment {
        version: 1,
        population_column: population.into(),
        condition_column: condition.into(),
        observed_conditions: observed,
        reference_conditions: references.clone(),
        comparison_conditions: comparisons.clone(),
        decision_source: source,
        min_reference_cells,
        population_references: Vec::new(),
    };
    if let Some(previous) = &old {
        anyhow::ensure!(previous.reference_conditions == experiment.reference_conditions
            && previous.comparison_conditions == experiment.comparison_conditions,
            "experiment design differs from saved {}; choose a new --out directory or review the saved design", file.display());
    } else {
        fs::create_dir_all(out)?;
        fs::write(&file, serde_yaml::to_string(&experiment)?)?;
    }
    Ok(Design {
        cases: comparisons,
        controls: references,
    })
}
pub fn record_references(out: &Path) -> Result<()> {
    let path = out.join("experiment.yaml");
    let mut experiment: Experiment = serde_yaml::from_str(&fs::read_to_string(&path)?)?;
    let tsv = fs::read_to_string(out.join("wt_reference_clusters.tsv"))?;
    experiment.population_references.clear();
    for line in tsv.lines().skip(1) {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 7 {
            bail!("invalid reference audit row: {line}");
        }
        experiment.population_references.push(Reference {
            population: fields[0].into(),
            case_cells: fields[1].parse()?,
            own_reference_cells: fields[2].parse()?,
            selected_population: (fields[3] != "NA").then(|| fields[3].to_owned()),
            selected_reference_cells: fields[4].parse()?,
            method: fields[5].into(),
            pearson_log1p_mean: if fields[6] == "NA" {
                None
            } else {
                Some(fields[6].parse()?)
            },
        });
    }
    fs::write(path, serde_yaml::to_string(&experiment)?)?;
    Ok(())
}

#[cfg(test)]
mod one_sided_tests {
    use super::{labels, parse_confirmation};

    #[test]
    fn confirmation_accepts_case_insensitive_answers_and_rejects_gibberish() {
        for answer in ["y", "Y", "yes", "YES", " Yes \n"] {
            assert_eq!(parse_confirmation(answer, false), Some(true));
        }
        for answer in ["n", "N", "no", "NO"] {
            assert_eq!(parse_confirmation(answer, true), Some(false));
        }
        assert_eq!(parse_confirmation("", true), Some(true));
        assert_eq!(parse_confirmation("", false), Some(false));
        assert_eq!(parse_confirmation("You do not need to even ask!", true), None);
    }

    #[test]
    fn empty_group_is_explicitly_representable() {
        assert!(labels("").is_empty());
        assert!(labels(" , , ").is_empty());
        assert_eq!(labels("DBA1_3, DBA1_2, DBA1_3"), vec!["DBA1_2", "DBA1_3"]);
    }
}
