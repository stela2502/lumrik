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
    io::stdin().read_line(&mut answer)?;
    Ok(labels(&answer))
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
    fs::create_dir_all(out)?;
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
    let (references, comparisons, source) = if cases.is_some() || controls.is_some() {
        anyhow::ensure!(
            cases.is_some() && controls.is_some(),
            "specify both --case and --control, or neither"
        );
        (
            labels(controls.unwrap()),
            labels(cases.unwrap()),
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
        let references = ask("Reference labels (comma-separated, exact)")?;
        let comparisons = ask("Comparison labels (comma-separated, exact)")?;
        (references, comparisons, "scientist".to_owned())
    } else {
        bail!(
            "experimental groups must be selected explicitly: run interactively or specify --case and --control (observed: {:?})",
            observed
        );
    };
    let ref_set: BTreeSet<_> = references.iter().collect();
    let case_set: BTreeSet<_> = comparisons.iter().collect();
    anyhow::ensure!(
        !ref_set.is_empty() && !case_set.is_empty() && ref_set.is_disjoint(&case_set),
        "reference and comparison labels must be nonempty and disjoint"
    );
    for name in ref_set.union(&case_set) {
        anyhow::ensure!(
            observed.contains_key(*name),
            "unknown condition label {name:?}"
        );
    }
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
    fs::write(&file, serde_yaml::to_string(&experiment)?)?;
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
