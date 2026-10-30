//! Marker annotation package and general-purpose, conservative H5AD obs importer.
use anyhow::{Context, Result, bail, ensure};
use clap::Args;
use hdf5::{File, types::VarLenUnicode};
use std::{collections::{BTreeMap, BTreeSet}, fs, path::{Path, PathBuf}, process::Command, str::FromStr};

const SUMMARIES: [&str; 5] = [
    "marker_expression_stats.tsv", "marker_summary.tsv", "marker_summary_general.tsv",
    "marker_summary_sparse_strong.tsv", "marker_summary_low_expression.tsv",
];
const TODO: &str = r#"# Cell-type annotation task

Read `stats/README.md` first: it documents the actual cluster-versus-rest
statistics, scoring formulas, rank directions, normalization caveats, and
limitations of each of the four marker summary tables.

`stats/marker_expression_stats.tsv` reports the union of the selected
high-variance genes and the genes in the four marker-summary lists.
Each row is one cluster/gene pair. `n_cells` counts all cells in the cluster;
`n_expressing` counts cells with expression > 0; `pct_expressing` is 100 times
n_expressing/n_cells. `mean_expressing` and q10/q25/median/q75/q90 are computed
ONLY among cells with expression > 0, using linearly interpolated empirical
quantiles. A gene with no positive cells has NA conditional statistics.
The expression values are log1p(10,000 * surviving_gene_UMI /
surviving_cell_UMI), as described by the GEX pipeline.
`gene_source` is `variable`, `marker_summary`, or `both`.
Do not interpret a conditional quantile as a whole-cluster quantile.

Infer likely cell identities using combinations of lineage markers, detection
prevalence, conditional expression intensity, and contrasting cluster markers.
Do not infer a definitive cell type from one marker. Distinguish immature
stages only where evidence supports them. Use `Unknown` for ambiguity.

## Required deliverable: cell_types.tsv

Return a UTF-8 tab-separated file, one row per cluster, with columns:

cluster\tcell_type\tconfidence\tsupporting_markers\tcontradicting_markers\treasoning

Copy cluster identifiers EXACTLY from the input (e.g. cluster_000).
Confidence: high, medium, or low. Marker lists: comma-separated.
No Markdown code fences. The scientist reviews the file before importing:

sc-analyse annotate --data analysis.h5ad --annotations cell_types.tsv

By default the importer maps TSV `cluster` to H5AD `obs/cluster` and imports
all other TSV columns into `obs`. Use --columns to select a subset.
"#;

#[derive(Args, Debug)]
pub struct AnnotateCli {
    /// Existing H5AD file to update in place.
    #[arg(long)]
    data: PathBuf,
    /// Tab-separated annotation table (one row per matching key).
    #[arg(long)]
    annotations: PathBuf,
    /// Existing obs column used to match rows (or _index for cell IDs).
    #[arg(long, default_value = "cluster")]
    match_by: String,
    /// Column in TSV containing matching keys.
    #[arg(long, default_value = "cluster")]
    key_column: String,
    /// Comma-separated TSV columns to import; default: all except key.
    #[arg(long, value_delimiter = ',')]
    columns: Option<Vec<String>>,
    /// Permit replacing existing obs columns.
    #[arg(long)]
    replace: bool,
}

pub fn run_annotate(cli: AnnotateCli) -> Result<()> {
    ensure!(cli.data.is_file(), "missing H5AD: {}", cli.data.display());
    import(&cli.data, &cli.annotations, &cli.match_by, &cli.key_column,
           cli.columns.as_deref(), cli.replace)
}

pub fn package(stats: &Path) -> Result<()> {
    for name in SUMMARIES {
        ensure!(stats.join(name).is_file(), "missing marker summary: {}", stats.join(name).display());
    }
    let root = stats.parent().context("stats directory has no parent")?;
    let todo = root.join("TODO.md");
    fs::write(&todo, TODO)?;
    let output = root.join("cell_type_annotation.tar.gz");
    let tmp = root.join("cell_type_annotation.tar.gz.tmp");
    let result = Command::new("tar")
        .current_dir(root).arg("-czf").arg(&tmp).arg("TODO.md").arg("stats/README.md")
        .args(SUMMARIES.map(|x| format!("stats/{x}")))
        .status().context("running tar to package annotation summaries")?;
    if !result.success() {
        let _ = fs::remove_file(&tmp);
        bail!("tar failed to package annotation summaries");
    }
    fs::rename(&tmp, &output)?;
    eprintln!("AI annotation package: {}", output.display());
    Ok(())
}

fn read_string_dataset(ds: &hdf5::Dataset) -> Result<Vec<String>> {
    if let Ok(values) = ds.read_raw::<VarLenUnicode>() {
        return Ok(values.iter().map(ToString::to_string).collect());
    }
    if let Ok(values) = ds.read_raw::<hdf5::types::VarLenAscii>() {
        return Ok(values.iter().map(ToString::to_string).collect());
    }
    macro_rules! numeric {
        ($($ty:ty),*) => { $(if let Ok(values) = ds.read_raw::<$ty>() {
            return Ok(values.iter().map(ToString::to_string).collect());
        })* };
    }
    numeric!(i8, i16, i32, i64, u8, u16, u32, u64);
    bail!("unsupported key dataset type: {}", ds.name())
}

fn read_obs_keys(file: &File, name: &str) -> Result<Vec<String>> {
    let obs = file.group("obs")?;
    if name == "_index" {
        let index_name = obs.attr("_index").ok()
            .and_then(|a| a.read_scalar::<VarLenUnicode>().ok())
            .map(|v| v.to_string()).unwrap_or_else(|| "_index".to_owned());
        return read_string_dataset(&obs.dataset(&index_name)?);
    }
    if let Ok(ds) = obs.dataset(name) { return read_string_dataset(&ds); }
    let group = obs.group(name).with_context(|| format!("missing obs column {name:?}"))?;
    let categories = read_string_dataset(&group.dataset("categories")?)?;
    let codes = group.dataset("codes")?.read_raw::<i32>()?;
    codes.into_iter().map(|code| {
        ensure!(code >= 0, "null categorical key in obs/{name}");
        categories.get(code as usize).cloned()
            .with_context(|| format!("invalid categorical code {code} in obs/{name}"))
    }).collect()
}

fn import(data: &Path, table: &Path, match_by: &str, key_col: &str,
          requested: Option<&[String]>, replace: bool) -> Result<()> {
    eprintln!("[1/4] Loading annotation table: {}", table.display());
    ensure!(table.is_file(), "missing annotation table: {}", table.display());
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(table)
        .with_context(|| format!("opening annotation table {}", table.display()))?;
    let headers = reader.headers()?.clone();
    let key_idx = headers.iter().position(|h| h == key_col)
        .with_context(|| format!("missing TSV key column {key_col:?}"))?;
    let selected: Vec<(usize, String)> = match requested {
        Some(names) => names.iter().map(|name| {
            let idx = headers.iter().position(|h| h == name)
                .with_context(|| format!("missing TSV column {name:?}"))?;
            ensure!(idx != key_idx, "cannot import key column {key_col:?}");
            Ok((idx, name.clone()))
        }).collect::<Result<_>>()?,
        None => headers.iter().enumerate().filter(|(i, _)| *i != key_idx)
            .map(|(i, name)| (i, name.to_owned())).collect(),
    };
    ensure!(!selected.is_empty(), "no annotation columns selected");
    let mut unique = BTreeSet::new();
    for (_, name) in &selected {
        ensure!(name != "_index" && name != match_by, "cannot overwrite matching/index column {name:?}");
        ensure!(unique.insert(name.clone()), "duplicate selected column {name:?}");
    }
    let mut rows = BTreeMap::<String, Vec<String>>::new();
    for record in reader.records() {
        let record = record?;
        let key = record.get(key_idx).unwrap_or("").trim();
        ensure!(!key.is_empty(), "empty annotation key");
        let values = selected.iter().map(|(idx, _)| record.get(*idx).unwrap_or("").to_owned()).collect();
        ensure!(rows.insert(key.to_owned(), values).is_none(), "duplicate annotation key {key:?}");
    }
    ensure!(!rows.is_empty(), "empty annotation table");
    eprintln!("      {} annotation keys; columns: {}", rows.len(),
        selected.iter().map(|(_, name)| name.as_str()).collect::<Vec<_>>().join(", "));
    eprintln!("[2/4] Loading H5AD and validating annotation keys: {}", data.display());
    let source = File::open(data)?;
    let keys = read_obs_keys(&source, match_by)?;
    let observed: BTreeSet<_> = keys.iter().cloned().collect();
    let proposed: BTreeSet<_> = rows.keys().cloned().collect();
    ensure!(observed == proposed, "key mismatch: missing {:?}; unexpected {:?}",
        observed.difference(&proposed).collect::<Vec<_>>(),
        proposed.difference(&observed).collect::<Vec<_>>());
    eprintln!("      {} cells; {} distinct keys; all annotation keys matched",
        keys.len(), observed.len());
    let obs = source.group("obs")?;
    for (_, name) in &selected {
        ensure!(replace || !obs.link_exists(name), "obs/{name} exists; use --replace");
    }
    drop(obs);
    drop(source);
    eprintln!("[3/4] Applying annotations directly to H5AD (no full-file copy)");
    eprintln!("      Writing {} annotation columns...", selected.len());
    let write_result = (|| -> Result<()> {
        let file = File::open_rw(data)?;
        let obs = file.group("obs")?;
        let mut order: Vec<String> = obs.attr("column-order")?.read_raw::<VarLenUnicode>()?
            .iter().map(ToString::to_string).collect();
        for (column_idx, (_, name)) in selected.iter().enumerate() {
            if obs.link_exists(name) { obs.unlink(name)?; }
            let values: Vec<VarLenUnicode> = keys.iter()
                .map(|key| VarLenUnicode::from_str(&rows[key][column_idx]))
                .collect::<std::result::Result<_, _>>()?;
            let ds = obs.new_dataset::<VarLenUnicode>().shape(values.len()).create(name.as_str())?;
            ds.write(&values)?;
            ds.new_attr::<VarLenUnicode>().shape(()).create("encoding-type")?
                .write_scalar(&VarLenUnicode::from_str("string-array")?)?;
            ds.new_attr::<VarLenUnicode>().shape(()).create("encoding-version")?
                .write_scalar(&VarLenUnicode::from_str("0.2.0")?)?;
            if !order.contains(name) { order.push(name.clone()); }
        }
        let order_vals: Vec<VarLenUnicode> = order.iter().map(|s| VarLenUnicode::from_str(s))
            .collect::<std::result::Result<_, _>>()?;
        let cname = std::ffi::CString::new("column-order")?;
        let rc = unsafe { hdf5_sys::h5a::H5Adelete(obs.id(), cname.as_ptr()) };
        ensure!(rc >= 0, "failed to replace /obs column-order");
        obs.new_attr::<VarLenUnicode>().shape(order_vals.len())
            .create("column-order")?.write(&order_vals)?;
        file.flush()?;
        Ok(())
    })();
    write_result.context("writing annotation columns in place (file may be partially modified on error)")?;
    eprintln!("      Saved in place: {}", data.display());
    eprintln!("[4/4] Annotation summary (cells per imported label)");
    for (column_idx, (_, name)) in selected.iter().enumerate() {
        let mut counts = BTreeMap::<&str, usize>::new();
        for key in &keys {
            *counts.entry(rows[key][column_idx].as_str()).or_default() += 1;
        }
        eprintln!("      {}:", name);
        for (label, count) in counts {
            eprintln!("        {:<40} {:>8}", if label.is_empty() { "<empty>" } else { label }, count);
        }
        eprintln!("        {:<40} {:>8}", "TOTAL", keys.len());
    }
    Ok(())
}
