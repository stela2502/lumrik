//! Conservative H5AD row splitter: never silently discard unsupported HDF5 content.
use anyhow::{Context, Result, bail, ensure};
use clap::Args;
use sc_anndata::{AnnData, Column, Matrix, Table};
use std::{collections::BTreeMap, fs, path::PathBuf, time::Instant};

#[derive(Args, Debug)]
pub struct SplitCli {
    /// Source H5AD; never modified.
    #[arg(long)]
    data: PathBuf,
    /// obs column whose distinct values define the outputs.
    #[arg(long)]
    by: String,
    /// Destination directory; existing files are never overwritten.
    #[arg(long, default_value = "split")]
    out: PathBuf,
}

fn subset_column(c: &Column, rows: &[usize]) -> Column {
    match c {
        Column::String(v) => Column::String(rows.iter().map(|&i| v[i].clone()).collect()),
        Column::F64(v) => Column::F64(rows.iter().map(|&i| v[i]).collect()),
        Column::I64(v) => Column::I64(rows.iter().map(|&i| v[i]).collect()),
        Column::Bool(v) => Column::Bool(rows.iter().map(|&i| v[i]).collect()),
        Column::NullableString(v) => Column::NullableString(rows.iter().map(|&i| v[i].clone()).collect()),
        Column::NullableI64(v) => Column::NullableI64(rows.iter().map(|&i| v[i]).collect()),
        Column::NullableBool(v) => Column::NullableBool(rows.iter().map(|&i| v[i]).collect()),
        Column::Categorical { categories, codes, ordered } => Column::Categorical {
            categories: categories.clone(), codes: rows.iter().map(|&i| codes[i]).collect(), ordered: *ordered,
        },
    }
}

fn subset_table(t: &Table, rows: &[usize]) -> Result<Table> {
    let mut result = Table::new(t.index_name.clone(), rows.iter().map(|&i| t.index[i].clone()).collect());
    for (name, column) in &t.columns {
        result.add_column(name.clone(), subset_column(column, rows))?;
    }
    Ok(result)
}

fn subset_matrix(m: &Matrix, rows: &[usize]) -> Matrix {
    match m {
        Matrix::DenseF32(a) => Matrix::DenseF32(a.select(ndarray::Axis(0), rows)),
        Matrix::CsrF32(a) => {
            let csr = a.to_csr();
            let mut indptr = Vec::with_capacity(rows.len() + 1);
            let mut indices = Vec::new();
            let mut data = Vec::new();
            indptr.push(0);
            for &r in rows {
                if let Some(row) = csr.outer_view(r) {
                    for (col, &value) in row.iter() {
                        indices.push(col);
                        data.push(value);
                    }
                }
                indptr.push(indices.len());
            }
            Matrix::CsrF32(sprs::CsMat::new((rows.len(), a.cols()), indptr, indices, data))
        }
    }
}

fn safe_name(value: &str) -> String {
    let name: String = value.chars().map(|c| if c.is_ascii_alphanumeric() || c == '-' || c == '_' { c } else { '_' }).collect();
    let name = name.trim_matches('_');
    if name.is_empty() { "unnamed".to_owned() } else { name.chars().take(90).collect() }
}

pub fn run_split(cli: SplitCli) -> Result<()> {
    let started = Instant::now();
    eprintln!("[1/4] Opening H5AD: {}", cli.data.display());
    ensure!(cli.data.is_file(), "Input H5AD not found: {}", cli.data.display());
    let source = hdf5::File::open(&cli.data)?;
    // sc_anndata currently drops these groups on write. Refuse, with actionable error,
    // instead of silently destroying data or producing misleading subsets.
    for name in ["uns", "obsp", "varm", "varp", "raw"] {
        if source.link_exists(name) {
            let group = source.group(name);
            if name == "raw" || group.as_ref().map_or(true, |g| g.member_names().map_or(true, |v| !v.is_empty())) {
                bail!("Cannot safely split {}: /{} contains data that the current H5AD writer cannot preserve. Source remains untouched; no outputs created.", cli.data.display(), name);
            }
        }
    }
    eprintln!("      Reading expression matrix and metadata...");
    let read = AnnData::read_h5ad(&cli.data)?;
    // sc_anndata currently warns on empty structural groups. Only suppress
    // those exact warnings when the corresponding group is verified empty.
    let unsupported: Vec<_> = read.warnings.iter().filter(|warning| {
        let name = warning.path.strip_prefix('/').unwrap_or("");
        if warning.reason != "not implemented yet" || !["uns", "obsp", "varm", "varp"].contains(&name) {
            return true;
        }
        match source.group(name) {
            Ok(group) => group.member_names().map_or(true, |members| !members.is_empty()),
            Err(_) => true,
        }
    }).collect();
    ensure!(unsupported.is_empty(), "Cannot safely split H5AD: unsupported content: {:?}. No outputs created.", unsupported);
    let adata = read.data;
    ensure!(adata.obs.declared_columns.len() == adata.obs.columns.len(), "Cannot safely split: some obs columns were not decoded");
    ensure!(adata.var.declared_columns.len() == adata.var.columns.len(), "Cannot safely split: some var columns were not decoded");
    eprintln!("      Loaded {} cells; no unsupported content", adata.obs.len());
    eprintln!("[2/4] Grouping cells by obs/{}", cli.by);
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for row in 0..adata.obs.len() {
        let label = if cli.by == "_index" { Some(adata.obs.index[row].clone()) }
            else { adata.obs.columns.get(&cli.by).and_then(|c| c.get_string(row)) };
        let label = label.with_context(|| format!("Missing split value in obs/{} at row {}", cli.by, row))?;
        groups.entry(label).or_default().push(row);
    }
    ensure!(!groups.is_empty(), "No cells to split");
    eprintln!("      Found {} groups covering {} cells", groups.len(), adata.obs.len());
    let mut destinations = Vec::new();
    for (label, rows) in &groups {
        let dest = cli.out.join(format!("{}.h5ad", safe_name(label)));
        ensure!(!dest.exists(), "Output already exists: {}. Nothing overwritten; choose a different --out.", dest.display());
        ensure!(!destinations.iter().any(|(other, _, _): &(PathBuf, String, Vec<usize>)| other == &dest), "Distinct values produce the same filename: {}. Rename labels first.", dest.display());
        destinations.push((dest, label.clone(), rows.clone()));
    }
    fs::create_dir_all(&cli.out).with_context(|| format!("creating {}", cli.out.display()))?;
    eprintln!("[3/4] Writing {} split H5AD files", destinations.len());
    let total_files = destinations.len();
    let mut written_cells = 0usize;
    for (dest, label, rows) in destinations {
        let result = AnnData {
            x: subset_matrix(&adata.x, &rows),
            obs: subset_table(&adata.obs, &rows)?,
            var: adata.var.clone(),
            obsm: adata.obsm.iter().map(|(name, m)| (name.clone(), subset_matrix(m, &rows))).collect(),
            layers: adata.layers.iter().map(|(name, m)| (name.clone(), subset_matrix(m, &rows))).collect(),
        };
        // Create-new semantics prevent clobbering even if another process raced us.
        let reservation = fs::OpenOptions::new().write(true).create_new(true).open(&dest)
            .with_context(|| format!("refusing to overwrite {}", dest.display()))?;
        drop(reservation);
        if let Err(error) = result.write_h5ad(&dest) {
            let _ = fs::remove_file(&dest);
            return Err(error).with_context(|| format!("writing subset {label:?}"));
        }
        written_cells += rows.len();
        eprintln!("      {:>6} cells | {} | written to {}", rows.len(), label, dest.display());
    }
    eprintln!("[4/4] Split summary");
    eprintln!("      Output files: {}", total_files);
    eprintln!("      Cells written: {} / {}", written_cells, adata.obs.len());
    eprintln!("      Elapsed: {:.2}s", started.elapsed().as_secs_f64());
    Ok(())
}
