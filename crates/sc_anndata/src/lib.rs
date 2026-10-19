//! AnnData/H5AD interchange for Lumrik.
//!
//! The contract is deliberately asymmetric:
//! - writing is strict: dimensions and supported values must be valid;
//! - reading is permissive: unsupported H5AD members are reported and skipped.

use anyhow::{Context, Result, bail};
use data_table::DataTable;
use hdf5::{
    File, Group,
    types::{VarLenAscii, VarLenUnicode},
};
use ndarray::Array2;
use sprs::CsMat;
use std::{collections::BTreeMap, path::Path, str::FromStr};

#[derive(Debug, Clone)]
pub enum Matrix {
    CsrF32(CsMat<f32>),
    DenseF32(Array2<f32>),
}

impl Matrix {
    pub fn shape(&self) -> (usize, usize) {
        match self {
            Self::CsrF32(x) => x.shape(),
            Self::DenseF32(x) => x.dim(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum Column {
    String(Vec<String>),
    F64(Vec<f64>),
    I64(Vec<i64>),
    Bool(Vec<bool>),
    /// Native AnnData categorical: codes are zero-based; -1 means missing.
    Categorical {
        categories: Vec<String>,
        codes: Vec<i64>,
        ordered: bool,
    },
}

impl Column {
    fn len(&self) -> usize {
        match self {
            Self::String(x) => x.len(),
            Self::F64(x) => x.len(),
            Self::I64(x) => x.len(),
            Self::Bool(x) => x.len(),
            Self::Categorical { codes, .. } => codes.len(),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Table {
    pub index_name: String,
    pub index: Vec<String>,
    pub columns: BTreeMap<String, Column>,
    /// Declared H5AD dataframe columns, including columns with unsupported encodings.
    pub declared_columns: BTreeMap<String, String>,
}

impl Table {
    pub fn new(index_name: impl Into<String>, index: Vec<String>) -> Self {
        Self {
            index_name: index_name.into(),
            index,
            columns: BTreeMap::new(),
            declared_columns: BTreeMap::new(),
        }
    }
    pub fn add_column(&mut self, name: impl Into<String>, column: Column) -> Result<()> {
        let name = name.into();
        if column.len() != self.index.len() {
            bail!(
                "column {name:?} has {} rows, expected {}",
                column.len(),
                self.index.len()
            );
        }
        self.declared_columns.insert(name.clone(), "decoded".into());
        self.columns.insert(name, column);
        Ok(())
    }
    /// Decode categorical values into the existing data_table::Factor representation.
    /// AnnData -1 becomes Factor's NaN; category ordering is preserved.
    pub fn factor(&self, name: &str) -> Result<Option<(data_table::Factor, Vec<f64>, bool)>> {
        let Some(Column::Categorical {
            categories,
            codes,
            ordered,
        }) = self.columns.get(name)
        else {
            return Ok(None);
        };
        let mut factor = data_table::Factor::new(name, false);
        for category in categories {
            if category.is_empty() || category.eq_ignore_ascii_case("NA") {
                bail!(
                    "factor {name:?}: category {category:?} conflicts with data_table missing-value convention"
                );
            }
            factor.level_to_index(category);
        }
        let values = codes
            .iter()
            .map(|&code| {
                if code == -1 {
                    Ok(f64::NAN)
                } else if code >= 0 && (code as usize) < categories.len() {
                    Ok(code as f64)
                } else {
                    bail!("factor {name:?}: invalid categorical code {code}")
                }
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Some((factor, values, *ordered)))
    }
    pub fn len(&self) -> usize {
        self.index.len()
    }
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }
    pub fn from_data_table(table: &DataTable, index_column: &str) -> Result<Self> {
        let index_col = table
            .headers
            .iter()
            .position(|x| x == index_column)
            .with_context(|| format!("DataTable has no index column {index_column:?}"))?;
        let index = if let Some(factor) = table.factors.get(index_column) {
            table
                .numeric_data
                .column(index_col)
                .iter()
                .map(|&x| factor.get_string(x))
                .collect()
        } else {
            table
                .numeric_data
                .column(index_col)
                .iter()
                .map(|x| x.to_string())
                .collect()
        };
        let mut out = Self::new(index_column, index);
        for (col, name) in table.headers.iter().enumerate() {
            if col == index_col || table.exclude.contains(name) {
                continue;
            }
            if let Some(factor) = table.factors.get(name) {
                out.add_column(
                    name.clone(),
                    Column::String(
                        table
                            .numeric_data
                            .column(col)
                            .iter()
                            .map(|&x| factor.get_string(x))
                            .collect(),
                    ),
                )?;
            } else {
                out.add_column(
                    name.clone(),
                    Column::F64(table.numeric_data.column(col).to_vec()),
                )?;
            }
        }
        Ok(out)
    }
}

#[derive(Debug, Clone)]
pub struct AnnData {
    pub x: Matrix,
    pub obs: Table,
    pub var: Table,
    pub obsm: BTreeMap<String, Matrix>,
    pub layers: BTreeMap<String, Matrix>,
}

impl AnnData {
    pub fn validate(&self) -> Result<()> {
        let (n_obs, n_var) = self.x.shape();
        if self.obs.len() != n_obs {
            bail!("X has {n_obs} observations but obs has {}", self.obs.len());
        }
        if self.var.len() != n_var {
            bail!("X has {n_var} variables but var has {}", self.var.len());
        }
        for (name, x) in &self.obsm {
            if x.shape().0 != n_obs {
                bail!("obsm/{name} has {} rows but obs has {n_obs}", x.shape().0);
            }
        }
        for (name, x) in &self.layers {
            if x.shape() != (n_obs, n_var) {
                bail!(
                    "layers/{name} has shape {:?}, expected ({n_obs}, {n_var})",
                    x.shape()
                );
            }
        }
        Ok(())
    }

    pub fn write_h5ad(&self, path: impl AsRef<Path>) -> Result<()> {
        self.validate()?;
        let file = File::create(path.as_ref())
            .with_context(|| format!("creating {}", path.as_ref().display()))?;
        attr_str(&file, "encoding-type", "anndata")?;
        attr_str(&file, "encoding-version", "0.1.0")?;
        write_matrix(&file, "X", &self.x)?;
        write_table(&file, "obs", &self.obs)?;
        write_table(&file, "var", &self.var)?;
        let obsm = file.create_group("obsm")?;
        attr_str(&obsm, "encoding-type", "dict")?;
        attr_str(&obsm, "encoding-version", "0.1.0")?;
        for (name, x) in &self.obsm {
            write_matrix(&obsm, name, x)?;
        }
        let layers = file.create_group("layers")?;
        attr_str(&layers, "encoding-type", "dict")?;
        attr_str(&layers, "encoding-version", "0.1.0")?;
        for (name, x) in &self.layers {
            write_matrix(&layers, name, x)?;
        }
        for name in ["uns", "obsp", "varm", "varp"] {
            let g = file.create_group(name)?;
            attr_str(&g, "encoding-type", "dict")?;
            attr_str(&g, "encoding-version", "0.1.0")?;
        }
        Ok(())
    }

    pub fn read_h5ad(path: impl AsRef<Path>) -> Result<ReadResult> {
        let file = File::open(path.as_ref())
            .with_context(|| format!("opening {}", path.as_ref().display()))?;
        let mut warnings = Vec::new();
        let x = read_matrix(&file, "X", &mut warnings)?.context("H5AD has no supported /X")?;
        let obs = read_table(&file, "obs", &mut warnings)?.context("H5AD has no supported /obs")?;
        let var = read_table(&file, "var", &mut warnings)?.context("H5AD has no supported /var")?;
        let mut obsm = BTreeMap::new();
        read_matrix_map(&file, "obsm", &mut obsm, &mut warnings)?;
        let mut layers = BTreeMap::new();
        read_matrix_map(&file, "layers", &mut layers, &mut warnings)?;
        for name in ["uns", "obsp", "varm", "varp"] {
            if file.link_exists(name) {
                warnings.push(ReadWarning::unsupported(
                    format!("/{name}"),
                    "not implemented yet",
                ));
            }
        }
        let data = Self {
            x,
            obs,
            var,
            obsm,
            layers,
        };
        data.validate()?;
        Ok(ReadResult { data, warnings })
    }
}

/// A format-neutral overview of the expression data held in AnnData.
/// The modality-specific interpretation belongs to downstream analysis crates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AnnotationSummary {
    pub name: String,
    pub data_type: &'static str,
    pub distinct_values: Option<usize>,
    pub counts: Option<Vec<(String, usize)>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetOverview {
    pub cells: usize,
    pub features: usize,
    pub obs: Vec<AnnotationSummary>,
    pub var: Vec<AnnotationSummary>,
}

impl AnnData {
    /// Summarize the loaded dataset without requiring a particular annotation.
    /// Up to `max_categories` categorical values are included with their counts.
    pub fn overview(&self, max_categories: usize) -> Result<DatasetOverview> {
        self.validate()?;
        let (cells, features) = self.x.shape();
        Ok(DatasetOverview {
            cells,
            features,
            obs: summarize_columns(&self.obs, max_categories),
            var: summarize_columns(&self.var, max_categories),
        })
    }
}

fn summarize_columns(table: &Table, max_categories: usize) -> Vec<AnnotationSummary> {
    table
        .columns
        .iter()
        .map(|(name, col)| {
            let (data_type, counts): (&'static str, Option<BTreeMap<String, usize>>) = match col {
                Column::String(values) => {
                    let mut counts = BTreeMap::new();
                    for value in values {
                        *counts.entry(value.clone()).or_insert(0) += 1;
                    }
                    ("string", Some(counts))
                }
                Column::Bool(values) => {
                    let mut counts = BTreeMap::new();
                    for value in values {
                        *counts.entry(value.to_string()).or_insert(0) += 1;
                    }
                    ("bool", Some(counts))
                }
                Column::I64(values) => {
                    let mut counts = BTreeMap::new();
                    for value in values {
                        *counts.entry(value.to_string()).or_insert(0) += 1;
                    }
                    ("integer", Some(counts))
                }
                Column::Categorical {
                    categories, codes, ..
                } => {
                    let mut counts = BTreeMap::new();
                    for &code in codes {
                        let label = if code == -1 {
                            "NA"
                        } else {
                            categories
                                .get(code as usize)
                                .map(String::as_str)
                                .unwrap_or("NA")
                        };
                        *counts.entry(label.to_string()).or_insert(0) += 1;
                    }
                    ("factor", Some(counts))
                }
                Column::F64(_) => ("float", None),
            };
            let distinct_values = counts.as_ref().map(BTreeMap::len);
            let counts = counts.and_then(|counts| {
                if counts.len() > max_categories {
                    return None;
                }
                let mut pairs: Vec<_> = counts.into_iter().collect();
                pairs.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
                Some(pairs)
            });
            AnnotationSummary {
                name: name.clone(),
                data_type,
                distinct_values,
                counts,
            }
        })
        .chain(
            table
                .declared_columns
                .iter()
                .filter(|(name, _)| !table.columns.contains_key(*name))
                .map(|(name, _reason)| AnnotationSummary {
                    name: name.clone(),
                    data_type: "unsupported",
                    distinct_values: None,
                    counts: None,
                }),
        )
        .collect()
}

#[derive(Debug, Clone)]
pub struct ReadWarning {
    pub path: String,
    pub reason: String,
}
impl ReadWarning {
    fn unsupported(path: impl Into<String>, reason: impl Into<String>) -> Self {
        Self {
            path: path.into(),
            reason: reason.into(),
        }
    }
}
#[derive(Debug, Clone)]
pub struct ReadResult {
    pub data: AnnData,
    pub warnings: Vec<ReadWarning>,
}

fn s(x: &str) -> VarLenUnicode {
    VarLenUnicode::from_str(x).expect("Rust str is valid UTF-8")
}
fn attr_str(loc: &hdf5::Location, name: &str, value: &str) -> Result<()> {
    loc.new_attr::<VarLenUnicode>()
        .shape(())
        .create(name)?
        .as_writer()
        .write_scalar(&s(value))?;
    Ok(())
}
fn attr_strings(loc: &hdf5::Location, name: &str, values: &[String]) -> Result<()> {
    let vals: Vec<_> = values.iter().map(|x| s(x)).collect();
    loc.new_attr::<VarLenUnicode>()
        .shape(vals.len())
        .create(name)?
        .write(&vals)?;
    Ok(())
}
fn write_strings(group: &Group, name: &str, values: &[String]) -> Result<()> {
    let vals: Vec<_> = values.iter().map(|x| s(x)).collect();
    let ds = group
        .new_dataset::<VarLenUnicode>()
        .shape(vals.len())
        .create(name)?;
    ds.write(&vals)?;
    attr_str(&ds, "encoding-type", "string-array")?;
    attr_str(&ds, "encoding-version", "0.2.0")?;
    Ok(())
}
fn write_table(file: &File, name: &str, table: &Table) -> Result<()> {
    let g = file.create_group(name)?;
    attr_str(&g, "encoding-type", "dataframe")?;
    attr_str(&g, "encoding-version", "0.2.0")?;
    attr_str(&g, "_index", &table.index_name)?;
    let names: Vec<String> = table.columns.keys().cloned().collect();
    attr_strings(&g, "column-order", &names)?;
    write_strings(&g, &table.index_name, &table.index)?;
    for (name, col) in &table.columns {
        match col {
            Column::String(x) => write_strings(&g, name, x)?,
            Column::Categorical {
                categories,
                codes,
                ordered,
            } => {
                if codes
                    .iter()
                    .any(|&code| code < -1 || (code >= 0 && code as usize >= categories.len()))
                {
                    bail!("/{}/{}: categorical code out of range", g.name(), name);
                }
                let cat = g.create_group(name)?;
                attr_str(&cat, "encoding-type", "categorical")?;
                attr_str(&cat, "encoding-version", "0.2.0")?;
                cat.new_attr::<u8>()
                    .create("ordered")?
                    .write_scalar(&u8::from(*ordered))?;
                write_strings(&cat, "categories", categories)?;
                cat.new_dataset_builder().with_data(codes).create("codes")?;
            }
            Column::F64(x) => {
                let ds = g
                    .new_dataset::<f64>()
                    .shape(x.len())
                    .create(name.as_str())?;
                ds.write(x)?;
                attr_str(&ds, "encoding-type", "array")?;
                attr_str(&ds, "encoding-version", "0.2.0")?;
            }
            Column::I64(x) => {
                let ds = g
                    .new_dataset::<i64>()
                    .shape(x.len())
                    .create(name.as_str())?;
                ds.write(x)?;
                attr_str(&ds, "encoding-type", "array")?;
                attr_str(&ds, "encoding-version", "0.2.0")?;
            }
            Column::Bool(x) => {
                let vals: Vec<u8> = x.iter().map(|v| u8::from(*v)).collect();
                let ds = g
                    .new_dataset::<u8>()
                    .shape(vals.len())
                    .create(name.as_str())?;
                ds.write(&vals)?;
                attr_str(&ds, "encoding-type", "array")?;
                attr_str(&ds, "encoding-version", "0.2.0")?;
            }
        }
    }
    Ok(())
}
fn write_matrix(parent: &Group, name: &str, matrix: &Matrix) -> Result<()> {
    match matrix {
        Matrix::DenseF32(x) => {
            let ds = parent.new_dataset::<f32>().shape(x.dim()).create(name)?;
            ds.write(x.view())?;
            attr_str(&ds, "encoding-type", "array")?;
            attr_str(&ds, "encoding-version", "0.2.0")?;
        }
        Matrix::CsrF32(x) => {
            let x = if x.is_csr() { x.clone() } else { x.to_csr() };
            let g = parent.create_group(name)?;
            attr_str(&g, "encoding-type", "csr_matrix")?;
            attr_str(&g, "encoding-version", "0.1.0")?;
            g.new_attr::<i64>()
                .shape(2)
                .create("shape")?
                .write(&[x.rows() as i64, x.cols() as i64])?;
            g.new_dataset::<f32>()
                .shape(x.data().len())
                .create("data")?
                .write(x.data())?;
            let indices: Vec<i32> = x.indices().iter().map(|&v| v as i32).collect();
            g.new_dataset::<i32>()
                .shape(indices.len())
                .create("indices")?
                .write(&indices)?;
            let indptr: Vec<i64> = x.indptr().raw_storage().iter().map(|&v| v as i64).collect();
            g.new_dataset::<i64>()
                .shape(indptr.len())
                .create("indptr")?
                .write(&indptr)?;
        }
    }
    Ok(())
}
fn read_matrix(
    parent: &Group,
    name: &str,
    warnings: &mut Vec<ReadWarning>,
) -> Result<Option<Matrix>> {
    if !parent.link_exists(name) {
        return Ok(None);
    }
    if let Ok(ds) = parent.dataset(name) {
        if ds.ndim() == 2 && ds.dtype()?.is::<f32>() {
            return Ok(Some(Matrix::DenseF32(ds.read_2d::<f32>()?)));
        }
        warnings.push(ReadWarning::unsupported(
            format!("/{name}"),
            format!(
                "dense dtype/dim {:?}/{} unsupported",
                ds.dtype()?,
                ds.ndim()
            ),
        ));
        return Ok(None);
    }
    let g = parent.group(name)?;
    let enc = read_attr_string(&g, "encoding-type").unwrap_or_default();
    if enc != "csr_matrix" && enc != "csc_matrix" {
        warnings.push(ReadWarning::unsupported(
            format!("/{name}"),
            format!("matrix encoding {enc:?} unsupported"),
        ));
        return Ok(None);
    }
    let shape: Vec<i64> = g.attr("shape")?.read_raw()?;
    if shape.len() != 2 {
        bail!("/{name}: malformed sparse shape");
    }
    let data: Vec<f32> = g.dataset("data")?.read_raw()?;
    let idx: Vec<i32> = g.dataset("indices")?.read_raw()?;
    let ptr: Vec<i64> = g.dataset("indptr")?.read_raw()?;
    // H5AD permits sparse indices in arbitrary order. Validate signed on-disk
    // coordinates before converting, then let sprs sort index/value pairs.
    let rows = usize::try_from(shape[0]).context("negative sparse matrix row count")?;
    let cols = usize::try_from(shape[1]).context("negative sparse matrix column count")?;
    let idx: Vec<usize> = idx
        .into_iter()
        .map(|x| usize::try_from(x).context("negative sparse matrix index"))
        .collect::<Result<_>>()?;
    let ptr: Vec<usize> = ptr
        .into_iter()
        .map(|x| usize::try_from(x).context("negative sparse matrix offset"))
        .collect::<Result<_>>()?;
    let sparse = if enc == "csr_matrix" {
        CsMat::new_from_unsorted((rows, cols), ptr, idx, data)
            .map_err(|(_, _, _, err)| anyhow::anyhow!("/{name}: invalid CSR matrix: {err:?}"))?
    } else {
        CsMat::new_from_unsorted_csc((rows, cols), ptr, idx, data)
            .map_err(|(_, _, _, err)| anyhow::anyhow!("/{name}: invalid CSC matrix: {err:?}"))?
            .to_csr()
    };
    Ok(Some(Matrix::CsrF32(sparse)))
}
fn read_matrix_map(
    file: &File,
    name: &str,
    out: &mut BTreeMap<String, Matrix>,
    warnings: &mut Vec<ReadWarning>,
) -> Result<()> {
    if !file.link_exists(name) {
        return Ok(());
    }
    let g = file.group(name)?;
    for member in g.member_names()? {
        match read_matrix(&g, &member, warnings) {
            Ok(Some(x)) => {
                out.insert(member, x);
            }
            Ok(None) => {}
            Err(e) => warnings.push(ReadWarning::unsupported(
                format!("/{name}/{member}"),
                e.to_string(),
            )),
        }
    }
    Ok(())
}
fn read_attr_string(loc: &hdf5::Location, name: &str) -> Option<String> {
    loc.attr(name)
        .ok()?
        .read_scalar::<VarLenUnicode>()
        .ok()
        .map(|x| x.as_str().to_owned())
}
fn read_attr_strings(loc: &hdf5::Location, name: &str) -> Option<Vec<String>> {
    let attr = loc.attr(name).ok()?;
    if let Ok(values) = attr.read_raw::<VarLenUnicode>() {
        return Some(values.into_iter().map(|v| v.as_str().to_owned()).collect());
    }
    attr.read_raw::<VarLenAscii>()
        .ok()
        .map(|values| values.into_iter().map(|v| v.as_str().to_owned()).collect())
}

fn read_strings(ds: &hdf5::Dataset) -> Result<Option<Vec<String>>> {
    if ds.dtype()?.is::<VarLenUnicode>() {
        return Ok(Some(
            ds.read_raw::<VarLenUnicode>()?
                .into_iter()
                .map(|v| v.as_str().to_owned())
                .collect(),
        ));
    }
    if ds.dtype()?.is::<VarLenAscii>() {
        return Ok(Some(
            ds.read_raw::<VarLenAscii>()?
                .into_iter()
                .map(|v| v.as_str().to_owned())
                .collect(),
        ));
    }
    Ok(None)
}

fn read_codes(ds: &hdf5::Dataset) -> Result<Option<Vec<i64>>> {
    macro_rules! codes {
        ($($ty:ty),*) => {$(if ds.dtype()?.is::<$ty>() {
            return Ok(Some(ds.read_raw::<$ty>()?.into_iter().map(i64::from).collect()));
        })*};
    }
    codes!(i8, i16, i32, i64, u8, u16, u32);
    Ok(None)
}

fn decode_categorical(
    codes: &hdf5::Dataset,
    categories: &hdf5::Dataset,
    ordered: bool,
) -> Result<Option<Column>> {
    let Some(codes) = read_codes(codes)? else {
        return Ok(None);
    };
    let Some(categories) = read_strings(categories)? else {
        return Ok(None);
    };
    for &code in &codes {
        if code < -1 || (code >= 0 && code as usize >= categories.len()) {
            bail!(
                "categorical code {code} out of range for {} categories",
                categories.len()
            );
        }
    }
    Ok(Some(Column::Categorical {
        categories,
        codes,
        ordered,
    }))
}

fn read_table(file: &File, name: &str, warnings: &mut Vec<ReadWarning>) -> Result<Option<Table>> {
    if !file.link_exists(name) {
        return Ok(None);
    }
    let g = file.group(name)?;
    let index_name = read_attr_string(&g, "_index").unwrap_or_else(|| "_index".into());
    let index = read_string_dataset(&g, &index_name)?
        .with_context(|| format!("/{name}/{index_name}: unsupported index"))?;
    let mut table = Table::new(index_name.clone(), index);
    let mut names = read_attr_strings(&g, "column-order").unwrap_or_default();
    // Older AnnData files may not declare column-order; discover members instead.
    for member in g.member_names()? {
        if member != index_name && member != "__categories" && !names.contains(&member) {
            names.push(member);
        }
    }
    for member in names {
        if member == index_name || member == "__categories" {
            continue;
        }
        let path = format!("/{name}/{member}");
        table
            .declared_columns
            .insert(member.clone(), "unsupported encoding".into());
        if !g.link_exists(&member) {
            warnings.push(ReadWarning::unsupported(path, "declared column is missing"));
            continue;
        }
        let col = if let Ok(ds) = g.dataset(&member) {
            if let Some(values) = read_strings(&ds)? {
                Some(Column::String(values))
            } else if ds.dtype()?.is::<f64>() {
                Some(Column::F64(ds.read_raw()?))
            } else if ds.dtype()?.is::<i64>() {
                Some(Column::I64(ds.read_raw()?))
            } else if ds.dtype()?.is::<f32>() {
                Some(Column::F64(
                    ds.read_raw::<f32>()?
                        .into_iter()
                        .map(|x| x as f64)
                        .collect(),
                ))
            } else if let Some(cats) = g
                .group("__categories")
                .ok()
                .and_then(|c| c.dataset(&member).ok())
            {
                decode_categorical(&ds, &cats, false)?
            } else {
                None
            }
        } else if let Ok(group) = g.group(&member) {
            if group.link_exists("codes") && group.link_exists("categories") {
                decode_categorical(
                    &group.dataset("codes")?,
                    &group.dataset("categories")?,
                    group
                        .attr("ordered")
                        .ok()
                        .and_then(|a| a.read_scalar::<u8>().ok())
                        .unwrap_or(0)
                        != 0,
                )?
            } else {
                None
            }
        } else {
            None
        };
        if let Some(col) = col {
            if col.len() != table.len() {
                bail!(
                    "/{name}/{member}: {} rows, expected {}",
                    col.len(),
                    table.len()
                );
            }
            table.columns.insert(member.clone(), col);
            table.declared_columns.insert(member, "decoded".into());
        } else {
            warnings.push(ReadWarning::unsupported(
                path,
                "column encoding not decoded",
            ));
        }
    }
    Ok(Some(table))
}
fn read_string_dataset(g: &Group, name: &str) -> Result<Option<Vec<String>>> {
    if !g.link_exists(name) {
        return Ok(None);
    }
    read_strings(&g.dataset(name)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reads_unsorted_h5ad_sparse_indices() -> Result<()> {
        for encoding in ["csr_matrix", "csc_matrix"] {
            let path = std::env::temp_dir().join(format!(
                "sc_anndata_unsorted_{}_{}.h5ad",
                encoding,
                std::process::id()
            ));
            {
                let file = File::create(&path)?;
                let group = file.create_group("X")?;
                attr_str(&group, "encoding-type", encoding)?;
                attr_str(&group, "encoding-version", "0.1.0")?;
                group
                    .new_attr::<i64>()
                    .shape(2)
                    .create("shape")?
                    .write(&[2, 3])?;
                // Both representations deliberately store indices out of order.
                let (ptr, idx, values): (Vec<i64>, Vec<i32>, Vec<f32>) = if encoding == "csr_matrix"
                {
                    (vec![0, 2, 3], vec![2, 0, 1], vec![5., 7., 11.])
                } else {
                    (vec![0, 2, 3, 3], vec![1, 0, 0], vec![5., 7., 11.])
                };
                group
                    .new_dataset_builder()
                    .with_data(&ptr)
                    .create("indptr")?;
                group
                    .new_dataset_builder()
                    .with_data(&idx)
                    .create("indices")?;
                group
                    .new_dataset_builder()
                    .with_data(&values)
                    .create("data")?;
            }
            let file = File::open(&path)?;
            let mut warnings = Vec::new();
            let matrix = read_matrix(&file, "X", &mut warnings)?.expect("sparse X");
            let Matrix::CsrF32(matrix) = matrix else {
                panic!("expected CSR")
            };
            let dense = matrix.to_dense();
            if encoding == "csr_matrix" {
                assert_eq!(dense[[0, 0]], 7.);
                assert_eq!(dense[[0, 2]], 5.);
                assert_eq!(dense[[1, 1]], 11.);
            } else {
                assert_eq!(dense[[0, 0]], 7.);
                assert_eq!(dense[[1, 0]], 5.);
                assert_eq!(dense[[0, 1]], 11.);
            }
            std::fs::remove_file(path)?;
        }
        Ok(())
    }

    #[test]
    fn roundtrip_minimal() -> Result<()> {
        let x = CsMat::new((2, 3), vec![0, 2, 3], vec![0, 2, 1], vec![1.0, 2.0, 3.0]);
        let mut obs = Table::new("_index", vec!["c1".into(), "c2".into()]);
        obs.add_column("cluster", Column::String(vec!["a".into(), "b".into()]))?;
        let var = Table::new("_index", vec!["g1".into(), "g2".into(), "g3".into()]);
        let a = AnnData {
            x: Matrix::CsrF32(x),
            obs,
            var,
            obsm: BTreeMap::new(),
            layers: BTreeMap::new(),
        };
        let p = std::env::temp_dir().join(format!("sc_anndata_{}.h5ad", std::process::id()));
        a.write_h5ad(&p)?;
        let b = AnnData::read_h5ad(&p)?;
        assert_eq!(b.data.x.shape(), (2, 3));
        std::fs::remove_file(p)?;
        Ok(())
    }
}

#[cfg(test)]
mod categorical_roundtrip_tests {
    use super::*;
    #[test]
    fn categorical_roundtrip_preserves_order_and_missing() -> Result<()> {
        let path = std::env::temp_dir().join(format!("lumrik_factor_{}.h5ad", std::process::id()));
        let mut obs = Table::new("_index", vec!["c1".into(), "c2".into(), "c3".into()]);
        obs.add_column(
            "cell_type",
            Column::Categorical {
                categories: vec!["HSC".into(), "MEP".into(), "Erythroid".into()],
                codes: vec![2, 0, -1],
                ordered: true,
            },
        )?;
        let var = Table::new("_index", vec!["TP53".into()]);
        let data = AnnData {
            x: Matrix::DenseF32(Array2::zeros((3, 1))),
            obs,
            var,
            obsm: BTreeMap::new(),
            layers: BTreeMap::new(),
        };
        data.write_h5ad(&path)?;
        let read = AnnData::read_h5ad(&path)?;
        let (factor, values, ordered) = read.data.obs.factor("cell_type")?.expect("factor");
        assert!(ordered);
        assert_eq!(factor.get_levels(), &["HSC", "MEP", "Erythroid"]);
        assert_eq!(values[0], 2.0);
        assert_eq!(values[1], 0.0);
        assert!(values[2].is_nan());
        std::fs::remove_file(path)?;
        Ok(())
    }
}
