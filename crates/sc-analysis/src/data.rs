use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use sprs::CsMat;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct CellAnnotations {
    headers: Vec<String>,
    columns: HashMap<String, Vec<String>>,
}

impl CellAnnotations {
    pub fn headers(&self) -> &[String] {
        &self.headers
    }
    pub fn column(&self, name: &str) -> Option<&[String]> {
        self.columns.get(name).map(Vec::as_slice)
    }
    pub fn groups(&self, name: &str) -> Result<BTreeMap<String, Vec<usize>>> {
        let values = self
            .column(name)
            .with_context(|| format!("cell annotation column {name:?} does not exist"))?;
        let mut groups = BTreeMap::<String, Vec<usize>>::new();
        for (cell_idx, value) in values.iter().enumerate() {
            groups.entry(value.clone()).or_default().push(cell_idx);
        }
        Ok(groups)
    }
    pub(crate) fn from_columns(columns: Vec<(&str, Vec<String>)>, n_cells: usize) -> Result<Self> {
        let mut headers = Vec::new();
        let mut map = HashMap::new();
        for (name, values) in columns {
            if values.len() != n_cells {
                bail!(
                    "annotation {name:?} has {} values for {n_cells} cells",
                    values.len()
                );
            }
            headers.push(name.to_string());
            map.insert(name.to_string(), values);
        }
        Ok(Self {
            headers,
            columns: map,
        })
    }
}

#[derive(Debug, Clone)]
pub struct SingleCellData {
    /// Sparse feature x cell matrix, CSR for efficient feature access.
    pub matrix: CsMat<f32>,
    pub features: Vec<String>,
    pub cells: Vec<String>,
    pub annotations: CellAnnotations,
    feature_lookup: HashMap<String, usize>,
}

impl SingleCellData {
    pub(crate) fn new(
        matrix: CsMat<f32>,
        features: Vec<String>,
        cells: Vec<String>,
        annotations: CellAnnotations,
    ) -> Result<Self> {
        if matrix.rows() != features.len() || matrix.cols() != cells.len() {
            bail!("matrix dimensions do not match features/cells");
        }
        let mut feature_lookup = HashMap::with_capacity(features.len());
        for (idx, feature) in features.iter().enumerate() {
            feature_lookup.entry(feature.clone()).or_insert(idx);
        }
        Ok(Self {
            matrix,
            features,
            cells,
            annotations,
            feature_lookup,
        })
    }

    pub fn from_mtx_dir(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        Self::from_mtx(
            path.join("expression.mtx.gz"),
            path.join("genes.tsv.gz"),
            path.join("cells.tsv.gz"),
        )
    }
    pub fn from_mtx(
        matrix: impl AsRef<Path>,
        features: impl AsRef<Path>,
        cells: impl AsRef<Path>,
    ) -> Result<Self> {
        let features = read_lines(features.as_ref())?;
        let (cell_ids, annotations) = read_cells(cells.as_ref())?;
        let matrix = read_matrix_market_csc(matrix.as_ref())?.to_csr();
        Self::new(matrix, features, cell_ids, annotations)
    }
    pub fn n_features(&self) -> usize {
        self.features.len()
    }
    pub fn n_cells(&self) -> usize {
        self.cells.len()
    }
    pub fn feature_index(&self, name: &str) -> Option<usize> {
        self.feature_lookup.get(name).copied()
    }
    pub fn detection_fraction(&self, feature_idx: usize, cells: &[usize]) -> Result<f64> {
        if cells.is_empty() {
            return Ok(0.0);
        }
        let row = self
            .matrix
            .outer_view(feature_idx)
            .with_context(|| format!("feature index {feature_idx} is out of range"))?;
        let wanted: HashSet<usize> = cells.iter().copied().collect();
        let detected = row
            .iter()
            .filter(|(cell_idx, value)| **value != 0.0 && wanted.contains(cell_idx))
            .count();
        Ok(detected as f64 / cells.len() as f64)
    }
}

pub(crate) fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}
pub(crate) fn read_lines(path: &Path) -> Result<Vec<String>> {
    open_reader(path)?
        .lines()
        .map(|line| line.with_context(|| format!("reading {}", path.display())))
        .collect()
}
fn read_cells(path: &Path) -> Result<(Vec<String>, CellAnnotations)> {
    let reader: Box<dyn Read> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(File::open(path)?))
    } else {
        Box::new(File::open(path)?)
    };
    let mut csv = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    let headers = csv.headers()?.iter().map(str::to_owned).collect::<Vec<_>>();
    let cell_col = headers
        .iter()
        .position(|x| x == "cell")
        .context("cells.tsv has no 'cell' column")?;
    let mut columns = headers
        .iter()
        .map(|h| (h.clone(), Vec::new()))
        .collect::<HashMap<_, _>>();
    let mut cell_ids = Vec::new();
    for record in csv.records() {
        let record = record?;
        cell_ids.push(record.get(cell_col).unwrap_or_default().to_owned());
        for (idx, header) in headers.iter().enumerate() {
            columns
                .get_mut(header)
                .unwrap()
                .push(record.get(idx).unwrap_or_default().to_owned());
        }
    }
    Ok((cell_ids, CellAnnotations { headers, columns }))
}

pub(crate) fn read_matrix_market_csc(path: &Path) -> Result<CsMat<f32>> {
    let mut lines = open_reader(path)?.lines();
    let banner = lines.next().context("empty Matrix Market file")??;
    if !banner.starts_with("%%MatrixMarket matrix coordinate") {
        bail!("unsupported Matrix Market banner: {banner}");
    }
    let dimensions = loop {
        let line = lines.next().context("Matrix Market dimensions missing")??;
        if !line.starts_with('%') && !line.trim().is_empty() {
            break line;
        }
    };
    let dims = dimensions.split_whitespace().collect::<Vec<_>>();
    if dims.len() != 3 {
        bail!("invalid Matrix Market dimensions: {dimensions}");
    }
    let rows: usize = dims[0].parse()?;
    let cols: usize = dims[1].parse()?;
    let nnz: usize = dims[2].parse()?;
    let mut indptr = Vec::with_capacity(cols + 1);
    let mut indices = Vec::with_capacity(nnz);
    let mut values = Vec::with_capacity(nnz);
    indptr.push(0);
    let mut current_col = 0usize;
    let mut previous_col = 0usize;
    for line in lines {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('%') {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 3 {
            bail!("invalid Matrix Market entry: {line}");
        }
        let row1: usize = fields[0].parse()?;
        let col1: usize = fields[1].parse()?;
        let value: f32 = fields[2].parse()?;
        if row1 == 0 || row1 > rows || col1 == 0 || col1 > cols {
            bail!("Matrix Market index out of range: {line}");
        }
        let col = col1 - 1;
        if !indices.is_empty() && col < previous_col {
            bail!(
                "Matrix Market entries are not column-major; expected Lumrik/10x column-major MEX output"
            );
        }
        while current_col < col {
            indptr.push(indices.len());
            current_col += 1;
        }
        indices.push(row1 - 1);
        values.push(value);
        previous_col = col;
    }
    while indptr.len() < cols + 1 {
        indptr.push(indices.len());
    }
    if indices.len() != nnz {
        bail!(
            "Matrix Market declared {nnz} entries but {} were read",
            indices.len()
        );
    }
    Ok(CsMat::new_csc((rows, cols), indptr, indices, values))
}
