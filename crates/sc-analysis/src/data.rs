use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use sprs::{CsMat, TriMat};
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

    /// Load sc-analysis' own pre-analysed export format.
    ///
    /// This is intentionally separate from the external 10x/BD MEX contract,
    /// which is owned by `scdata`.
    pub fn from_mtx_dir(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let features = Self::read_lines(&path.join("genes.tsv.gz"))?;
        let (cells, annotations) = Self::read_cells(&path.join("cells.tsv.gz"))?;
        let matrix = Self::read_matrix_market(&path.join("expression.mtx.gz"))?;
        Self::new(matrix, features, cells, annotations)
    }

    fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        if path.extension().and_then(|x| x.to_str()) == Some("gz") {
            Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }

    fn read_lines(path: &Path) -> Result<Vec<String>> {
        Self::open_reader(path)?
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
        let mut csv = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(reader);
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
                    .expect("annotation header initialized above")
                    .push(record.get(idx).unwrap_or_default().to_owned());
            }
        }
        Ok((cell_ids, CellAnnotations { headers, columns }))
    }

    fn read_matrix_market(path: &Path) -> Result<CsMat<f32>> {
        let mut lines = Self::open_reader(path)?.lines();
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
        let mut triplets = TriMat::<f32>::with_capacity((rows, cols), nnz);
        for line in lines {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('%') {
                continue;
            }
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() < 3 {
                bail!("invalid Matrix Market entry: {line}");
            }
            let row: usize = fields[0].parse()?;
            let col: usize = fields[1].parse()?;
            let value: f32 = fields[2].parse()?;
            if row == 0 || row > rows || col == 0 || col > cols {
                bail!("Matrix Market index out of range: {line}");
            }
            triplets.add_triplet(row - 1, col - 1, value);
        }
        if triplets.nnz() != nnz {
            bail!(
                "Matrix Market declared {nnz} entries but {} were read",
                triplets.nnz()
            );
        }
        Ok(triplets.to_csr())
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

