use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use int_to_dna::IntToDna;
use mapping_info::MappingInfo;

use crate::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MexFeature {
    pub source_row: usize,
    pub id: String,
    pub name: String,
    pub feature_type: String,
}

#[derive(Debug, Clone)]
pub struct MexFeatureIndex {
    features: Vec<MexFeature>,
    name_to_id: HashMap<String, u64>,
    source_row_to_feature: HashMap<usize, u64>,
}

impl MexFeatureIndex {
    pub fn from_dir(dir: impl AsRef<Path>, feature_type: &str) -> Result<Self> {
        let dir = dir.as_ref();
        let path = find_text_file(dir, &["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])?;
        let reader = text_lines(&path)?;

        let mut features = Vec::new();
        let mut name_to_id = HashMap::new();
        let mut source_row_to_feature = HashMap::new();

        for (source_row, line) in reader.enumerate() {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                bail!(
                    "Malformed 10x feature line {} in {}",
                    source_row + 1,
                    path.display()
                );
            }

            let this_type = fields.get(2).copied().unwrap_or("Gene Expression");
            if this_type != feature_type {
                continue;
            }

            let feature_id = features.len() as u64;
            let feature = MexFeature {
                source_row,
                id: fields[0].to_string(),
                name: fields[1].to_string(),
                feature_type: this_type.to_string(),
            };

            source_row_to_feature.insert(source_row, feature_id);
            name_to_id.insert(feature.id.clone(), feature_id);
            name_to_id.entry(feature.name.clone()).or_insert(feature_id);
            features.push(feature);
        }

        if features.is_empty() {
            bail!(
                "No features of type {:?} found in {}",
                feature_type,
                path.display()
            );
        }

        Ok(Self {
            features,
            name_to_id,
            source_row_to_feature,
        })
    }

    pub fn features(&self) -> &[MexFeature] {
        &self.features
    }

    pub fn feature_for_source_row(&self, zero_based_row: usize) -> Option<u64> {
        self.source_row_to_feature.get(&zero_based_row).copied()
    }

    pub fn validate_compatible(&self, other: &Self) -> Result<()> {
        if self.features.len() != other.features.len() {
            bail!(
                "10x matrices contain different selected feature counts: {} != {}",
                self.features.len(),
                other.features.len()
            );
        }

        for (a, b) in self.features.iter().zip(&other.features) {
            if a.id != b.id || a.name != b.name || a.feature_type != b.feature_type {
                bail!("10x feature definitions differ: {:?} != {:?}", a, b);
            }
        }

        Ok(())
    }
}

impl FeatureIndex for MexFeatureIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        &self.features[feature_id as usize].name
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let f = &self.features[feature_id as usize];
        format!("{}\t{}\t{}", f.id, f.name, f.feature_type)
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.features.len() as u64).collect()
    }
}

/// Load one selected 10x feature type into `Scdata`.
///
/// Returns the sparse counts, the selected feature index, and the barcode
/// sequence length. Matrix rows belonging to other feature types are ignored.
pub fn load_mtx_feature_matrix(
    dir: impl AsRef<Path>,
    feature_type: &str,
    threads: usize,
) -> Result<(Scdata, MexFeatureIndex, usize)> {
    let dir = dir.as_ref();
    let index = MexFeatureIndex::from_dir(dir, feature_type)?;
    let barcodes = read_mtx_barcodes(dir)?;

    let cell_barcode_len = barcodes
        .first()
        .map(|(barcode, _)| {
            barcode
                .split_once('-')
                .map(|(seq, _)| seq)
                .unwrap_or(barcode)
                .len()
        })
        .unwrap_or(0);

    let mut cells = Scdata::new(threads.max(1), MatrixValueType::Integer);
    let mut report = MappingInfo::new(None, 0.0, 0);

    let matrix_path = find_text_file(dir, &["matrix.mtx.gz", "matrix.mtx"])?;
    let mut lines = text_lines(&matrix_path)?;
    let mut header_seen = false;
    let mut dims_seen = false;

    for line in lines.by_ref() {
        let line = line.with_context(|| format!("reading {}", matrix_path.display()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if !header_seen {
            header_seen = true;
            if !trimmed.starts_with("%%MatrixMarket") {
                bail!("{} is not a MatrixMarket file", matrix_path.display());
            }
            continue;
        }

        if trimmed.starts_with('%') {
            continue;
        }

        if !dims_seen {
            let dims: Vec<_> = trimmed.split_whitespace().collect();
            if dims.len() < 3 {
                bail!(
                    "Malformed MatrixMarket dimensions in {}",
                    matrix_path.display()
                );
            }
            let n_cols: usize = dims[1].parse()?;
            if n_cols != barcodes.len() {
                bail!(
                    "Matrix columns ({n_cols}) do not match barcode count ({}) in {}",
                    barcodes.len(),
                    dir.display()
                );
            }
            dims_seen = true;
            continue;
        }

        let mut p = trimmed.split_whitespace();
        let row: usize = p.next().context("missing MatrixMarket row")?.parse()?;
        let col: usize = p.next().context("missing MatrixMarket column")?.parse()?;
        let value: f64 = p.next().context("missing MatrixMarket value")?.parse()?;

        if row == 0 || col == 0 || col > barcodes.len() {
            bail!("Out-of-range MatrixMarket coordinate: row={row}, col={col}");
        }
        if value < 0.0 || value.fract() != 0.0 {
            bail!("Feature count must be a non-negative integer, got {value}");
        }

        let Some(feature_id) = index.feature_for_source_row(row - 1) else {
            continue;
        };

        let count = value as u32;
        if count == 0 {
            continue;
        }

        let cell_id = barcodes[col - 1].1;
        cells.try_insert_value(
            &cell_id,
            GeneUmiHash(feature_id, 0),
            count as f32,
            &mut report,
        );
    }

    cells.set_matrix_order(
        index.ordered_feature_ids(),
        barcodes.iter().map(|(_, cell_id)| *cell_id).collect(),
    );

    Ok((cells, index, cell_barcode_len))
}

/// Return the packed cell ids from a 10x-style MEX export.
///
/// `path` may point either at the MEX directory itself or at a Nelrune
/// analysis directory containing `exonic/`. Both compressed and plain TSV
/// barcode files are accepted. Consumers should use this API instead of
/// interpreting `barcodes.tsv[.gz]` themselves.
pub fn read_mtx_cell_ids(path: impl AsRef<Path>) -> Result<HashSet<u64>> {
    let dir = resolve_mtx_dir(path.as_ref())?;
    Ok(read_mtx_barcodes(&dir)?
        .into_iter()
        .map(|(_, id)| id)
        .collect())
}

fn resolve_mtx_dir(path: &Path) -> Result<std::path::PathBuf> {
    if barcode_path(path).is_some() {
        return Ok(path.to_path_buf());
    }

    let exonic = path.join("exonic");
    if barcode_path(&exonic).is_some() {
        return Ok(exonic);
    }

    bail!(
        "{} is neither a MEX directory nor an analysis directory containing exonic/barcodes.tsv[.gz]",
        path.display()
    )
}

fn find_text_file(dir: &Path, names: &[&str]) -> Result<std::path::PathBuf> {
    names
        .iter()
        .map(|name| dir.join(name))
        .find(|path| path.is_file())
        .with_context(|| format!("none of {} found in {}", names.join(", "), dir.display()))
}

fn barcode_path(dir: &Path) -> Option<std::path::PathBuf> {
    ["barcodes.tsv.gz", "barcodes.tsv"]
        .into_iter()
        .map(|name| dir.join(name))
        .find(|path| path.is_file())
}

/// Return barcode labels together with stable internal cell ids in MEX column order.
///
/// DNA barcodes retain Lumrik's packed `IntToDna` identity so imported Lumrik/10x
/// matrices still match molecule-level cell ids. Numeric vendor ids are used as-is.
/// Any other external label receives a deterministic internal id while the original
/// label is preserved for analysis/reporting. Consumers must never decode these ids
/// back into barcode strings; the label in this return value is authoritative.
pub fn read_mtx_barcodes(dir: &Path) -> Result<Vec<(String, u64)>> {
    let path = barcode_path(dir).with_context(|| {
        format!(
            "no barcodes.tsv.gz or barcodes.tsv found in {}",
            dir.display()
        )
    })?;
    let reader = text_lines(&path)?;
    let mut out = Vec::new();
    let mut ids = HashMap::<u64, String>::new();

    for line in reader {
        let barcode = line.with_context(|| format!("reading {}", path.display()))?;
        let barcode = barcode.trim().to_string();
        if barcode.is_empty() {
            continue;
        }

        let id = external_cell_id(&barcode);
        if let Some(previous) = ids.insert(id, barcode.clone()) {
            if previous != barcode {
                bail!(
                    "cell labels {:?} and {:?} map to the same internal id {id}",
                    previous,
                    barcode
                );
            }
        }
        out.push((barcode, id));
    }

    Ok(out)
}

fn external_cell_id(label: &str) -> u64 {
    let sequence = label.split_once('-').map(|(seq, _)| seq).unwrap_or(label);

    if !sequence.is_empty()
        && sequence
            .bytes()
            .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        return IntToDna::new(sequence.as_bytes()).into_u64();
    }

    if let Ok(id) = label.parse::<u64>() {
        return id;
    }

    // FNV-1a gives arbitrary external labels a stable process-independent id.
    // Collisions are checked explicitly by `read_mtx_barcodes`.
    let mut hash = 0xcbf29ce484222325_u64;
    for byte in label.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x100000001b3);
    }
    hash
}

fn text_lines(path: &Path) -> Result<Box<dyn Iterator<Item = std::io::Result<String>>>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file)).lines()))
    } else {
        Ok(Box::new(BufReader::new(file).lines()))
    }
}
