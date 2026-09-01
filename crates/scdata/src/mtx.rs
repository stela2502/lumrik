use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use int_to_str::IntToStr;
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
        let path = dir.as_ref().join("features.tsv.gz");
        let reader = gz_lines(&path)?;

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

    let matrix_path = dir.join("matrix.mtx.gz");
    let mut lines = gz_lines(&matrix_path)?;
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

    Ok((cells, index, cell_barcode_len))
}

pub fn read_mtx_cell_ids(dir: impl AsRef<Path>) -> Result<HashSet<u64>> {
    let dir = dir.as_ref();
    Ok(read_mtx_barcodes(dir)?
        .into_iter()
        .map(|(_, id)| id)
        .collect())
}

fn read_mtx_barcodes(dir: &Path) -> Result<Vec<(String, u64)>> {
    let path = dir.join("barcodes.tsv.gz");
    let reader = gz_lines(&path)?;
    let mut out = Vec::new();

    for line in reader {
        let barcode = line.with_context(|| format!("reading {}", path.display()))?;
        let barcode = barcode.trim().to_string();
        if barcode.is_empty() {
            continue;
        }

        let sequence = barcode
            .split_once('-')
            .map(|(seq, _)| seq)
            .unwrap_or(&barcode);

        let id = IntToStr::new(sequence.as_bytes()).into_u64();
        out.push((barcode, id));
    }

    Ok(out)
}

fn gz_lines(path: &Path) -> Result<impl Iterator<Item = std::io::Result<String>>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let decoder = GzDecoder::new(file);
    Ok(BufReader::new(decoder).lines())
}
