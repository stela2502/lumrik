use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use mapping_info::MappingInfo;
use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};
use stela_int_to_str::IntToStr;

use crate::dataset::{GuideDataset, GuideObservation};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GuideFeature {
    pub source_row: usize,
    pub id: String,
    pub name: String,
    pub feature_type: String,
}

#[derive(Debug, Clone)]
pub struct GuideFeatureIndex {
    guides: Vec<GuideFeature>,
    name_to_id: HashMap<String, u64>,
    source_row_to_guide: HashMap<usize, u32>,
}

impl GuideFeatureIndex {
    pub fn from_10x_dir(dir: impl AsRef<Path>, feature_type: &str) -> Result<Self> {
        let path = dir.as_ref().join("features.tsv.gz");
        let reader = gz_lines(&path)?;

        let mut guides = Vec::new();
        let mut name_to_id = HashMap::new();
        let mut source_row_to_guide = HashMap::new();

        for (source_row, line) in reader.enumerate() {
            let line = line.with_context(|| format!("reading {}", path.display()))?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 2 {
                bail!("Malformed 10x feature line {} in {}", source_row + 1, path.display());
            }

            let this_type = fields.get(2).copied().unwrap_or("Gene Expression");
            if this_type != feature_type {
                continue;
            }

            let local_id = guides.len() as u32;
            let feature = GuideFeature {
                source_row,
                id: fields[0].to_string(),
                name: fields[1].to_string(),
                feature_type: this_type.to_string(),
            };

            // 10x matrices index rows by the full features.tsv row.
            source_row_to_guide.insert(source_row, local_id);

            // Prefer the feature id, but also accept the display name.
            name_to_id.insert(feature.id.clone(), local_id as u64);
            name_to_id.entry(feature.name.clone()).or_insert(local_id as u64);
            guides.push(feature);
        }

        if guides.is_empty() {
            bail!(
                "No features of type {:?} found in {}",
                feature_type,
                path.display()
            );
        }

        Ok(Self {
            guides,
            name_to_id,
            source_row_to_guide,
        })
    }

    pub fn guides(&self) -> &[GuideFeature] {
        &self.guides
    }

    pub fn guide_for_source_row(&self, zero_based_row: usize) -> Option<u32> {
        self.source_row_to_guide.get(&zero_based_row).copied()
    }

    pub fn validate_compatible(&self, other: &Self) -> Result<()> {
        if self.guides.len() != other.guides.len() {
            bail!(
                "Raw and filtered matrices contain different guide counts: {} != {}",
                self.guides.len(),
                other.guides.len()
            );
        }

        for (a, b) in self.guides.iter().zip(&other.guides) {
            if a.id != b.id || a.name != b.name || a.feature_type != b.feature_type {
                bail!(
                    "Raw/filtered guide definitions differ: {:?} != {:?}",
                    a,
                    b
                );
            }
        }
        Ok(())
    }
}

impl FeatureIndex for GuideFeatureIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        &self.guides[feature_id as usize].name
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let f = &self.guides[feature_id as usize];
        format!("{}\t{}\t{}", f.id, f.name, f.feature_type)
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.guides.len() as u64).collect()
    }
}

#[derive(Debug, Clone)]
pub struct TenxGuideInput {
    pub raw_dir: PathBuf,
    pub filtered_dir: PathBuf,
    pub feature_type: String,
    pub threads: usize,
}

impl TenxGuideInput {
    pub fn new(raw_dir: PathBuf, filtered_dir: PathBuf) -> Self {
        Self {
            raw_dir,
            filtered_dir,
            feature_type: "CRISPR Guide Capture".to_string(),
            threads: 1,
        }
    }

    pub fn indexes(&self) -> Result<(GuideFeatureIndex, GuideFeatureIndex)> {
        let raw = GuideFeatureIndex::from_10x_dir(&self.raw_dir, &self.feature_type)?;
        let filtered = GuideFeatureIndex::from_10x_dir(&self.filtered_dir, &self.feature_type)?;
        raw.validate_compatible(&filtered)?;
        Ok((raw, filtered))
    }

    pub fn filtered_cell_ids(&self) -> Result<HashSet<u64>> {
        Ok(read_barcodes(&self.filtered_dir)?
            .into_iter()
            .map(|(_, id)| id)
            .collect())
    }

    /// Load only the raw droplets that are NOT present in the filtered matrix.
    /// This is the empirical ambient/background population.
    pub fn load_background(&self, index: &GuideFeatureIndex) -> Result<GuideDataset> {
        let filtered = self.filtered_cell_ids()?;
        load_guide_dataset(
            &self.raw_dir,
            index,
            self.threads,
            |cell_id| !filtered.contains(&cell_id),
        )
    }

    /// Load the called-cell matrix after the ambient model has been fitted.
    pub fn load_filtered(&self, index: &GuideFeatureIndex) -> Result<GuideDataset> {
        load_guide_dataset(&self.filtered_dir, index, self.threads, |_| true)
    }
}

fn load_guide_dataset<F>(
    dir: &Path,
    index: &GuideFeatureIndex,
    threads: usize,
    keep_cell: F,
) -> Result<GuideDataset>
where
    F: Fn(u64) -> bool,
{
    let barcodes = read_barcodes(dir)?;
    let selected: Vec<bool> = barcodes.iter().map(|(_, id)| keep_cell(*id)).collect();

    let mut cell_ids = Vec::new();
    let mut barcode_by_id = HashMap::new();
    for ((barcode, id), keep) in barcodes.iter().zip(&selected) {
        if *keep {
            cell_ids.push(*id);
            barcode_by_id.insert(*id, barcode.clone());
        }
    }

    let mut cells = Scdata::new(threads.max(1), MatrixValueType::Integer);
    let mut by_guide = vec![Vec::new(); index.guides.len()];
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
                bail!("Malformed MatrixMarket dimensions in {}", matrix_path.display());
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
            bail!("Guide count must be a non-negative integer, got {value}");
        }
        if !selected[col - 1] {
            continue;
        }

        let Some(guide_id) = index.guide_for_source_row(row - 1) else {
            continue;
        };

        let count = value as u32;
        if count == 0 {
            continue;
        }

        let cell_id = barcodes[col - 1].1;

        // One MatrixMarket entry already represents the aggregated UMI count.
        // A fixed synthetic UMI is therefore enough: scdata stores the supplied
        // numeric value in CellData::total_reads.
        cells.try_insert_value(
            &cell_id,
            GeneUmiHash(guide_id as u64, 0),
            count as f32,
            &mut report,
        );

        by_guide[guide_id as usize].push(GuideObservation {
            cell_id,
            guide_id,
            count,
        });
    }

    Ok(GuideDataset {
        cells,
        by_guide,
        cell_ids,
        barcode_by_id,
    })
}

fn read_barcodes(dir: &Path) -> Result<Vec<(String, u64)>> {
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
