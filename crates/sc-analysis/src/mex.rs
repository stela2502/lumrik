use crate::data::{CellAnnotations, SingleCellData, open_reader, read_matrix_market_csc};
use anyhow::{Context, Result, bail};
use std::io::BufRead;
use std::path::{Path, PathBuf};

fn find(dir: &Path, names: &[&str]) -> Result<PathBuf> {
    names
        .iter()
        .map(|n| dir.join(n))
        .find(|p| p.is_file())
        .with_context(|| format!("none of {} found in {}", names.join(", "), dir.display()))
}

pub(crate) fn load_mex(dir: &Path) -> Result<SingleCellData> {
    let matrix_path = find(dir, &["matrix.mtx.gz", "matrix.mtx"])?;
    let features_path = find(
        dir,
        &[
            "features.tsv.gz",
            "features.tsv",
            "genes.tsv.gz",
            "genes.tsv",
        ],
    )?;
    let barcodes_path = find(dir, &["barcodes.tsv.gz", "barcodes.tsv"])?;
    let matrix = read_matrix_market_csc(&matrix_path)?.to_csr();

    let mut features = Vec::new();
    for line in open_reader(&features_path)?.lines() {
        let line = line?;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.is_empty() {
            continue;
        }
        // Standard MEX: id, symbol, feature type. Prefer the human-readable symbol.
        features.push(fields.get(1).copied().unwrap_or(fields[0]).to_string());
    }
    let cells = open_reader(&barcodes_path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?;
    if matrix.rows() != features.len() {
        bail!(
            "MEX matrix rows {} != feature rows {}",
            matrix.rows(),
            features.len()
        );
    }
    if matrix.cols() != cells.len() {
        bail!(
            "MEX matrix columns {} != barcodes {}",
            matrix.cols(),
            cells.len()
        );
    }
    let annotations = CellAnnotations::from_columns(vec![("cell", cells.clone())], cells.len())?;
    SingleCellData::new(matrix, features, cells, annotations)
}
