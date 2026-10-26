//! Inspect a real H5AD dataset before requesting AI interpretation.
use anyhow::{Context, Result, bail};
use sc_anndata::{AnnData, Column, DatasetOverview};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, PartialEq, Eq)]
pub struct DatasetSummary {
    pub cells: usize,
    pub genes: usize,
    pub grouping: String,
    pub groups: BTreeMap<String, usize>,
    pub available_annotations: Vec<String>,
    pub warnings: Vec<String>,
}

/// Inspect an H5AD without assuming that any particular cell annotation exists.
pub fn inspect_overview(path: impl AsRef<Path>) -> Result<(DatasetOverview, Vec<String>)> {
    let read = AnnData::read_h5ad(path.as_ref())
        .with_context(|| format!("reading {}", path.as_ref().display()))?;
    let overview = read.data.overview(40)?;
    let warnings = read
        .warnings
        .into_iter()
        .map(|w| format!("{}: {}", w.path, w.reason))
        .collect();
    Ok((overview, warnings))
}

pub fn inspect(path: impl AsRef<Path>, grouping: &str) -> Result<DatasetSummary> {
    let read = AnnData::read_h5ad(path.as_ref())
        .with_context(|| format!("reading {}", path.as_ref().display()))?;
    summarize(
        &read.data,
        grouping,
        read.warnings
            .into_iter()
            .map(|w| format!("{}: {}", w.path, w.reason))
            .collect(),
    )
}

pub fn summarize(data: &AnnData, grouping: &str, warnings: Vec<String>) -> Result<DatasetSummary> {
    data.validate()?;
    let column = data.obs.columns.get(grouping).with_context(|| {
        format!(
            "grouping {grouping:?} missing; available: {:?}",
            data.obs.columns.keys().collect::<Vec<_>>()
        )
    })?;
    let mut groups = BTreeMap::new();
    if matches!(column, Column::F64(_)) {
        bail!("grouping {grouping:?} is continuous; select a categorical annotation");
    }
    for row in 0..column.len() {
        let label = column.get_string(row).unwrap_or_else(|| "NA".to_owned());
        *groups.entry(label).or_insert(0) += 1;
    }
    let (cells, genes) = data.x.shape();
    Ok(DatasetSummary {
        cells,
        genes,
        grouping: grouping.into(),
        groups,
        available_annotations: data.obs.columns.keys().cloned().collect(),
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use sc_anndata::{Matrix, Table};
    #[test]
    fn grouping_counts_and_dimensions_are_exact() {
        let mut obs = Table::new("cell", vec!["c1".into(), "c2".into(), "c3".into()]);
        obs.add_column(
            "celltype",
            Column::String(vec![
                "erythroid".into(),
                "myeloid".into(),
                "erythroid".into(),
            ]),
        )
        .unwrap();
        let var = Table::new("gene", vec!["TP53".into(), "CDKN1A".into()]);
        let data = AnnData {
            x: Matrix::DenseF32(ndarray::Array2::zeros((3, 2))),
            obs,
            var,
            obsm: Default::default(),
            layers: Default::default(),
        };
        let result = summarize(&data, "celltype", vec![]).unwrap();
        assert_eq!((result.cells, result.genes), (3, 2));
        assert_eq!(result.groups["erythroid"], 2);
        assert_eq!(result.groups["myeloid"], 1);
        assert!(summarize(&data, "missing", vec![]).is_err());
    }
}

#[cfg(test)]
mod categorical_integration_tests {
    use super::*;
    use sc_anndata::{Matrix, Table};
    #[test]
    fn native_celltype_factor_is_groupable() -> Result<()> {
        let mut obs = Table::new("_index", vec!["c1".into(), "c2".into(), "c3".into()]);
        obs.add_column(
            "cell_type",
            Column::Categorical {
                categories: vec!["HSC".into(), "MEP".into(), "Erythroid".into()],
                codes: vec![2, 0, -1],
                ordered: true,
            },
        )?;
        let data = AnnData {
            x: Matrix::DenseF32(ndarray::Array2::zeros((3, 1))),
            obs,
            var: Table::new("_index", vec!["TP53".into()]),
            obsm: Default::default(),
            layers: Default::default(),
        };
        let result = summarize(&data, "cell_type", vec![])?;
        assert_eq!(result.groups["Erythroid"], 1);
        assert_eq!(result.groups["HSC"], 1);
        assert_eq!(result.groups["NA"], 1);
        Ok(())
    }
}
