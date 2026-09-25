// thrust_data.rs
use std::fs;

//use gtf_splice_index::{GeneId, RefBlock, SpliceIndex, Strand, TranscriptId};
//use snp_index::{Genome, SnpIndex, VcfReadOptions};

//use rand::rngs::SmallRng;
//use rand::{Rng, SeedableRng};

use std::collections::{HashMap, HashSet};
use std::path::Path;

use crate::{FeatureIndex, MatrixValueType, Scdata};

use mapping_info::MappingInfo;

use crate::sparse_matrix::scdata;

pub struct QuantData {
    data: HashMap<String, Scdata>,
}

#[derive(Debug, Clone)]
pub struct CellAccounting {
    pub exonic_cells: usize,
    pub intronic_cells: usize,
    pub exonic_or_intronic_cells: usize,
    pub exonic_umis_before_filter: usize,
    pub exonic_umis_after_filter: usize,
    pub exonic_thresholds: Vec<(usize, usize)>,
}

impl std::fmt::Display for CellAccounting {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Cell accounting")?;
        writeln!(f, "---------------")?;
        writeln!(f, "exonic cells before cutoff: {}", self.exonic_cells)?;
        writeln!(f, "intronic cells before cutoff: {}", self.intronic_cells)?;
        writeln!(
            f,
            "cells with exonic or intronic evidence: {}",
            self.exonic_or_intronic_cells
        )?;
        writeln!(
            f,
            "unique exonic UMIs before filter: {}",
            self.exonic_umis_before_filter
        )?;
        writeln!(
            f,
            "unique exonic UMIs after filter: {}",
            self.exonic_umis_after_filter
        )?;
        for (threshold, count) in &self.exonic_thresholds {
            writeln!(f, "exonic cells with >= {threshold} UMIs: {count}")?;
        }
        Ok(())
    }
}

impl Default for QuantData {
    fn default() -> Self {
        Self::standard()
    }
}

impl QuantData {
    pub const EXONIC: &'static str = "exonic";
    pub const INTRONIC: &'static str = "intronic";
    pub const SNP_REF: &'static str = "ref";
    pub const SNP_ALT: &'static str = "alt";
    pub const STANDARD_DATASETS: &'static [&'static str] = &[
        Self::EXONIC,
        Self::INTRONIC,
        Self::SNP_REF,
        Self::SNP_ALT,
    ];

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.data.keys().map(String::as_str)
    }

    /// Non-destructive cell accounting before the canonical GEX cutoff.
    pub fn cell_accounting(&self) -> CellAccounting {
        const THRESHOLDS: [usize; 7] = [1, 10, 50, 100, 200, 400, 1_000];
        let exonic = self.data(Self::EXONIC).cell_ids();
        let intronic = self.data(Self::INTRONIC).cell_ids();
        let mut union = exonic.clone();
        union.extend(intronic.iter().copied());
        CellAccounting {
            exonic_cells: exonic.len(),
            intronic_cells: intronic.len(),
            exonic_or_intronic_cells: union.len(),
            exonic_umis_before_filter: self.data(Self::EXONIC).total_umis(),
            exonic_umis_after_filter: self.data(Self::EXONIC).total_umis(),
            exonic_thresholds: THRESHOLDS
                .into_iter()
                .map(|threshold| (threshold, self.data(Self::EXONIC).cells_with_min_umis(threshold)))
                .collect(),
        }
    }

    /// Return exonic cell ids with at least `min_umis` unique UMIs.
    ///
    /// This is the explicit fixed-cutoff counterpart to sc-beacon cell calling.
    pub fn cells_with_min_exonic_umis(&self, min_umis: usize) -> HashSet<u64> {
        self.data(Self::EXONIC)
            .cell_umi_counts()
            .into_iter()
            .filter_map(|(cell_id, umi_count)| {
                ((umi_count as usize) >= min_umis).then_some(cell_id)
            })
            .collect()
    }

    /// Apply one caller-selected cell set to every owned matrix.
    ///
    /// The caller owns the biological decision of which cells are canonical.
    /// QuantData owns applying that decision consistently to all matrices.
    pub fn finalize_for_cells(
        &mut self,
        keep: &HashSet<u64>,
        indexes: &HashMap<String, &dyn FeatureIndex>,
    ) -> Result<(), String> {
        for (name, data) in &mut self.data {
            if data.is_empty() {
                continue;
            }
            let index = indexes
                .get(name)
                .copied()
                .ok_or_else(|| format!("missing feature index for QuantData dataset '{name}'"))?;
            data.finalize_for_cells(keep, index);
        }
        Ok(())
    }

    /// Write every owned matrix using the caller-provided biological row/index
    /// definition for that named matrix.
    ///
    /// QuantData owns the output orchestration. Callers do not write its
    /// individual Scdata objects themselves.
    pub fn write_sparse<P: AsRef<Path>>(
        &mut self,
        base: P,
        indexes: &HashMap<String, &dyn FeatureIndex>,
    ) -> Result<(), String> {
        self.write_sparse_impl(base.as_ref(), indexes, None)
    }

    /// As `write_sparse`, but render cell barcodes using the known biological
    /// barcode length.
    pub fn write_sparse_with_cell_len<P: AsRef<Path>>(
        &mut self,
        base: P,
        indexes: &HashMap<String, &dyn FeatureIndex>,
        cell_barcode_len: usize,
    ) -> Result<(), String> {
        self.write_sparse_impl(base.as_ref(), indexes, Some(cell_barcode_len))
    }

    /// Write the complete sparse stack and a caller-selected filtered view.
    ///
    /// Each Scdata owns its own observed barcode and feature universe. The raw
    /// view therefore contains exactly the cells/features present in that
    /// dataset. The filtered view applies the same allowed-cell set to every
    /// dataset, while Scdata naturally keeps only cells it actually contains.
    pub fn write_raw_and_filtered_for_cells<P: AsRef<Path>>(
        &mut self,
        base: P,
        keep: &HashSet<u64>,
        indexes: &HashMap<String, &dyn FeatureIndex>,
        cell_barcode_len: Option<usize>,
    ) -> Result<CellAccounting, String> {
        let accounting = self.cell_accounting();
        let base = base.as_ref();

        // Raw is written first because filtering is intentionally destructive.
        // Finalizing each Scdata against its own cells creates a valid sparse
        // matrix without manufacturing empty barcode columns from other data.
        for (name, data) in &mut self.data {
            if data.is_empty() {
                continue;
            }
            let index = indexes
                .get(name)
                .copied()
                .unwrap_or_else(|| panic!("missing feature index for QuantData dataset '{name}'"));
            let observed = data.cell_ids();
            data.finalize_for_cells(&observed, index);
        }
        self.write_sparse_impl(&base.join("raw"), indexes, cell_barcode_len)?;

        self.finalize_for_cells(keep, indexes)?;
        self.write_sparse_impl(&base.join("filtered"), indexes, cell_barcode_len)?;

        Ok(accounting)
    }

    /// Apply the caller-selected canonical cells and write all named matrices.
    pub fn write_sparse_for_cells<P: AsRef<Path>>(
        &mut self,
        base: P,
        keep: &HashSet<u64>,
        indexes: &HashMap<String, &dyn FeatureIndex>,
        cell_barcode_len: Option<usize>,
    ) -> Result<(), String> {
        self.finalize_for_cells(keep, indexes)?;
        self.write_sparse_impl(base.as_ref(), indexes, cell_barcode_len)
    }

    fn write_sparse_impl(
        &mut self,
        base: &Path,
        indexes: &HashMap<String, &dyn FeatureIndex>,
        cell_barcode_len: Option<usize>,
    ) -> Result<(), String> {
        fs::create_dir_all(base)
            .map_err(|e| format!("failed to create {}: {e}", base.display()))?;

        for (name, data) in &mut self.data {
            if data.is_empty() {
                continue;
            }
            let index = indexes
                .get(name)
                .copied()
                .ok_or_else(|| format!("missing feature index for QuantData dataset '{name}'"))?;
            let out = base.join(name);
            match cell_barcode_len {
                Some(len) => data
                    .write_sparse_with_cell_len(&out, index, len)
                    .map_err(|e| format!("writing QuantData dataset '{name}' failed: {e}"))?,
                None => data
                    .write_sparse(&out, index)
                    .map_err(|e| format!("writing QuantData dataset '{name}' failed: {e}"))?,
            };
        }
        Ok(())
    }

    pub fn new(names: &[&str]) -> Self {
        let mut data = HashMap::new();
        for name in names {
            if data
                .insert((*name).to_string(), Scdata::new(1, MatrixValueType::Real))
                .is_some()
            {
                panic!("QuantData dataset name '{name}' was supplied more than once");
            }
        }
        Self { data }
    }

    pub fn standard() -> Self {
        Self::new(Self::STANDARD_DATASETS)
    }

    /// Add a completed externally-produced dataset to this quantification bundle.
    ///
    /// Dataset names are unique identities inside QuantData. Supplying a name
    /// that is already present is a programmer error and therefore panics.
    pub fn add_dataset(&mut self, name: impl Into<String>, data: Scdata) {
        let name = name.into();
        if self.data.contains_key(&name) {
            panic!("QuantData already contains dataset '{name}'");
        }
        self.data.insert(name, data);
    }

    fn data(&self, name: &str) -> &Scdata {
        self.data.get(name)
            .unwrap_or_else(|| panic!("QuantData has no dataset named '{name}'"))
    }

    fn data_mut(&mut self, name: &str) -> &mut Scdata {
        self.data.get_mut(name)
            .unwrap_or_else(|| panic!("QuantData has no dataset named '{name}'"))
    }

    pub fn merge_named(&mut self, target: &str, source: &str) -> MappingInfo {
        if target == source {
            panic!("cannot merge QuantData dataset '{source}' into itself");
        }

        let source_data = self.data.remove(source)
            .unwrap_or_else(|| panic!("QuantData has no dataset named '{source}'"));
        let merge = self.data_mut(target).merge(&source_data);
        self.data.insert(source.to_string(), source_data);
        merge
    }

    /// Add one observation to a named matrix. Run accounting stays with the caller.
    pub fn try_insert(
        &mut self,
        name: &str,
        cell: &u64,
        value: crate::GeneUmiHash,
        weight: f32,
        report: &mut MappingInfo,
    ) -> crate::ScdataInsertState {
        let state = self.data_mut(name).try_insert(cell, value, weight);
        report.report(format!("{name}_{}", state.as_str()));
        report.ok_reads += 1;
        if state == crate::ScdataInsertState::Duplicate {
            report.pcr_duplicates += 1;
            report.local_dup += 1;
        }
        state
    }

    pub fn is_empty(&self, name: &str) -> bool { self.data(name).is_empty() }
    pub fn cell_ids(&self, name: &str) -> HashSet<u64> { self.data(name).cell_ids() }
    pub fn total_umis(&self, name: &str) -> usize { self.data(name).total_umis() }
    pub fn cell_umi_counts(&self, name: &str) -> Vec<(u64, u32)> { self.data(name).cell_umi_counts() }
    pub fn observed_feature_ids(&self, name: &str) -> HashSet<u64> { self.data(name).observed_feature_ids() }
    pub fn dimensions(&self, name: &str) -> (usize, usize, usize) { self.data(name).dimensions() }

    /// Read truth matrices from a truth output directory.
    ///
    /// Expected layout:
    ///
    /// - `<base>/exonic`
    /// - `<base>/intronic`
    /// - `<base>/ref`
    /// - `<base>/alt`
    pub fn from_path<P, G, S>(base: P, gene_index: &G, snp_index: &S) -> Result<Self, String>
    where
        P: AsRef<std::path::Path>,
        G: FeatureIndex,
        S: FeatureIndex,
    {
        let base = base.as_ref();

        let mut data = HashMap::new();
        data.insert(Self::EXONIC.to_string(), Scdata::read_matrix_market(base.join(Self::EXONIC), gene_index)
            .map_err(|e| format!("failed to read exonic truth from {:?}: {e}", base.join(Self::EXONIC)))?);
        data.insert(Self::INTRONIC.to_string(), Scdata::read_matrix_market(base.join(Self::INTRONIC), gene_index)
            .map_err(|e| format!("failed to read intronic truth from {:?}: {e}", base.join(Self::INTRONIC)))?);
        data.insert(Self::SNP_REF.to_string(), Scdata::read_matrix_market(base.join(Self::SNP_REF), snp_index)
            .map_err(|e| format!("failed to read SNP ref truth from {:?}: {e}", base.join(Self::SNP_REF)))?);
        data.insert(Self::SNP_ALT.to_string(), Scdata::read_matrix_market(base.join(Self::SNP_ALT), snp_index)
            .map_err(|e| format!("failed to read SNP alt truth from {:?}: {e}", base.join(Self::SNP_ALT)))?);
        Ok(Self { data })
    }

    pub fn merge(&mut self, other: &Self) -> MappingInfo {
        let mut report = MappingInfo::new(None, 0.0, usize::MAX);
        for (name, incoming) in &other.data {
            let local = self.data.get_mut(name).unwrap_or_else(|| {
                panic!("cannot merge QuantData: incoming dataset '{name}' does not exist locally")
            });
            report.merge(&local.merge(incoming));
        }
        report
    }


}

impl QuantData {
    /// Compare two quantification bundles.
    ///
    /// Returns `Ok(())` if all matrices match, otherwise returns the first
    /// useful discrepancy message.
    pub fn compare(&self, other: &Self) -> Result<(), String> {
        for (name, local) in &self.data {
            let Some(incoming) = other.data.get(name) else {
                return Err(format!("{name}: dataset missing from other QuantData"));
            };
            local.compare(incoming, name)?;
        }
        Ok(())
    }

    /// Return true if all quantification matrices match.
    pub fn equals(&self, other: &Self) -> bool {
        self.compare(other).is_ok()
    }
}
