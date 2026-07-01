use std::collections::{HashMap, HashSet};

use core::fmt;
use rayon::prelude::*;

use mapping_info::MappingInfo;

use crate::cell_data::CellData;
use crate::cell_data::GeneUmiHash;
use crate::{FeatureIndex, MatrixValueType};

/// Sparse single-cell count store.
///
/// Cells are partitioned into 256 buckets by the top byte of the cell id.
/// Each bucket stores `CellData` objects keyed by the full cell id.
pub struct Scdata {
    /// 256 buckets, indexed by top 8 bits of the cell id.
    pub(crate) data: [HashMap<u64, CellData>; u8::MAX as usize + 1],

    /// Cached ordered feature ids currently observed in retained cells.
    pub(crate) feature_ids_with_data: Vec<u64>,

    /// Cached number of stored feature entries across all retained cells.
    pub(crate) total_feature_data_entries: usize,

    /// Ordered list of cells that define export column order.
    pub(crate) export_cell_ids: Vec<u64>,

    /// True once export selection/filtering has been evaluated.
    checked: bool,

    /// Preferred thread count for internal parallel helpers.
    pub num_threads: usize,

    /// Numeric matrix value type.
    pub(crate) value_type: MatrixValueType,
}

/// Sparse single-cell UMI count store.
///
/// # Example
///
/// ```no_run
/// use std::collections::HashMap;
/// use mapping_info::MappingInfo;
/// use std::path::PathBuf;
/// use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};
///
/// struct SimpleIndex {
///     names: Vec<String>,
///     ids: HashMap<String, u64>,
/// }
///
/// impl SimpleIndex {
///     fn new(names: Vec<&str>) -> Self {
///         let names: Vec<String> = names.into_iter().map(|s| s.to_string()).collect();
///         let ids = names
///             .iter()
///             .enumerate()
///             .map(|(i, n)| (n.clone(), i as u64))
///             .collect();
///         Self { names, ids }
///     }
/// }
///
/// impl FeatureIndex for SimpleIndex {
///     fn feature_name(&self, feature_id: u64) -> &str {
///         &self.names[feature_id as usize]
///     }
///
///     fn feature_id(&self, name: &str) -> Option<u64> {
///         self.ids.get(name).copied()
///     }
///
///     fn ordered_feature_ids(&self) -> Vec<u64> {
///         (0..self.names.len() as u64).collect()
///     }
///
///     fn to_10x_feature_line(&self, feature_id: u64) -> String {
///         let name = self.feature_name(feature_id);
///         format!("{name}\t{name}\tGene Expression")
///     }
/// }
///
/// let mut data = Scdata::new(1, MatrixValueType::Integer);
/// let mut report = MappingInfo::new(None, 0.0, 0);
/// let index = SimpleIndex::new(vec!["GeneA", "GeneB"]);
///
/// let gene_a = index.feature_id("GeneA").unwrap();
///
/// data.try_insert(&1_u64, GeneUmiHash(gene_a, 100), 0.0, &mut report);
/// data.try_insert(&1_u64, GeneUmiHash(gene_a, 101), 0.0, &mut report);
///
/// data.finalize_for_export(0, &index);
/// let out = PathBuf::from("example_sparse_out");
/// let _ = data.write_sparse(&out, &index);
///
/// assert_eq!(data.passing_cells(), 1);
/// ```
impl Default for Scdata {
    /// Create an integer-valued matrix with one default worker thread.
    fn default() -> Self {
        Self::new(1, MatrixValueType::Integer)
    }
}

impl fmt::Display for Scdata {
    /// Print a short structural summary of the current matrix state.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Scdata Summary")?;
        writeln!(f, "===============")?;
        writeln!(f, "Checked: {}", self.checked)?;
        writeln!(f, "Export Cells: {}", self.export_cell_ids.len())?;
        writeln!(f, "Matrix Value Type: {:?}", self.value_type)?;
        writeln!(
            f,
            "Features with Data: [{} entries]",
            self.feature_ids_with_data.len()
        )?;

        let total_cells: usize = self.data.iter().map(|m| m.len()).sum();
        let non_empty_buckets = self.data.iter().filter(|m| !m.is_empty()).count();

        writeln!(f, "Total Cells: {}", total_cells)?;
        writeln!(f, "Non-empty Buckets: {}", non_empty_buckets)?;
        Ok(())
    }
}

impl Scdata {
    /// Create a new empty sparse matrix store.
    pub fn new(num_threads: usize, value_type: MatrixValueType) -> Self {
        let data = std::array::from_fn(|_| HashMap::<u64, CellData>::new());

        Self {
            data,
            feature_ids_with_data: Vec::new(),
            total_feature_data_entries: 0,
            export_cell_ids: Vec::new(),
            checked: false,
            num_threads,
            value_type,
        }
    }

    /// Get the current matrix value type, or update it if a new one is provided.
    pub fn value_type(&mut self, value_type: Option<MatrixValueType>) -> &MatrixValueType {
        if let Some(val) = value_type {
            self.value_type = val;
        }
        &self.value_type
    }

    /// Map a cell id to one of the 256 storage buckets.
    #[inline]
    fn to_key(&self, name: &u64) -> usize {
        (*name >> 56) as usize
    }

    /// Invalidate cached export-derived state after matrix mutation.
    fn invalidate_export_cache(&mut self) {
        self.checked = false;
        self.export_cell_ids.clear();
        self.feature_ids_with_data.clear();
        self.total_feature_data_entries = 0;
    }

    /// Iterate over all currently stored cells across all buckets.
    pub(crate) fn values(&self) -> impl Iterator<Item = &CellData> {
        self.data.iter().flat_map(|map| map.values())
    }

    /// Return all currently stored cell ids.
    pub fn keys(&self) -> Vec<u64> {
        self.data
            .iter()
            .flat_map(|map| map.keys())
            .copied()
            .collect()
    }

    /// Check whether no cells are currently stored.
    pub fn is_empty(&self) -> bool {
        self.data.iter().all(|m| m.is_empty())
    }

    /// Return the total number of stored cells.
    pub fn len(&self) -> usize {
        self.data.iter().map(|m| m.len()).sum()
    }

    /// Get read-only access to one cell by id.
    pub fn get(&self, key: &u64) -> Option<&CellData> {
        let index = self.to_key(key);
        self.data[index].get(key)
    }

    /// Merge another matrix into this one in parallel, bucket by bucket.
    pub fn merge(&mut self, other: &Scdata) {
        if other.is_empty() {
            return;
        }

        self.invalidate_export_cache();

        self.data
            .par_iter_mut()
            .enumerate()
            .for_each(|(index, self_bucket)| {
                if let Some(other_bucket) = other.data.get(index) {
                    for (cell_name, other_cell) in other_bucket {
                        match self_bucket.entry(*cell_name) {
                            std::collections::hash_map::Entry::Occupied(mut entry) => {
                                entry.get_mut().merge(other_cell);
                            }
                            std::collections::hash_map::Entry::Vacant(entry) => {
                                entry.insert(other_cell.clone());
                            }
                        }
                    }
                }
            });
    }

    /// Merge another matrix into this one by consuming the source in a single thread.
    pub fn merge_single_thread(&mut self, mut other: Scdata) {
        if other.is_empty() {
            return;
        }

        self.invalidate_export_cache();

        for other_cell in other.data.iter_mut().flat_map(|map| map.values_mut()) {
            let index = self.to_key(&other_cell.name);
            match self.data[index].entry(other_cell.name) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().merge(other_cell);
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert(std::mem::take(other_cell));
                }
            }
        }
    }

    /// Insert one unique feature/UMI observation for a cell.
    ///
    /// Returns `false` if the exact feature/UMI pair already existed.
    pub fn try_insert(
        &mut self,
        name: &u64,
        data: GeneUmiHash,
        _value: f32,
        report: &mut MappingInfo,
    ) -> bool {
        let index = self.to_key(name);
        self.invalidate_export_cache();

        let cell_info = self.data[index]
            .entry(*name)
            .or_insert_with(|| CellData::new(*name));

        report.ok_reads += 1;

        if !cell_info.add(data) {
            report.pcr_duplicates += 1;
            report.local_dup += 1;
            false
        } else {
            true
        }
    }

    /// Insert or accumulate a numeric value for one feature/UMI observation.
    pub fn try_insert_value(
        &mut self,
        name: &u64,
        data: GeneUmiHash,
        value: f32,
        report: &mut MappingInfo,
    ) -> bool {
        let index = self.to_key(name);
        self.invalidate_export_cache();

        let cell_info = self.data[index]
            .entry(*name)
            .or_insert_with(|| CellData::new(*name));

        report.ok_reads += 1;

        cell_info.add_value(data, value)
    }

    pub(crate) fn check_sparse_export_ready(&self) -> Result<(), String> {
        if !self.checked {
            return Err("Sparse export requires finalize_for_export(...) first.".to_string());
        }

        if self.export_cell_ids.is_empty() {
            eprintln!("Sparse export failed: no export cells available.",);
        }

        Ok(())
    }

    /// Return the number of cells currently selected for export.
    pub fn passing_cells(&self) -> usize {
        self.export_cell_ids.len()
    }

    /// Return all currently stored cell ids.
    pub fn cell_ids(&self) -> HashSet<u64> {
        self.data
            .iter()
            .flat_map(|bucket| bucket.keys().copied())
            .collect()
    }

    /// Return the set of cell ids with at least `min_count` UMIs.
    fn passing_cell_set_by_umi(&self, min_count: usize) -> HashSet<u64> {
        let keys = self.keys();
        if keys.is_empty() {
            return HashSet::new();
        }

        let n_threads = self.num_threads.max(1);
        let chunk_size = keys.len() / n_threads + 1;

        let passing: Vec<u64> = keys
            .par_chunks(chunk_size)
            .flat_map(|chunk| {
                let mut keep = Vec::<u64>::with_capacity(chunk.len());
                for key in chunk {
                    if let Some(cell) = self.get(key)
                        && cell.total_umis() >= min_count
                    {
                        keep.push(*key);
                    }
                }
                keep
            })
            .collect();

        passing.into_iter().collect()
    }

    /// Rebuild cached export metadata from the currently retained cells.
    ///
    /// The feature id cache is collected from the per-cell `total_reads` keys.
    /// The total number of exported sparse entries is the sum of
    /// `cell.total_reads.len()` across all retained cells.
    fn rebuild_feature_ids_with_data<I: FeatureIndex>(&mut self, index: &I) {
        let (observed_feature_ids, total_entries) = self
            .data
            .par_iter()
            .map(|bucket| {
                let mut local_ids = HashSet::with_capacity(64);
                let mut local_entries = 0usize;

                for cell in bucket.values() {
                    for feature_id in cell.total_reads.keys() {
                        local_ids.insert(*feature_id);
                    }

                    local_entries += cell.total_reads.len();
                }

                (local_ids, local_entries)
            })
            .reduce(
                || (HashSet::new(), 0usize),
                |(mut ids_a, n_a), (ids_b, n_b)| {
                    ids_a.extend(ids_b);
                    (ids_a, n_a + n_b)
                },
            );

        let feature_ids_with_data = index
            .ordered_feature_ids()
            .into_iter()
            .filter(|fid| observed_feature_ids.contains(fid))
            .collect();

        self.feature_ids_with_data = feature_ids_with_data;
        self.total_feature_data_entries = total_entries;
    }

    /// Restrict the matrix to a predefined set of cell ids.
    ///
    /// Also establishes deterministic export column order.
    fn restrict_to_cells(&mut self, keep: &HashSet<u64>) {
        for bucket in &mut self.data {
            bucket.retain(|cell_id, _| keep.contains(cell_id));
        }

        let mut export_cell_ids: Vec<u64> = self
            .data
            .iter()
            .flat_map(|bucket| bucket.keys().copied())
            .collect();

        export_cell_ids.sort_unstable();

        self.export_cell_ids = export_cell_ids;
    }

    /// Prepare the object for export.
    ///
    /// This applies the UMI cutoff, retains only passing cells,
    /// establishes deterministic export order, rebuilds export caches,
    /// and marks the object as checked.
    pub fn finalize_for_export<I: FeatureIndex>(&mut self, min_total_umis: usize, index: &I) {
        let keep = self.passing_cell_set_by_umi(min_total_umis);
        self.restrict_to_cells(&keep);
        self.rebuild_feature_ids_with_data(index);
        self.checked = true;
    }

    /// Prepare the object for export.
    ///
    /// Here a given set of cells are marked for export
    pub fn finalize_for_cells<I: FeatureIndex>(&mut self, keep: &HashSet<u64>, index: &I) {
        self.restrict_to_cells(keep);
        self.rebuild_feature_ids_with_data(index);
        self.checked = true;
    }

    /// Return the ordered export cell ids.
    pub fn export_cell_ids(&self) -> &[u64] {
        &self.export_cell_ids
    }

    pub fn observed_feature_ids(&self) -> HashSet<u64> {
        self.data
            .par_iter()
            .map(|bucket| {
                let mut local = HashSet::new();

                for cell_data in bucket.values() {
                    local.extend(cell_data.total_reads.keys().copied());
                }

                local
            })
            .reduce(HashSet::new, |mut a, b| {
                a.extend(b);
                a
            })
    }

    pub fn retain_features(&mut self, keep: &HashSet<u64>) {
        self.data.par_iter_mut().for_each(|bucket| {
            for cell_data in bucket.values_mut() {
                cell_data
                    .total_reads
                    .retain(|feature_id, _| keep.contains(feature_id));

                cell_data.seen.retain(|gene_umi| keep.contains(&gene_umi.0));

                cell_data
                    .multimapper
                    .retain(|feature_id, _| keep.contains(feature_id));
            }

            bucket.retain(|_, cell_data| !cell_data.total_reads.is_empty());
        });

        self.checked = false;
        self.feature_ids_with_data.clear();
        self.total_feature_data_entries = 0;
        self.export_cell_ids.clear();
    }

    #[cfg(feature = "compare")]
    /// quickly compare two Scdata objects.
    pub fn compare(&self, other: &Self, label: &str) -> Result<(), String> {
        if self.export_cell_ids != other.export_cell_ids {
            return Err(format!("{label}: export cells differ"));
        }

        if self.feature_ids_with_data != other.feature_ids_with_data {
            return Err(format!("{label}: feature ids differ"));
        }

        if self.total_feature_data_entries != other.total_feature_data_entries {
            return Err(format!(
                "{label}: total feature entries differ: left={} right={}",
                self.total_feature_data_entries, other.total_feature_data_entries
            ));
        }

        for bucket_id in 0..self.data.len() {
            let left = &self.data[bucket_id];
            let right = &other.data[bucket_id];

            if left.len() != right.len() {
                return Err(format!(
                    "{label}: bucket {bucket_id} size differs: left={} right={}",
                    left.len(),
                    right.len()
                ));
            }

            for (cell_id, left_cell) in left {
                let Some(right_cell) = right.get(cell_id) else {
                    return Err(format!(
                        "{label}: cell {cell_id} exists on left but not right"
                    ));
                };

                left_cell.compare(right_cell, label, *cell_id)?;
            }

            for cell_id in right.keys() {
                if !left.contains_key(cell_id) {
                    return Err(format!(
                        "{label}: cell {cell_id} exists on right but not left"
                    ));
                }
            }
        }

        Ok(())
    }
}
