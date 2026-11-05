use std::collections::{hash_map::Entry, HashMap, HashSet};
use crate::umi_evidence::UmiEvidence;

use crate::ambient_rna_detect::AmbientRnaDetect;
use crate::cell_data::GeneUmiHash;

use core::fmt;
use int_to_dna::int_to_dna::IntToDna;
use run_stats::RunStats;

/// Struct describing multimapping reads.
#[derive(Clone, Debug, Default)]
pub struct MultiMapper {
    pub features: HashSet<u64>,
    pub positions: Vec<usize>,
}

/// Cell-local feature store in a global feature space.
///
/// Logic:
/// - `seen` is the registry of already observed `(feature_id, umi)` pairs
/// - `total_reads` stores the accumulated per-feature value
/// - duplicate `(feature_id, umi)` pairs are rejected
#[derive(Clone, Debug, Default)]
pub struct CellData<T: UmiEvidence = ()> {
    pub name: u64,

    /// Registry of already seen feature_id + umi pairs.
    pub seen: HashMap<GeneUmiHash, T>,

    /// feature_id -> accumulated value / count
    pub total_reads: HashMap<u64, f32>,

    /// Multimapped sequences.
    pub multimapper: HashMap<u64, MultiMapper>,
}

impl<'a, T: UmiEvidence> IntoIterator for &'a CellData<T> {
    type Item = (&'a u64, &'a f32);
    type IntoIter = std::collections::hash_map::Iter<'a, u64, f32>;

    fn into_iter(self) -> Self::IntoIter {
        self.total_reads.iter()
    }
}

impl<T: UmiEvidence> fmt::Display for CellData<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "CellData Summary")?;
        writeln!(f, "----------------")?;
        writeln!(f, "Cell Name (ID): {}", self.name)?;
        writeln!(f, "Total UMIs: {}", self.total_umis())?;
        writeln!(f, "Seen Feature/UMI Pairs: {}", self.seen.len())?;

        writeln!(f, "\nPer-Feature Totals:")?;
        for (feature_id, count) in &self.total_reads {
            writeln!(f, "  {} => {}", feature_id, count)?;
        }

        writeln!(f, "\nSeen Feature/UMI Registry:")?;
        for key in self.seen.keys() {
            writeln!(f, "  {:?}", key)?;
        }

        writeln!(f, "\nMultimappers ({} entries):", self.multimapper.len())?;
        for (seq_id, mm) in &self.multimapper {
            writeln!(
                f,
                "  Seq {}: Features {:?}, Positions {:?}",
                seq_id, mm.features, mm.positions
            )?;
        }

        Ok(())
    }
}

impl CellData<()> {
    /// Create an empty cell container.
    pub fn new(name: u64) -> Self {
        Self {
            name,
            ..Default::default()
        }
    }

}

impl<T: UmiEvidence> CellData<T> {
    pub fn new_with_evidence(name: u64) -> Self {
        Self { name, seen: HashMap::new(), total_reads: HashMap::new(), multimapper: HashMap::new() }
    }

    /// Split this cell into ambient and non-ambient fractions.
    ///
    /// The decision is made per seen feature/UMI observation.
    pub fn split_ambient(&self, ambient: &AmbientRnaDetect) -> (Self, Self) {
        let mut non_ambient = Self::new_with_evidence(self.name);
        let mut ambient_data = Self::new_with_evidence(self.name);

        for (gh, evidence) in &self.seen {
            if ambient.is_ambient(gh) {
                ambient_data.add_evidence(*gh, evidence.clone());
            } else {
                non_ambient.add_evidence(*gh, evidence.clone());
            }
        }

        (non_ambient, ambient_data)
    }

    /// Merge another cell into this one.
    ///
    /// Only previously unseen `(feature_id, umi)` pairs are inserted.
    pub fn merge(&mut self, other: &CellData<T>) {
        let mut mapping = RunStats::new(None, 0.0, 0);
        self.merge_with_mapping(other, &mut mapping);
    }

    /// Merge another cell while accounting accepted and duplicate UMIs.
    ///
    /// The decision is made exactly where the `(feature_id, umi)` pair is
    /// inserted into `seen`: successful inserts are accepted, failed inserts
    /// are UMI-level duplicates.
    pub fn merge_with_mapping(&mut self, other: &CellData<T>, mapping: &mut RunStats) {
        for (gh, evidence) in &other.seen {
            match self.seen.entry(*gh) {
                Entry::Occupied(mut entry) => {
                    entry.get_mut().merge(evidence);
                    mapping.pcr_duplicates += 1;
                }
                Entry::Vacant(entry) => {
                    entry.insert(evidence.clone());
                    mapping.ok_reads += 1;
                    *self.total_reads.entry(gh.0).or_insert(0.0) += 1.0;
                }
            }
        }
    }

    /// Merge another cell by adding its already-accumulated numeric feature values.
    ///
    /// This is intended for disjoint real-valued evidence layers (for example
    /// unique and fractional multimapper TE counts). It preserves f32 totals
    /// instead of reconstructing each observation as 1.0.
    pub fn merge_values(&mut self, other: CellData<T>) {
        for (key, evidence) in other.seen {
            match self.seen.entry(key) {
                Entry::Occupied(mut entry) => entry.get_mut().merge(&evidence),
                Entry::Vacant(entry) => { entry.insert(evidence); }
            }
        }
        for (feature_id, value) in other.total_reads {
            *self.total_reads.entry(feature_id).or_insert(0.0) += value;
        }
        self.multimapper.extend(other.multimapper);
    }

    /// Insert one unique feature/UMI observation.
    ///
    /// Returns `false` if this exact feature/UMI pair was already seen.
    pub fn add(&mut self, fh: GeneUmiHash) -> bool {
        self.add_evidence(fh, T::default())
    }

    /// Insert or merge evidence for a feature/UMI without counting duplicates twice.
    pub fn add_evidence(&mut self, fh: GeneUmiHash, evidence: T) -> bool {
        match self.seen.entry(fh) {
            Entry::Occupied(mut entry) => {
                entry.get_mut().merge(&evidence);
                false
            }
            Entry::Vacant(entry) => {
                entry.insert(evidence);
                *self.total_reads.entry(fh.0).or_insert(0.0) += 1.0;
                true
            }
        }
    }

    pub fn add_value(&mut self, fh: GeneUmiHash, value: f32) -> bool {
        if self.seen.contains_key(&fh) { return false; }
        self.seen.insert(fh, T::default());
        *self.total_reads.entry(fh.0).or_insert(0.0) += value;
        true
    }

    /// better API - Number of unique feature/UMI observations in this cell.
    pub fn total_umis(&self) -> usize {
        self.seen.len()
    }

    /// Accumulated value for one feature id.
    pub fn total_umis_4_gene_id(&self, feature_id: &u64) -> f32 {
        *self.total_reads.get(feature_id).unwrap_or(&0.0)
    }

    /// Mean value helper.
    ///
    /// In the current model this is just the total feature value divided
    /// by the number of unique UMIs seen for that feature.
    pub fn get_mean_for_gene(&self, feature_id: &u64) -> Option<f32> {
        let total = *self.total_reads.get(feature_id)?;
        let n = self.seen.keys().filter(|gh| gh.0 == *feature_id).count();

        if n > 0 { Some(total / n as f32) } else { None }
    }

    /// Dense row export helper.
    pub fn to_str_for_feature_ids(&self, feature_ids: &[u64], _length: usize) -> String {
        let mut data = Vec::<String>::with_capacity(feature_ids.len() + 4);

        data.push(IntToDna::u8_array_to_str(&self.name.to_le_bytes()).to_string());

        let mut total = 0.0f32;
        let mut max = 0.0f32;
        let mut max_id: Option<u64> = None;

        for id in feature_ids {
            let n = *self.total_reads.get(id).unwrap_or(&0.0);

            if n > max {
                max = n;
                max_id = Some(*id);
            }

            data.push(n.to_string());
            total += n;
        }

        let mut dist2max = f32::MAX;
        for id in feature_ids {
            let n = *self.total_reads.get(id).unwrap_or(&0.0);
            let diff = max - n;

            if diff > 0.0 && diff < dist2max {
                dist2max = diff;
            }
        }

        data.push(
            max_id
                .map(|id| id.to_string())
                .unwrap_or_else(|| "na".to_string()),
        );

        data.push(if total > 0.0 {
            (max / total).to_string()
        } else {
            "0".to_string()
        });

        data.push(total.to_string());

        data.push(if max > 0.0 && dist2max.is_finite() {
            (dist2max / max).to_string()
        } else {
            "0".to_string()
        });

        data.join("\t")
    }

    /// quick compare function
    pub fn compare(&self, other: &Self, label: &str, cell_id: u64) -> Result<(), String> {
        const EPS: f32 = 1e-6;

        if self.name != other.name {
            return Err(format!(
                "{label}: cell {cell_id} name differs: left={} right={}",
                self.name, other.name
            ));
        }

        if self.seen.len() != other.seen.len()
            || self.seen.keys().any(|key| !other.seen.contains_key(key)) {
            return Err(format!("{label}: cell {cell_id} distinct UMI keys differ"));
        }
        if self.total_reads.len() != other.total_reads.len() {
            let only_left: Vec<_> = self
                .total_reads
                .keys()
                .filter(|k| !other.total_reads.contains_key(*k))
                .take(10)
                .copied()
                .collect();

            let only_right: Vec<_> = other
                .total_reads
                .keys()
                .filter(|k| !self.total_reads.contains_key(*k))
                .take(10)
                .copied()
                .collect();

            return Err(format!(
                "{label}: cell {cell_id} total_reads feature count differs: left={} right={}
  sample only-in-left: {:?}
  sample only-in-right: {:?}",
                self.total_reads.len(),
                other.total_reads.len(),
                only_left,
                only_right
            ));
        }
        for (feature_id, left_value) in &self.total_reads {
            let Some(right_value) = other.total_reads.get(feature_id) else {
                return Err(format!(
                    "{label}: cell {cell_id} feature {feature_id} exists on left but not right"
                ));
            };

            if (*left_value - *right_value).abs() > EPS {
                return Err(format!(
                    "{label}: cell {cell_id} feature {feature_id} differs: left={left_value} right={right_value}"
                ));
            }
        }

        for feature_id in other.total_reads.keys() {
            if !self.total_reads.contains_key(feature_id) {
                return Err(format!(
                    "{label}: cell {cell_id} feature {feature_id} exists on right but not left"
                ));
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell_data::GeneUmiHash;

    #[test]
    fn test_new_cell() {
        let cell = CellData::new(42);

        assert_eq!(cell.name, 42);
        assert_eq!(cell.total_umis(), 0);
        assert!(cell.seen.is_empty());
        assert!(cell.total_reads.is_empty());
    }

    #[test]
    fn test_add_unique_feature_umi() {
        let mut cell = CellData::new(1);

        let fh = GeneUmiHash(10, 123);
        let inserted = cell.add(fh);

        assert!(inserted);
        assert_eq!(cell.total_umis(), 1);
        assert!(cell.seen.contains_key(&fh));
        assert_eq!(cell.total_reads.get(&10), Some(&1.0));
    }

    #[test]
    fn test_add_duplicate_rejected() {
        let mut cell = CellData::new(1);

        let fh = GeneUmiHash(10, 123);

        assert!(cell.add(fh));
        assert!(!cell.add(fh));

        assert_eq!(cell.total_umis(), 1);
        assert_eq!(cell.total_reads.get(&10), Some(&1.0));
    }

    #[test]
    fn test_add_value_new_only() {
        let mut cell = CellData::new(1);

        let fh = GeneUmiHash(5, 111);

        assert!(cell.add_value(fh, 2.0));
        assert!(!cell.add_value(fh, 3.0));

        assert!(cell.seen.contains_key(&fh));
        assert_eq!(cell.total_reads.get(&5), Some(&2.0));
        assert_eq!(cell.total_umis(), 1);
    }

    #[test]
    fn test_total_umis() {
        let mut cell = CellData::new(1);

        cell.add(GeneUmiHash(1, 10));
        cell.add(GeneUmiHash(1, 11));
        cell.add(GeneUmiHash(2, 12));

        assert_eq!(cell.total_umis(), 3);
    }

    #[test]
    fn test_merge_cells() {
        let mut a = CellData::new(1);
        let mut b = CellData::new(1);

        let fh1 = GeneUmiHash(1, 10);
        let fh2 = GeneUmiHash(2, 20);

        a.add(fh1);
        b.add(fh2);

        a.merge(&b);

        assert_eq!(a.total_umis(), 2);
        assert!(a.seen.contains_key(&fh1));
        assert!(a.seen.contains_key(&fh2));
        assert_eq!(a.total_reads.get(&1), Some(&1.0));
        assert_eq!(a.total_reads.get(&2), Some(&1.0));
    }

    #[test]
    fn test_total_umis_for_feature() {
        let mut cell = CellData::new(1);

        cell.add(GeneUmiHash(5, 1));
        cell.add(GeneUmiHash(5, 2));
        cell.add(GeneUmiHash(3, 1));

        assert_eq!(cell.total_umis_4_gene_id(&5), 2.0);
        assert_eq!(cell.total_umis_4_gene_id(&3), 1.0);
        assert_eq!(cell.total_umis_4_gene_id(&9), 0.0);
    }

    #[test]
    fn test_mean_for_feature() {
        let mut cell = CellData::new(1);

        cell.add_value(GeneUmiHash(1, 1), 2.0);
        cell.add_value(GeneUmiHash(1, 2), 4.0);

        let mean = cell.get_mean_for_gene(&1).unwrap();

        assert!((mean - 3.0).abs() < 1e-6);
    }
}

#[cfg(test)]
mod evidence_tests {
    use super::*;
    use crate::umi_evidence::SnpUmiEvidence;

    #[test]
    fn duplicate_umi_merges_evidence_without_increasing_count() {
        let mut cell = CellData::<SnpUmiEvidence>::new_with_evidence(42);
        let key = GeneUmiHash(7, 99);
        let mut ref_support = SnpUmiEvidence::default();
        ref_support.add_ref(17);
        let mut alt_support = SnpUmiEvidence::default();
        alt_support.add_alt(17);
        assert!(cell.add_evidence(key, ref_support));
        assert!(!cell.add_evidence(key, alt_support));
        assert_eq!(cell.total_umis_4_gene_id(&7), 1.0);
        assert_eq!(cell.seen[&key].conflicts().count(), 1);
    }

    #[test]
    fn merge_across_cells_merges_evidence() {
        let key = GeneUmiHash(7, 99);
        let mut left = CellData::<SnpUmiEvidence>::new_with_evidence(42);
        let mut right = CellData::<SnpUmiEvidence>::new_with_evidence(42);
        let mut a = SnpUmiEvidence::default(); a.add_ref(17);
        let mut b = SnpUmiEvidence::default(); b.add_alt(17);
        left.add_evidence(key, a);
        right.add_evidence(key, b);
        left.merge(&right);
        assert_eq!(left.total_umis(), 1);
        assert_eq!(left.seen[&key].conflicts().count(), 1);
    }
}
