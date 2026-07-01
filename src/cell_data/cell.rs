use std::collections::{HashMap, HashSet};

use crate::ambient_rna_detect::AmbientRnaDetect;
use crate::cell_data::GeneUmiHash;

use core::fmt;
use int_to_str::int_to_str::IntToStr;

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
pub struct CellData {
    pub name: u64,

    /// Registry of already seen feature_id + umi pairs.
    pub seen: HashSet<GeneUmiHash>,

    /// feature_id -> accumulated value / count
    pub total_reads: HashMap<u64, f32>,

    /// Multimapped sequences.
    pub multimapper: HashMap<u64, MultiMapper>,
}

impl fmt::Display for CellData {
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
        for key in &self.seen {
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

impl CellData {
    /// Create an empty cell container.
    pub fn new(name: u64) -> Self {
        Self {
            name,
            ..Default::default()
        }
    }

    /// Split this cell into ambient and non-ambient fractions.
    ///
    /// The decision is made per seen feature/UMI observation.
    pub fn split_ambient(&self, ambient: &AmbientRnaDetect) -> (Self, Self) {
        let mut non_ambient = Self::new(self.name);
        let mut ambient_data = Self::new(self.name);

        for gh in &self.seen {
            if ambient.is_ambient(gh) {
                ambient_data.add(*gh);
            } else {
                non_ambient.add(*gh);
            }
        }

        (non_ambient, ambient_data)
    }

    /// Merge another cell into this one.
    ///
    /// Only previously unseen `(feature_id, umi)` pairs are inserted.
    pub fn merge(&mut self, other: &CellData) {
        for gh in &other.seen {
            if !self.seen.insert(*gh) {
                #[cfg(debug_assertions)]
                eprintln!(
                    "Warning: duplicate GeneUmiHash {:?} found during merge for cell {}",
                    gh, self.name
                );
                continue;
            }

            *self.total_reads.entry(gh.0).or_insert(0.0) += 1.0;
        }
    }

    /// Insert one unique feature/UMI observation.
    ///
    /// Returns `false` if this exact feature/UMI pair was already seen.
    pub fn add(&mut self, fh: GeneUmiHash) -> bool {
        if !self.seen.insert(fh) {
            #[cfg(debug_assertions)]
            eprintln!(
                "Warning: duplicate GeneUmiHash {:?} detected in cell {} — ignoring.",
                fh, self.name
            );
            return false;
        }

        *self.total_reads.entry(fh.0).or_insert(0.0) += 1.0;
        true
    }

    /// Insert one unique feature/UMI observation with a numeric value.
    ///
    /// Duplicate feature/UMI pairs are rejected.
    pub fn add_value(&mut self, fh: GeneUmiHash, value: f32) -> bool {
        if !self.seen.insert(fh) {
            #[cfg(debug_assertions)]
            eprintln!(
                "Warning: duplicate GeneUmiHash {:?} detected in cell {} — ignoring value insert.",
                fh, self.name
            );
            return false;
        }

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
        let n = self.seen.iter().filter(|gh| gh.0 == *feature_id).count();

        if n > 0 { Some(total / n as f32) } else { None }
    }

    /// Dense row export helper.
    pub fn to_str_for_feature_ids(&self, feature_ids: &[u64], _length: usize) -> String {
        let mut data = Vec::<String>::with_capacity(feature_ids.len() + 4);

        data.push(IntToStr::u8_array_to_str(&self.name.to_le_bytes()).to_string());

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

    #[cfg(feature = "compare")]
    /// quick compare function
    pub fn compare(&self, other: &Self, label: &str, cell_id: u64) -> Result<(), String> {
        const EPS: f32 = 1e-6;

        if self.name != other.name {
            return Err(format!(
                "{label}: cell {cell_id} name differs: left={} right={}",
                self.name, other.name
            ));
        }

        if self.seen != other.seen {
            let only_left: Vec<_> = self
                .seen
                .difference(&other.seen)
                .take(10)
                .cloned()
                .collect();

            let only_right: Vec<_> = other
                .seen
                .difference(&self.seen)
                .take(10)
                .cloned()
                .collect();

            return Err(format!(
                "{label}: cell {cell_id} distinct molecule count differs: left={} right={}
  sample only-in-left (up to 10): {:?}
  sample only-in-right (up to 10): {:?}",
                self.seen.len(),
                other.seen.len(),
                only_left,
                only_right
            ));
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
        assert!(cell.seen.contains(&fh));
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

        assert!(cell.seen.contains(&fh));
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
        assert!(a.seen.contains(&fh1));
        assert!(a.seen.contains(&fh2));
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
