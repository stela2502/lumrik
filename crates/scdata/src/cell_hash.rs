use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::ops::{Deref, DerefMut};

use rayon::prelude::*;
use run_stats::RunStats;

/// Sparse cell-indexed storage split into 1024 buckets.
///
/// Cell ids are expected to come from `IntToDna::into_u64()`. That encoding
/// stores the beginning of the packed 2-bit sequence in the least-significant
/// bits of the `u64`, so the low 10 bits (the first five barcode bases) select
/// the bucket regardless of whether the full cell barcode is 16, 27, or 32 bases.
const CELL_HASH_BITS: u32 = 10;
const CELL_HASH_BUCKETS: usize = 1 << CELL_HASH_BITS;
const CELL_HASH_MASK: u64 = (CELL_HASH_BUCKETS as u64) - 1;

#[derive(Debug)]
pub struct CellHash<T> {
    data: [HashMap<u64, T>; CELL_HASH_BUCKETS],
}

impl<T> Default for CellHash<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> CellHash<T> {
    pub fn new() -> Self {
        Self {
            data: std::array::from_fn(|_| HashMap::new()),
        }
    }

    /// Bucket selected by the first five bases encoded by `IntToDna`.
    #[inline]
    pub fn bucket_index(cell_id: u64) -> usize {
        (cell_id & CELL_HASH_MASK) as usize
    }

    #[inline]
    pub fn get_cell(&self, cell_id: &u64) -> Option<&T> {
        self.data[Self::bucket_index(*cell_id)].get(cell_id)
    }

    #[inline]
    pub fn get_cell_mut(&mut self, cell_id: &u64) -> Option<&mut T> {
        self.data[Self::bucket_index(*cell_id)].get_mut(cell_id)
    }

    #[inline]
    pub fn entry_cell(&mut self, cell_id: u64) -> Entry<'_, u64, T> {
        self.data[Self::bucket_index(cell_id)].entry(cell_id)
    }

    pub fn cell_count(&self) -> usize {
        self.data.iter().map(HashMap::len).sum()
    }

    pub fn cells_are_empty(&self) -> bool {
        self.data.iter().all(HashMap::is_empty)
    }

    pub fn buckets(&self) -> &[HashMap<u64, T>; CELL_HASH_BUCKETS] {
        &self.data
    }

    pub fn buckets_mut(&mut self) -> &mut [HashMap<u64, T>; CELL_HASH_BUCKETS] {
        &mut self.data
    }

    /// Consume this CellHash and transform each bucket independently in Rayon.
    ///
    /// CellHash owns the bucket geometry, so callers can build a completely fresh
    /// cell tree without touching any persistent destination state. This is useful
    /// for batch-local model construction: expensive interpretation happens against
    /// bounded batch state first, and only the resulting compact tree is reconciled
    /// with historical state afterwards.
    pub fn par_map_owned<U, F>(self, map: F) -> CellHash<U>
    where
        T: Send,
        U: Send,
        F: Fn(HashMap<u64, T>) -> HashMap<u64, U> + Sync + Send,
    {
        let buckets: Vec<HashMap<u64, U>> = self.data.into_par_iter().map(map).collect();
        let data: [HashMap<u64, U>; CELL_HASH_BUCKETS] = buckets
            .try_into()
            .unwrap_or_else(|_| unreachable!("CellHash bucket count is invariant"));
        CellHash { data }
    }

    /// Consume this CellHash while reading the corresponding buckets of a second
    /// CellHash. Useful when batch-local work needs a compact read-only persistent
    /// gate with the same cell partitioning.
    pub fn par_map_owned_with<U, V, F>(self, other: &CellHash<U>, map: F) -> CellHash<V>
    where
        T: Send,
        U: Sync,
        V: Send,
        F: Fn(HashMap<u64, T>, &HashMap<u64, U>) -> HashMap<u64, V> + Sync + Send,
    {
        let buckets: Vec<HashMap<u64, V>> = self
            .data
            .into_par_iter()
            .zip(other.data.par_iter())
            .map(|(self_bucket, other_bucket)| map(self_bucket, other_bucket))
            .collect();
        let data: [HashMap<u64, V>; CELL_HASH_BUCKETS] = buckets
            .try_into()
            .unwrap_or_else(|_| unreachable!("CellHash bucket count is invariant"));
        CellHash { data }
    }

    /// Consume this CellHash while reading corresponding buckets from another
    /// CellHash, and reduce one thread-local RunStats per Rayon bucket.
    pub fn par_map_owned_with_stats<U, V, F>(
        self,
        other: &CellHash<U>,
        map: F,
    ) -> (CellHash<V>, RunStats)
    where
        T: Send,
        U: Sync,
        V: Send,
        F: Fn(HashMap<u64, T>, &HashMap<u64, U>, &mut RunStats) -> HashMap<u64, V> + Sync + Send,
    {
        let bucket_results: Vec<(HashMap<u64, V>, RunStats)> = self
            .data
            .into_par_iter()
            .zip(other.data.par_iter())
            .map(|(self_bucket, other_bucket)| {
                let mut local_stats = RunStats::new(None, 0.0, 0);
                let mapped = map(self_bucket, other_bucket, &mut local_stats);
                (mapped, local_stats)
            })
            .collect();

        let mut stats = RunStats::new(None, 0.0, 0);
        let buckets: Vec<HashMap<u64, V>> = bucket_results
            .into_iter()
            .map(|(bucket, local_stats)| {
                stats.merge(&local_stats);
                bucket
            })
            .collect();
        let data: [HashMap<u64, V>; CELL_HASH_BUCKETS] = buckets
            .try_into()
            .unwrap_or_else(|_| unreachable!("CellHash bucket count is invariant"));
        (CellHash { data }, stats)
    }

    /// Merge corresponding cell buckets in parallel while collecting mapping telemetry.
    ///
    /// `CellHash` owns the partitioning invariant: a cell id can only occur in one
    /// corresponding bucket in both stores. Callers therefore only define the
    /// per-cell merge semantics; the bucket-level Rayon fan-out and RunStats
    /// reduction live here once for every CellHash consumer.
    pub fn par_merge_with_mapping<U, F>(&mut self, other: &CellHash<U>, merge: F) -> RunStats
    where
        T: Send,
        U: Sync,
        F: Fn(&mut HashMap<u64, T>, &HashMap<u64, U>, &mut RunStats) + Sync + Send,
    {
        self.data
            .par_iter_mut()
            .zip(other.data.par_iter())
            .map(|(self_bucket, other_bucket)| {
                let mut local_mapping = RunStats::new(None, 0.0, 0);
                merge(self_bucket, other_bucket, &mut local_mapping);
                local_mapping
            })
            .reduce(
                || RunStats::new(None, 0.0, 0),
                |mut mapping, local_mapping| {
                    mapping.merge(&local_mapping);
                    mapping
                },
            )
    }

    /// Consuming variant of `par_merge_with_mapping`.
    ///
    /// This is useful for compact batch deltas whose contents should be moved into
    /// persistent storage rather than cloned. Each source bucket is consumed by the
    /// Rayon worker that owns the corresponding mutable destination bucket.
    pub fn par_merge_owned_with_mapping<U, F>(&mut self, other: CellHash<U>, merge: F) -> RunStats
    where
        T: Send,
        U: Send,
        F: Fn(&mut HashMap<u64, T>, HashMap<u64, U>, &mut RunStats) + Sync + Send,
    {
        self.data
            .par_iter_mut()
            .zip(other.data.into_par_iter())
            .map(|(self_bucket, other_bucket)| {
                let mut local_mapping = RunStats::new(None, 0.0, 0);
                merge(self_bucket, other_bucket, &mut local_mapping);
                local_mapping
            })
            .reduce(
                || RunStats::new(None, 0.0, 0),
                |mut mapping, local_mapping| {
                    mapping.merge(&local_mapping);
                    mapping
                },
            )
    }
}

impl<T> Deref for CellHash<T> {
    type Target = [HashMap<u64, T>; CELL_HASH_BUCKETS];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T> DerefMut for CellHash<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl<T> IntoIterator for CellHash<T> {
    type Item = HashMap<u64, T>;
    type IntoIter = std::array::IntoIter<HashMap<u64, T>, CELL_HASH_BUCKETS>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a CellHash<T> {
    type Item = &'a HashMap<u64, T>;
    type IntoIter = std::slice::Iter<'a, HashMap<u64, T>>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut CellHash<T> {
    type Item = &'a mut HashMap<u64, T>;
    type IntoIter = std::slice::IterMut<'a, HashMap<u64, T>>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter_mut()
    }
}
