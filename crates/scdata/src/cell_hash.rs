use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::ops::{Deref, DerefMut};

/// Sparse cell-indexed storage split into 256 buckets.
///
/// Cell ids are expected to come from `IntToStr::into_u64()`. That encoding
/// stores the first packed 4-base byte (`u8_encoded[0]`) in the least
/// significant byte of the `u64`, so those first four barcode bases select the
/// bucket regardless of whether the full cell barcode is 16, 27, or 32 bases.
#[derive(Debug)]
pub struct CellHash<T> {
    data: [HashMap<u64, T>; 256],
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

    /// Bucket selected by the first four bases encoded by `IntToStr`.
    #[inline]
    pub fn bucket_index(cell_id: u64) -> usize {
        cell_id as u8 as usize
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

    pub fn buckets(&self) -> &[HashMap<u64, T>; 256] {
        &self.data
    }

    pub fn buckets_mut(&mut self) -> &mut [HashMap<u64, T>; 256] {
        &mut self.data
    }
}

impl<T> Deref for CellHash<T> {
    type Target = [HashMap<u64, T>; 256];

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
    type IntoIter = std::array::IntoIter<HashMap<u64, T>, 256>;

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

#[cfg(test)]
mod tests {
    use super::CellHash;
    use int_to_str::IntToStr;

    fn assert_first_encoded_byte_selects_bucket(barcode: &[u8]) {
        let encoded = IntToStr::new(barcode);
        let cell_id = encoded.into_u64();
        assert_eq!(
            CellHash::<()>::bucket_index(cell_id),
            encoded.u8_encoded[0] as usize
        );
    }

    #[test]
    fn bucket_is_first_int_to_str_byte_for_16bp_cell() {
        assert_first_encoded_byte_selects_bucket(b"GCTCGTAGCTAGCTAG");
    }

    #[test]
    fn bucket_is_first_int_to_str_byte_for_27bp_cell() {
        assert_first_encoded_byte_selects_bucket(b"GCTCGTAGCTAGCTAGCTAGCTAGCTA");
    }

    #[test]
    fn first_four_bases_control_bucket_independent_of_barcode_length() {
        let short = IntToStr::new(b"GCTCGTAGCTAGCTAG");
        let long = IntToStr::new(b"GCTCGTAGCTAGCTAGCTAGCTAGCTA");
        let other = IntToStr::new(b"ACTCGTAGCTAGCTAG");

        assert_eq!(
            CellHash::<()>::bucket_index(short.into_u64()),
            CellHash::<()>::bucket_index(long.into_u64())
        );
        assert_ne!(
            CellHash::<()>::bucket_index(short.into_u64()),
            CellHash::<()>::bucket_index(other.into_u64())
        );
    }
}
