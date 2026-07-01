/// One 8bp 2bit table hit.
///
/// This is not a feature/sample. It is one kmer position inside one feature sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TagEntry {
    /// Index into `FastTagMapper.features`.
    pub feature_index: usize,

    /// Position of this 8mer inside the original feature/tag sequence.
    pub tag_pos: usize,
}
