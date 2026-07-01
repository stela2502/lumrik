pub trait CellIdGenerator {
    /// return the cell sequence for the cell id.
    fn cell_seq_for_index(&self, index: u64) -> Option<Vec<u8>>;

    /// Returns the allocation index for a valid exported cell sequence.
    fn cell_index_for_seq(&self, cell_seq: &[u8]) -> Option<u64>;
}
