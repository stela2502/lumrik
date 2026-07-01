use super::{CellIdGenerator, RhapsodyWhitelist, TenxWhitelist};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SingleCellSystem {
    Rhapsody(RhapsodyWhitelist),
    Tenx(TenxWhitelist),
}

impl CellIdGenerator for SingleCellSystem {
    fn cell_seq_for_index(&self, allocation_index: u64) -> Option<Vec<u8>> {
        match self {
            Self::Rhapsody(x) => x.cell_seq_for_index(allocation_index),
            Self::Tenx(x) => x.cell_seq_for_index(allocation_index),
        }
    }

    fn cell_index_for_seq(&self, cell_seq: &[u8]) -> Option<u64> {
        match self {
            Self::Rhapsody(x) => x.cell_index_for_seq(cell_seq),
            Self::Tenx(x) => x.cell_index_for_seq(cell_seq),
        }
    }
}

pub type Range = (usize, usize);
