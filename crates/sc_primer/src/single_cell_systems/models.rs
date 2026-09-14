use super::{CellIdGenerator, RhapsodyWhitelist, TenxWhitelist};
use crate::whitelist_hash::RuntimeWhitelistHash;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SingleCellSystem {
    Rhapsody(RhapsodyWhitelist),
    Tenx(TenxWhitelist),
    Whitelist(RuntimeWhitelistHash),
}

impl CellIdGenerator for SingleCellSystem {
    fn cell_seq_for_index(&self, allocation_index: u64) -> Option<Vec<u8>> {
        match self {
            Self::Rhapsody(x) => x.cell_seq_for_index(allocation_index),
            Self::Tenx(x) => x.cell_seq_for_index(allocation_index),
            Self::Whitelist(x) => x.sequence(allocation_index as usize),
        }
    }

    fn cell_index_for_seq(&self, cell_seq: &[u8]) -> Option<u64> {
        match self {
            Self::Rhapsody(x) => x.cell_index_for_seq(cell_seq),
            Self::Tenx(x) => x.cell_index_for_seq(cell_seq),
            Self::Whitelist(x) => x.best_match(cell_seq).map(|(index, _)| index as u64),
        }
    }
}

pub type Range = (usize, usize);
