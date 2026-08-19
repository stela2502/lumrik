use std::collections::HashMap;

use scdata::Scdata;

#[derive(Debug, Clone, Copy)]
pub struct GuideObservation {
    pub cell_id: u64,
    pub guide_id: u32,
    pub count: u32,
}

/// Two views of the same sparse guide data:
///
/// * `cells` is the cell-major sparse representation owned by `scdata`.
/// * `by_guide` is the complementary guide-major view used by model fitting.
///
/// The guide-major view deliberately stores only non-zero entries.
pub struct GuideDataset {
    pub cells: Scdata,
    pub by_guide: Vec<Vec<GuideObservation>>,
    pub cell_ids: Vec<u64>,
    pub barcode_by_id: HashMap<u64, String>,
}

impl GuideDataset {
    pub fn n_cells(&self) -> usize {
        self.cell_ids.len()
    }

    pub fn n_guides(&self) -> usize {
        self.by_guide.len()
    }

    pub fn cell_total(&self, cell_id: u64) -> u32 {
        self.cells
            .get(&cell_id)
            .map(|cell| {
                cell.total_reads
                    .values()
                    .filter(|v| **v > 0.0)
                    .map(|v| *v as u32)
                    .sum()
            })
            .unwrap_or(0)
    }

    pub fn observations_for_cell(&self, cell_id: u64) -> Vec<GuideObservation> {
        let Some(cell) = self.cells.get(&cell_id) else {
            return Vec::new();
        };

        let mut out: Vec<_> = cell
            .total_reads
            .iter()
            .filter_map(|(&guide_id, &value)| {
                if value <= 0.0 {
                    None
                } else {
                    Some(GuideObservation {
                        cell_id,
                        guide_id: guide_id as u32,
                        count: value.round() as u32,
                    })
                }
            })
            .collect();

        out.sort_unstable_by_key(|x| x.guide_id);
        out
    }
}
