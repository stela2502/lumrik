use std::collections::HashMap;
use anyhow::Result;
use int_to_str::IntToStr;
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

    pub fn from_scdata(
        cells: Scdata,
        cell_barcode_len: usize,
    ) -> Result<Self> {
        let n_features = cells.n_features();
        let mut by_guide =
            vec![Vec::new(); n_features];

        let cell_ids =
            cells.export_cell_ids();

        let barcode_by_id = cell_ids
            .iter()
            .map(|&cell_id| {
                (
                    cell_id,
                    IntToStr::from_u64(cell_id)
                        .to_string(cell_barcode_len),
                )
            })
            .collect();

        for cell_id in cell_ids {
            let Some(cell) =
                cells.get(cell_id)
            else {
                continue;
            };

            for (feature_id, value) in cell {
                let guide_id =
                    *feature_id as u32;

                if guide_id as usize >= n_features {
                    panic!(
                        "feature id {guide_id} exceeds feature count {n_features}"
                    );
                }

                let count =
                    value.round() as u32;

                if count == 0 {
                    continue;
                }

                by_guide[guide_id as usize].push(
                    GuideObservation {
                        cell_id : *cell_id,
                        guide_id,
                        count,
                    }
                );
            }
        }

        Ok(Self {
            cell_ids: cell_ids.to_vec(),
            cells,
            by_guide,
            barcode_by_id,
        })
    }
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
