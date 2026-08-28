use std::collections::HashMap;

use anyhow::Result;
use int_to_str::IntToStr;
use scdata::{FeatureIndex, Scdata};

#[derive(Debug, Clone, Copy)]
pub struct GuideObservation {
    pub cell_id: u64,
    /// Dense model-local feature slot. Translate back with `feature_id()`.
    pub guide_id: u32,
    pub count: u32,
}

/// Lightweight model view over lumrik-native [`Scdata`].
///
/// `Scdata` remains the canonical cell-major store. `by_guide` is only the
/// complementary sparse view needed by the statistical model. Feature IDs are
/// not assumed to be dense: `feature_ids[guide_id]` maps the model-local slot
/// back to the real [`FeatureIndex`] id.
pub struct GuideDataset<'a> {
    pub cells: &'a Scdata,
    pub by_guide: Vec<Vec<GuideObservation>>,
    pub cell_ids: Vec<u64>,
    pub feature_ids: Vec<u64>,
    pub cell_barcode_len: usize,
}

impl<'a> GuideDataset<'a> {
    pub fn from_scdata<I: FeatureIndex>(
        cells: &'a Scdata,
        feature_index: &I,
        cell_barcode_len: usize,
    ) -> Result<Self> {
        let mut cell_ids = cells.keys();
        cell_ids.sort_unstable();
        Self::from_scdata_with_cells(cells, cell_ids, feature_index, cell_barcode_len)
    }

    /// Build a model view while preserving an explicit canonical cell set.
    /// Cells with zero feature counts remain present in `cell_ids` even though
    /// Scdata has no sparse entry for them.
    pub fn from_scdata_with_cells<I: FeatureIndex>(
        cells: &'a Scdata,
        mut cell_ids: Vec<u64>,
        feature_index: &I,
        cell_barcode_len: usize,
    ) -> Result<Self> {
        cell_ids.sort_unstable();
        cell_ids.dedup();
        let feature_ids = feature_index.ordered_feature_ids();
        let feature_to_guide: HashMap<u64, u32> = feature_ids
            .iter()
            .enumerate()
            .map(|(guide_id, feature_id)| (*feature_id, guide_id as u32))
            .collect();

        let mut by_guide = vec![Vec::new(); feature_ids.len()];

        for cell_id in &cell_ids {
            let Some(cell) = cells.get(cell_id) else {
                continue;
            };

            for (&feature_id, &value) in &cell.total_reads {
                let Some(&guide_id) = feature_to_guide.get(&feature_id) else {
                    continue;
                };

                let count = value.round() as u32;
                if count == 0 {
                    continue;
                }

                by_guide[guide_id as usize].push(GuideObservation {
                    cell_id: *cell_id,
                    guide_id,
                    count,
                });
            }
        }

        Ok(Self {
            cells,
            by_guide,
            cell_ids,
            feature_ids,
            cell_barcode_len,
        })
    }

    pub fn n_cells(&self) -> usize {
        self.cell_ids.len()
    }

    pub fn n_guides(&self) -> usize {
        self.feature_ids.len()
    }

    pub fn feature_id(&self, guide_id: u32) -> u64 {
        self.feature_ids[guide_id as usize]
    }

    pub fn barcode(&self, cell_id: u64) -> String {
        IntToStr::from_u64(cell_id).to_string(self.cell_barcode_len)
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

        let feature_to_guide: HashMap<u64, u32> = self
            .feature_ids
            .iter()
            .enumerate()
            .map(|(guide_id, feature_id)| (*feature_id, guide_id as u32))
            .collect();

        let mut out: Vec<_> = cell
            .total_reads
            .iter()
            .filter_map(|(&feature_id, &value)| {
                if value <= 0.0 {
                    return None;
                }

                let guide_id = *feature_to_guide.get(&feature_id)?;
                Some(GuideObservation {
                    cell_id,
                    guide_id,
                    count: value.round() as u32,
                })
            })
            .collect();

        out.sort_unstable_by_key(|x| x.guide_id);
        out
    }
}
