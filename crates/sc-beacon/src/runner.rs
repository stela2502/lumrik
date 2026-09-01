use std::collections::HashSet;

use anyhow::Result;
use mapping_info::MappingInfo;
use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};

use crate::model::fit_mixture;
use crate::{
    AmbientModel, BackgroundConfig, BeaconResult, CallConfig, CellGuideAssignments, FitConfig,
    GuideCalls, GuideDataset, MultiGuideGapStats,
};

/// Run Beacon directly on lumrik-native sparse data.
///
/// `filtered` contains feature counts for retained cells and `background`
/// contains non-cell/background droplets. `cells` is the complete canonical
/// retained-cell set, including cells with zero feature counts. Cell splitting
/// is deliberately owned by the caller; Beacon only performs inference.
/// The supplied `FeatureIndex` defines the real feature ids and their canonical
/// order; Beacon uses dense model-local slots internally without assuming ids
/// are contiguous.
pub fn run_from_scdata<I>(
    filtered: &Scdata,
    background: &Scdata,
    cells: &HashSet<u64>,
    cell_barcode_len: usize,
    feature_index: &I,
    background_config: &BackgroundConfig,
    fit_config: &FitConfig,
    call_config: &CallConfig,
) -> Result<BeaconResult>
where
    I: FeatureIndex,
{
    let background = GuideDataset::from_scdata(background, feature_index, cell_barcode_len)?;

    let filtered = GuideDataset::from_scdata_with_cells(
        filtered,
        cells.iter().copied().collect(),
        feature_index,
        cell_barcode_len,
    )?;

    let ambient = AmbientModel::fit(&background, background_config)?;
    let fitted = fit_mixture(&filtered, &ambient, fit_config)?;
    let calls = GuideCalls::from_model(&fitted, call_config);

    let assignments = CellGuideAssignments::new(feature_index, &filtered, &calls);

    let multi_gap_stats = MultiGuideGapStats::collect(&assignments);

    // Native posterior matrix: same cell ids and real feature ids as the input
    // Scdata. Missing cell/feature combinations are implicit zeros.
    let mut posteriors = Scdata::new(rayon::current_num_threads().max(1), MatrixValueType::Real);
    let mut report = MappingInfo::new(None, 0.0, 0);

    for call in &calls.flat {
        if call.posterior.probability <= 0.0 {
            continue;
        }

        posteriors.try_insert_value(
            &call.cell_id,
            GeneUmiHash(call.feature_id, 0),
            call.posterior.probability as f32,
            &mut report,
        );
    }

    Ok(BeaconResult {
        ambient,
        fitted,
        calls,
        assignments,
        multi_gap_stats,
        posteriors,
    })
}
