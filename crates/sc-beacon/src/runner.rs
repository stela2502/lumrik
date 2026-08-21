use std::collections::HashSet;

use anyhow::Result;
use scdata::{FeatureIndex, Scdata};

use crate::{
    AmbientModel,
    BackgroundConfig,
    CallConfig,
    CellGuideAssignments,
    FitConfig,
    GuideCalls,
    GuideDataset,
    MultiGuideGapStats,
    BeaconResult
};

use crate::model::{fit_mixture};


pub fn run_from_scdata<I>(
    raw: Scdata,
    cells: &HashSet<u64>,
    cell_barcode_len: usize,
    feature_index: &I,
    background_config: &BackgroundConfig,
    fit_config: &FitConfig,
    call_config: &CallConfig,
) -> Result<(BeaconResult, Scdata)>
where
    I: FeatureIndex,
{
    /*
     * ------------------------------------------------------------
     * Split raw counts into called cells and ambient/background
     * droplets.
     * ------------------------------------------------------------
     */
    let (
        filtered_scdata,
        background_scdata,
    ) = raw.split_by_cells(cells);

    /*
     * ------------------------------------------------------------
     * Convert Scdata into Beacon's internal representation.
     * ------------------------------------------------------------
     */
    let background =
        GuideDataset::from_scdata(
            background_scdata,
            cell_barcode_len
        )?;

    let filtered =
        GuideDataset::from_scdata(
            filtered_scdata,
            cell_barcode_len
        )?;

    /*
     * ------------------------------------------------------------
     * Ambient model
     * ------------------------------------------------------------
     */
    let ambient =
        AmbientModel::fit(
            &background,
            background_config,
        )?;

    /*
     * ------------------------------------------------------------
     * Ambient + true-feature mixture
     * ------------------------------------------------------------
     */
    let fitted =
        fit_mixture(
            &filtered,
            &ambient,
            fit_config,
        )?;

    /*
     * ------------------------------------------------------------
     * Calls
     * ------------------------------------------------------------
     */
    let calls =
        GuideCalls::from_model(
            &fitted,
            call_config,
        );

    /*
     * ------------------------------------------------------------
     * Cell-level assignments
     * ------------------------------------------------------------
     */
    let assignments =
        CellGuideAssignments::new(
            feature_index,
            &filtered,
            &calls,
        );

    /*
     * ------------------------------------------------------------
     * Multi-feature QC
     * ------------------------------------------------------------
     */
    let multi_gap_stats =
        MultiGuideGapStats::collect(
            &assignments,
        );

    Ok(
        (BeaconResult {
            ambient,
            fitted,
            calls,
            assignments,
            multi_gap_stats,
        },
        filtered.cells
        )
    )
}