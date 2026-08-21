use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

use scdata::FeatureIndex;

use crate::{
    AmbientModel,
    CellGuideAssignments,
    FittedModel,
    GuideCalls,
    MultiGuideGapStats,
    MultiGuideGapStatsTable,
};

pub struct BeaconResult {
    pub ambient: AmbientModel,
    pub fitted: FittedModel,
    pub calls: GuideCalls,
    pub assignments: CellGuideAssignments,
    pub multi_gap_stats: Vec<MultiGuideGapStats>,
}

impl BeaconResult {
    pub fn write<P, I>(
        &self,
        out: P,
        feature_index: &I,
        call_tag_len: usize,
    ) -> Result<()>
    where
        P: AsRef<Path>,
        I: FeatureIndex,
    {
        let out = out.as_ref().to_path_buf();
        //let filtered = feature_index.ordered_feature_ids().into_iter().map(|id| feature_index.feature_name(id).to_owned()).collect();

        fs::create_dir_all(&out)
            .with_context(|| {
                format!(
                    "creating Beacon output directory {}",
                    out.display()
                )
            })?;

        /*
         * ------------------------------------------------------------
         * Ambient model
         * ------------------------------------------------------------
         */
        self.ambient
            .write_table(
                &out,
                feature_index,
            )
            .context(
                "writing Beacon ambient model",
            )?;

        /*
         * ------------------------------------------------------------
         * Fitted model
         * ------------------------------------------------------------
         */
        self.fitted
            .write_table(
                &out,
                feature_index,
            )
            .context(
                "writing Beacon fitted model",
            )?;

        /*
         * ------------------------------------------------------------
         * Feature calls
         * ------------------------------------------------------------
         */
        self.calls
            .write_table(
                &out,
                feature_index,
                call_tag_len,
            )
            .context(
                "writing Beacon feature calls",
            )?;

        /*
         * ------------------------------------------------------------
         * Cell-level assignments
         * ------------------------------------------------------------
         */
        self.assignments
            .write_table(&out)
            .context(
                "writing Beacon cell assignments",
            )?;

        /*
         * ------------------------------------------------------------
         * Multi-feature statistics
         * ------------------------------------------------------------
         */
        self.multi_gap_stats
            .write_table(&out)
            .context(
                "writing Beacon multi-feature statistics",
            )?;

        /*
         * ------------------------------------------------------------
         * Human-readable summary
         * ------------------------------------------------------------
         */
        let main_log =
            self.multi_gap_stats
                .print_assignment_summary(
                    &self.assignments,
                );

        let feature_summary =
            self.multi_gap_stats
                .primary_guide_counts(
                    &self.assignments,
                    100.0,
                );

        let run_log = format!(
            "{}\n{}",
            main_log,
            feature_summary,
        );

        let log_path =
            out.join("sc_beacon.log");

        fs::write(
            &log_path,
            &run_log,
        )
        .with_context(|| {
            format!(
                "writing {}",
                log_path.display()
            )
        })?;

        Ok(())
    }
}