pub mod background;
pub mod binary_counts;
pub mod caller;
pub mod cell_guide_assignments;
pub mod dataset;
pub mod guide_stats;
pub mod model;
pub mod stats;

mod beacon_result;
mod reporting;
mod runner;
mod utils;

pub use background::{AmbientModel, BackgroundConfig};
pub use beacon_result::BeaconResult;
pub use binary_counts::{
    BinaryCountFit, BinaryCountFitConfig, fit_binary_counts, fit_binary_counts_with_config,
};
pub use caller::{CallConfig, GuideCall, GuideCalls};
pub use cell_guide_assignments::{CellGuideAssignment, CellGuideAssignments, GuideEvidence};
pub use dataset::{GuideDataset, GuideObservation};
pub use guide_stats::{MultiGuideGapStats, MultiGuideGapStatsTable};
pub use model::{FitConfig, FittedModel, GuideExpressionModel};
pub use runner::run_from_scdata;
