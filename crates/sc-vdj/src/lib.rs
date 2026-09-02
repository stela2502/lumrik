//! Posterior single-cell V(D)J analyzer for Nelrune.
//! Runs after normal mapping/quantification and combines the retained BAM with
//! the complete per-cell expression matrix.
pub mod align;
pub mod audit;
pub mod bam;
pub mod gex;
pub mod identity;
pub mod junction;
pub mod mapper;
pub mod output;
pub mod posterior;
pub mod reference;
pub mod score;
pub mod sequence;
pub mod sterile;
pub mod types;

pub use bam::{
    read_bam, read_bam_filtered, read_bam_receptor_evidence,
    read_bam_receptor_evidence_with_progress, AuxTagIdentityResolver, BamEvidenceStats,
    BamIdentityResolver, NelruneIdentityResolver, QnameIdentityResolver, RoutedBamEvidence,
};
pub use gex::{ExpressionMatrix, LongTsvExpression};
pub use identity::{
    LightChainStatus, PackedRecombinationId, ReceptorRole, RecombinationMeasurements,
};
pub use junction::{JunctionInput, JunctionMeasurement};
pub use mapper::{VdjMapper, VdjMapperConfig};
pub use posterior::{
    BamReadEvidence, CellVdjSummary, GermlineSegmentSupport, PosteriorAnalyzer, PosteriorConfig,
    ReadSegmentAlignment, RearrangementCall, RearrangementSupportingRead, RecombinationStage,
};
pub use reference::{VdjReference, VdjReferenceBuilder};
pub use score::{MarkerContribution, RecombinationActivityEvidence};
pub use sterile::{SterileBin, SterileProfile, SupportedInterval};
pub use types::{Chain, Orientation, SegmentHit, SegmentKind, VdjCandidate, VdjSegment};
