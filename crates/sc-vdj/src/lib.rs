//! Posterior single-cell V(D)J analyzer for Nelrune.
//! Runs after normal mapping/quantification and combines the retained BAM with
//! the complete per-cell expression matrix.
pub mod align;
pub mod audit;
pub mod bam;
pub mod gex;
pub mod mapper;
pub mod output;
pub mod posterior;
pub mod reference;
pub mod score;
pub mod sequence;
pub mod sterile;
pub mod types;

pub use bam::{
    bam_shard_for_cell, read_bam, read_bam_filtered, shard_bam_receptor_evidence,
    AuxTagIdentityResolver, BamIdentityResolver, BamShardStats, NelruneIdentityResolver,
    QnameIdentityResolver,
};
pub use gex::{ExpressionMatrix, LongTsvExpression};
pub use mapper::{VdjMapper, VdjMapperConfig};
pub use posterior::{
    BamReadEvidence, CellVdjSummary, GermlineSegmentSupport, PosteriorAnalyzer, PosteriorConfig,
    ReadSegmentAlignment, RearrangementCall, RearrangementSupportingRead, RecombinationStage,
};
pub use reference::{VdjReference, VdjReferenceBuilder};
pub use score::{decode_evidence_code, DevelopmentEvidence, DevelopmentProgram};
pub use sterile::{SterileBin, SterileProfile, SupportedInterval};
pub use types::{Chain, Orientation, SegmentHit, SegmentKind, VdjCandidate, VdjSegment};
