//! Clean single-cell V(D)J reconstruction for Lumrik/Nelrune.
//!
//! The crate has four biological layers:
//! - `index`: immutable V/D/J/C reference and genomic overlap lookup
//! - `cellrep`: BAM-derived per-cell sequence evidence
//! - `recombination`: receptor assembly and V(D)J junction inference
//! - `output`: AIRR-compatible serialization
//! `runner` is deliberately glue only.

pub mod cellrep;
pub mod index;
pub mod output;
pub mod recombination;
pub mod runner;

pub use cellrep::{
    AlignmentGeometry, BamFeatureEvidence, BamFeatureSequenceParts, CellEvidence, CellEvidenceVdj,
    EvidenceId, MapperEvidence, SequencePart,
};
pub use index::{Chain, SegmentId, SegmentKind, Strand, VdjIndex, VdjIndexBuilder, VdjSegment};
pub use recombination::{
    ConstantRegionEvidence, JunctionStructure, ProductivityStatus, ReceptorRole, Recombination,
    RecombinationId,
};
pub use runner::{
    BamIdentityResolver, BamIngestProgress, NelruneIdentityResolver, VdjRunner, VdjRunnerConfig,
};
