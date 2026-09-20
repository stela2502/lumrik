//! Fast seed-and-verify supplemental-feature mapper.
//!
//! `FastTagMapper` uses exact 16-bp seeds only to nominate feature positions,
//! then verifies the complete candidate against a concatenated forward/reverse
//! one-hot reference. Numeric FASTQ qualities can weight the final score.

pub mod builtin_tags;
pub mod cli;
pub mod fast_mapper;
pub mod feature_entry;
pub mod feature_index;
pub mod locus_mapper;
pub mod map_status;

pub use builtin_tags::{BuiltinTagSet, HUMAN_SAMPLE_TAGS, MOUSE_SAMPLE_TAGS};
pub use cli::FastMapperCli;
pub use fast_mapper::{AlignmentStrand, FastAlignment, FastTagMapper};
pub use feature_entry::FeatureEntry;
pub use feature_index::FastTagFeatureIndex;
pub use locus_mapper::FastLocusMapper;
pub use map_status::MapStatus;
