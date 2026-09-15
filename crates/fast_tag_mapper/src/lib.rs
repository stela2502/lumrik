//! Fast exact supplemental-feature mapper.
//!
//! `FastTagMapper` indexes exact 16-bp seeds for dynamic features such as BD
//! sample tags, guides, antibodies, GFP, and other sequences that should not
//! require rebuilding the genomic mapper index.

pub mod builtin_tags;
pub mod cli;
pub mod fast_mapper;
pub mod feature_entry;
pub mod feature_index;
pub mod locus_mapper;
pub mod map_status;

pub use builtin_tags::{BuiltinTagSet, HUMAN_SAMPLE_TAGS, MOUSE_SAMPLE_TAGS};
pub use cli::FastMapperCli;
pub use fast_mapper::FastTagMapper;
pub use feature_entry::FeatureEntry;
pub use feature_index::FastTagFeatureIndex;
pub use locus_mapper::FastLocusMapper;
pub use map_status::MapStatus;
