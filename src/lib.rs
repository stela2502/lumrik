//! Fast BD/FASTA feature-id mapper.
//!
//! Layer model:
//! - `FeatureEntry`: clean FASTA/BD name + Scdata feature id
//! - `TagEntry`: one 8bp 2bit table entry pointing to a feature and a position
//! - public hot path: `map_feature_id(...) -> Option<u64>`
//!
//! A hit is accepted only when internal `min_hits` is surpassed.

pub mod builtin_tags;
pub mod cli;
pub mod fast_mapper;
pub mod feature_entry;
pub mod feature_index;
pub mod tag_entry;

pub use builtin_tags::{BuiltinTagSet, HUMAN_SAMPLE_TAGS, MOUSE_SAMPLE_TAGS};
pub use cli::FastMapperCli;
pub use fast_mapper::{encode_seq_positions_with_int_to_str, FastTagMapper, MapStatus, Slot};
pub use feature_entry::FeatureEntry;
pub use feature_index::FastTagFeatureIndex;
pub use tag_entry::TagEntry;
