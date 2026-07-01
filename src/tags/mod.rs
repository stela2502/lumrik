pub mod cli;
pub mod fast_tag_mapper;

pub use cli::{BuiltinTagSet, TagCli, TagRead};
pub use fast_tag_mapper::{FastTagMapper, TagCall, TagEntry};
