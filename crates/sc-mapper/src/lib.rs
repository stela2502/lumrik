pub mod cli;
pub mod core;
pub mod process;
pub mod traits;

pub use cli::{MapperKind, StreamingMapperCli};
pub use core::{MappingCall, MappingClass, StreamingMapper};
pub use process::{Bwa, Minimap2, Star};

pub use traits::ExternalMapper;
