pub mod collector;
pub mod config;
pub mod read_group;

pub use collector::{BamCollector, BamCollectorHandle};

pub use config::BamCollectorConfig;
