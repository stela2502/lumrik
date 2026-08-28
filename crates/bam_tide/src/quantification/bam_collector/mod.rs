pub mod config;
pub mod collector;


pub use collector::{
    BamCollector,
    BamCollectorHandle,
};

pub use config::BamCollectorConfig;