use anyhow::Result;

use crate::core::{StreamingMapper};


pub trait ExternalMapper: std::fmt::Debug + Send + Sync {

    fn check(&self) -> Result<()>;
    fn spawn(&self) -> Result<StreamingMapper>;

    fn command_preview(&self) -> String;
}