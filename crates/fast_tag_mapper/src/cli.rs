use clap::{Args};
use std::path::PathBuf;

use crate::{BuiltinTagSet, FastTagMapper};

#[derive(Debug, Clone, Args)]
pub struct FastMapperCli {
    /// Built-in BD sample tags
    #[arg(long, value_enum)]
    pub species: Option<BuiltinTagSet>,

    /// Additional FASTA file containing feature tags
    #[arg(long)]
    pub tags: Option<PathBuf>,

    /// Minimum supporting 8-mer hits
    #[arg(long, default_value_t = 4)]
    pub min_hits: u32,
}

impl FastMapperCli {
    pub fn mapper(&self) -> anyhow::Result<FastTagMapper> {
        let mut mapper = match self.species {
            Some(BuiltinTagSet::Mouse) => FastTagMapper::mouse_samples(),
            Some(BuiltinTagSet::Human) => FastTagMapper::human_samples(),
            None => FastTagMapper::new(),
        };

        if let Some(path) = &self.tags {
            mapper.load_fasta(path)?;
        }

        Ok(mapper.with_min_hits(self.min_hits))
    }
}