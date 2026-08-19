
use anyhow::{Context, Result};
use std::path::PathBuf;
use clap::{Args, ValueEnum};

use crate::Bwa;
use crate::Minimap2;
use crate::Star;

use crate::core::{StreamingMapper,MapperLaunch};
use crate::traits::ExternalMapper;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MapperKind {
    Minimap2,
    Star,
    Bwa,
}

#[derive(Debug, Clone, Args)]
pub struct StreamingMapperCli {
    /// External mapper to use.
    #[arg(long = "mapper", value_enum)]
    pub mapper: MapperKind,

    /// Mapper binary.
    ///
    /// If omitted, defaults to:
    /// minimap2, STAR, or bwa.
    #[arg(long = "mapper-bin")]
    pub mapper_bin: Option<PathBuf>,

    /// Mapper index path.
    ///
    /// For minimap2 this is the reference FASTA or .mmi index.
    /// For STAR this is the genomeDir.
    /// For BWA this is the indexed reference prefix / FASTA.
    #[arg(long = "mapper-index")]
    pub mapper_index: PathBuf,

    /// Additional mapper options.
    ///
    /// These are appended after the mapper defaults.
    ///
    /// Example:
    /// --mapper-options "-t 8 --secondary=yes"
    #[arg(long = "mapper-options")]
    pub mapper_options: Option<String>,

    /// Number of mapper threads.
    ///
    /// This is translated into mapper-specific thread arguments.
    #[arg(long = "mapper-threads", default_value_t = 1)]
    pub mapper_threads: usize,

    /// Input layout.
    ///
    /// Single-end mappers receive one FASTQ stream.
    /// Paired-end mappers receive two FASTQ streams where supported.
    #[arg(long = "mapper-paired", default_value_t = false)]
    pub mapper_paired: bool,

    /// Keep secondary / supplementary / multimapping records where supported.
    #[arg(long = "mapper-keep-multimappers", default_value_t = true)]
    pub mapper_keep_multimappers: bool,
}

impl StreamingMapperCli {
    pub fn launch(&self) -> Result<MapperLaunch> {
        let mapper_bin = self.mapper_bin.clone().unwrap_or_else(|| {
            match self.mapper {
                MapperKind::Minimap2 => "minimap2".into(),
                MapperKind::Star => "STAR".into(),
                MapperKind::Bwa => "bwa".into(),
            }
        });

        Ok(MapperLaunch {
            mapper_bin,
            index: self.mapper_index.clone(),
            threads: self.mapper_threads,
            paired: self.mapper_paired,
            options: match &self.mapper_options {
                Some(x) => split_mapper_options(x)
                    .with_context(|| format!("Failed to parse --mapper-options: {x:?}"))?,
                None => Vec::new(),
            },
        })
    }

    pub fn from_cli(&self) -> Result<StreamingMapper> {
        let launch = self.launch()?;

        let mapper: Box<dyn ExternalMapper> = match self.mapper {
            MapperKind::Minimap2 => Box::new(Minimap2::from_launch(launch)),
            MapperKind::Star => Box::new(Star::from_launch(launch)),
            MapperKind::Bwa => Box::new(Bwa::from_launch(launch)),
        };

        mapper.check()?;
        mapper.spawn()
    }
}

fn split_mapper_options(options: &str) -> Result<Vec<String>> {
    shell_words::split(options)
        .with_context(|| format!("Invalid mapper options: {options:?}"))
}