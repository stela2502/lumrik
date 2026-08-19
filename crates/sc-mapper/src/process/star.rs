use anyhow::{Context, Result};

use crate::core::{MapperLaunch, StreamingMapper};
use crate::process::{check_binary, MapperProcess};
use crate::traits::ExternalMapper;

#[derive(Debug, Clone)]
pub struct Star {
    launch: MapperLaunch,
}

impl Star {
    pub fn from_launch(launch: MapperLaunch) -> Self {
        Self { launch }
    }
}

impl ExternalMapper for Star {
    fn check(&self) -> Result<()> {
        check_binary(&self.launch.mapper_bin, "STAR").with_context(|| {
            format!(
                "failed to validate STAR binary: {}",
                self.launch.mapper_bin.display()
            )
        })
    }

    fn spawn(&self) -> Result<StreamingMapper> {
        /*
         * STAR mapping command:
         *
         * STAR <user-options>
         *      --genomeDir <index>
         *      --readFilesIn <r1> [r2]
         *
         * Unlike minimap2/BWA, STAR does not treat "-"
         * as stdin for --readFilesIn.
         *
         * Therefore:
         *
         *   single-end -> one named FIFO
         *   paired-end -> two named FIFOs
         *
         * MapperProcess owns the FIFO writers and streams FASTQ
         * records into them while STAR reads them as normal files.
         */

        let mut args = self.launch.options.clone();

        args.push("--genomeDir".to_string());
        args.push(
            self.launch
                .index
                .to_string_lossy()
                .to_string(),
        );

        args.push("--readFilesIn".to_string());

        let process = if self.launch.paired {
            MapperProcess::spawn_paired_fifo(
                &self.launch.mapper_bin,
                &args,
            )
            .context(
                "failed to spawn STAR with paired FIFO input",
            )?
        } else {
            MapperProcess::spawn_single_fifo(
                &self.launch.mapper_bin,
                &args,
            )
            .context(
                "failed to spawn STAR with single FIFO input",
            )?
        };

        Ok(StreamingMapper::new(Box::new(process)))
    }

    fn command_preview(&self) -> String {
        let reads = if self.launch.paired {
            "<R1_FIFO> <R2_FIFO>"
        } else {
            "<R1_FIFO>"
        };

        format!(
            "{} {} --genomeDir {} --readFilesIn {}",
            self.launch.mapper_bin.display(),
            self.launch.options.join(" "),
            self.launch.index.display(),
            reads,
        )
    }
}