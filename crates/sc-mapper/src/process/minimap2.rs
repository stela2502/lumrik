use anyhow::{Context, Result};

use crate::core::{MapperLaunch, StreamingMapper};
use crate::process::{check_binary, MapperProcess};
use crate::traits::ExternalMapper;

#[derive(Debug, Clone)]
pub struct Minimap2 {
    launch: MapperLaunch,
}

impl Minimap2 {
    pub fn from_launch(launch: MapperLaunch) -> Self {
        Self { launch }
    }
}

impl ExternalMapper for Minimap2 {
    fn check(&self) -> Result<()> {
        check_binary(&self.launch.mapper_bin, "minimap2").with_context(|| {
            format!(
                "failed to validate minimap2 binary: {}",
                self.launch.mapper_bin.display()
            )
        })
    }

    fn spawn(&self) -> Result<StreamingMapper> {
        // minimap2 <user-options> <index> <r1> [r2]
        //
        // User must provide presets, e.g.:
        // --mapper-options "-ax map-ont"
        // --mapper-options "-ax sr"
        let mut args = self.launch.options.clone();

        args.push(self.launch.index.to_string_lossy().to_string());

        let process = if self.launch.paired {
            MapperProcess::spawn_paired_fifo(&self.launch.mapper_bin, &args, None)
                .context("failed to spawn minimap2 with paired FIFO input")?
        } else {
            args.push("-".to_string());

            MapperProcess::spawn_single_stdin(&self.launch.mapper_bin, &args, None)
                .context("failed to spawn minimap2 with single-end stdin input")?
        };

        Ok(StreamingMapper::new(Box::new(process)))
    }

    fn command_preview(&self) -> String {
        format!(
            "{} {} {} {}",
            self.launch.mapper_bin.display(),
            self.launch.options.join(" "),
            self.launch.index.display(),
            if self.launch.paired { "<R1_FIFO> <R2_FIFO>" } else { "-" }
        )
    }
}