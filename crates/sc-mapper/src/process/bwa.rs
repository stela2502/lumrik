use anyhow::{Context, Result};

use crate::core::{MapperLaunch, StreamingMapper};
use crate::process::{check_binary, MapperProcess};
use crate::traits::ExternalMapper;

#[derive(Debug, Clone)]
pub struct Bwa {
    launch: MapperLaunch,
}

impl Bwa {
    pub fn from_launch(launch: MapperLaunch) -> Self {
        Self { launch }
    }
}

impl ExternalMapper for Bwa {
    fn check(&self) -> Result<()> {
        check_binary(&self.launch.mapper_bin, "bwa").with_context(|| {
            format!(
                "failed to validate bwa binary: {}",
                self.launch.mapper_bin.display()
            )
        })
    }

    fn spawn(&self) -> Result<StreamingMapper> {
        // bwa <user-options> <index> <r1> [r2]
        //
        // User must provide the subcommand, e.g.:
        // --mapper-options "mem -t 8"
        let mut args = self.launch.options.clone();

        args.push(self.launch.index.to_string_lossy().to_string());

        let process = if self.launch.paired {
            MapperProcess::spawn_paired_fifo(&self.launch.mapper_bin, &args)
                .context("failed to spawn bwa with paired FIFO input")?
        } else {
            args.push("-".to_string());

            MapperProcess::spawn_single_stdin(&self.launch.mapper_bin, &args)
                .context("failed to spawn bwa with single-end stdin input")?
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