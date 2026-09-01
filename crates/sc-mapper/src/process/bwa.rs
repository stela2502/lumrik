use anyhow::{Context, Result};

use crate::core::{MapperLaunch, StreamingMapper};
use crate::process::{MapperProcess, check_binary, remove_option};
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
                "failed to validate BWA binary: {}",
                self.launch.mapper_bin.display()
            )
        })
    }

    fn spawn(&self) -> Result<StreamingMapper> {
        // bwa mem <options> <index> <r1> [r2]
        //
        // The user controls mapping behavior.
        // sc-mapper controls the thread count and input paths.

        let mut args = self.launch.options.clone();

        remove_option(&mut args, "-t", Some(1));

        args.extend([
            "mem".into(),
            "-t".into(),
            self.launch.threads.to_string(),
            self.launch.index.to_string_lossy().into_owned(),
        ]);

        let process = if self.launch.paired {
            MapperProcess::spawn_paired_fifo(&self.launch.mapper_bin, &args, None)
                .context("failed to spawn BWA with paired FIFO input")?
        } else {
            args.push("-".into());

            MapperProcess::spawn_single_stdin(&self.launch.mapper_bin, &args, None)
                .context("failed to spawn BWA with single-end stdin input")?
        };

        Ok(StreamingMapper::new(Box::new(process)))
    }

    fn command_preview(&self) -> String {
        format!(
            "{} mem {} -t {} {} {}",
            self.launch.mapper_bin.display(),
            self.launch.options.join(" "),
            self.launch.threads,
            self.launch.index.display(),
            if self.launch.paired {
                "<R1_FIFO> <R2_FIFO>"
            } else {
                "-"
            },
        )
    }
}
