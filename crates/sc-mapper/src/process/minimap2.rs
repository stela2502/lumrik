use anyhow::{Context, Result};

use crate::core::{MapperLaunch, StreamingMapper};
use crate::process::{
    check_binary,
    has_option,
    remove_option,
    MapperProcess,
};
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
        // minimap2 <options> <index> <r1> [r2]
        //
        // The user controls the mapping preset, e.g.:
        //
        //   -x sr
        //   -x map-ont
        //
        // sc-mapper controls threading and requires SAM output.

        let mut args = self.launch.options.clone();

        // sc-mapper owns the thread count.
        remove_option(&mut args, "-t", Some(1));

        // PAF output is incompatible with StreamingMapper.
        remove_option(&mut args, "-c", Some(0));

        args.extend([
            "-t".into(),
            self.launch.threads.to_string(),
        ]);

        // Force SAM output unless the user already requested it.
        if !has_option(&args, "-a") {
            args.push("-a".into());
        }

        args.push(
            self.launch
                .index
                .to_string_lossy()
                .into_owned(),
        );

        let process = if self.launch.paired {
            MapperProcess::spawn_paired_fifo(
                &self.launch.mapper_bin,
                &args,
                None,
            )
            .context(
                "failed to spawn minimap2 with paired FIFO input",
            )?
        } else {
            args.push("-".into());

            MapperProcess::spawn_single_stdin(
                &self.launch.mapper_bin,
                &args,
                None,
            )
            .context(
                "failed to spawn minimap2 with single-end stdin input",
            )?
        };

        Ok(StreamingMapper::new(Box::new(process)))
    }

    fn command_preview(&self) -> String {
        format!(
            "{} {} -t {} -a {} {}",
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