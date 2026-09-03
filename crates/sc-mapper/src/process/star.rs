use anyhow::{Context, Result, bail};

use crate::core::{MapperLaunch, StreamingMapper};
use crate::process::{MapperProcess, check_binary, has_option, remove_option};
use crate::traits::ExternalMapper;

use rust_htslib::bam::{Header, header::HeaderRecord};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Star {
    launch: MapperLaunch,
    header: Header,
}

impl Star {
    pub fn from_launch(launch: MapperLaunch) -> Result<Self> {
        let header = Self::header_from_star_index(&launch.index)?;

        Ok(Self { launch, header })
    }

    fn header_from_star_index(genome_dir: &Path) -> Result<Header> {
        let names_path = genome_dir.join("chrName.txt");

        let lengths_path = genome_dir.join("chrLength.txt");

        let names = fs::read_to_string(&names_path).with_context(|| {
            format!(
                "reading STAR chromosome names from {}",
                names_path.display()
            )
        })?;

        let lengths = fs::read_to_string(&lengths_path).with_context(|| {
            format!(
                "reading STAR chromosome lengths from {}",
                lengths_path.display()
            )
        })?;

        let names: Vec<&str> = names.lines().filter(|line| !line.is_empty()).collect();

        let lengths: Vec<u64> = lengths
            .lines()
            .filter(|line| !line.is_empty())
            .map(|line| {
                line.parse::<u64>()
                    .with_context(|| format!("invalid STAR chromosome length `{line}`"))
            })
            .collect::<Result<_>>()?;

        if names.len() != lengths.len() {
            bail!(
                "STAR index is inconsistent: {} chromosome names but {} lengths",
                names.len(),
                lengths.len(),
            );
        }

        let mut header = Header::new();

        for (name, length) in names.iter().zip(lengths.iter()) {
            let mut sq = HeaderRecord::new(b"SQ");

            sq.push_tag(b"SN", name);

            sq.push_tag(b"LN", length);

            header.push_record(&sq);
        }

        Ok(header)
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

        remove_option(&mut args, "--runThreadN", Some(1));
        remove_option(&mut args, "--genomeDir", Some(1));
        remove_option(
            &mut args,
            "--readFilesIn",
            Some(if self.launch.paired { 2 } else { 1 }),
        );
        remove_option(&mut args, "--outSAMtype", None);
        remove_option(&mut args, "--outStd", Some(1));

        /*
         * MapperProcess runs FIFO-based mappers in a private temporary
         * working directory so concurrent STAR processes do not collide
         * on files such as ./_STARtmp. Resolve a relative genome index
         * against the caller's working directory before changing cwd.
         */
        let genome_dir = if self.launch.index.is_absolute() {
            self.launch.index.clone()
        } else {
            std::env::current_dir()
                .context("getting current directory for relative STAR genome index")?
                .join(&self.launch.index)
        };

        args.extend([
            "--runThreadN".into(),
            self.launch.threads.min(4).to_string(),
            "--genomeDir".into(),
            genome_dir.to_string_lossy().into_owned(),
            "--outSAMtype".into(),
            "BAM".into(),
            "Unsorted".into(),
            "--outStd".into(),
            "BAM_Unsorted".into(),
        ]);

        if !has_option(&args, "--outBAMcompression") {
            args.extend(["--outBAMcompression".into(), "0".into()]);
        }

        if !has_option(&args, "--outSAMattributes") {
            args.extend([
                "--outSAMattributes".into(),
                "NH".into(),
                "HI".into(),
                "AS".into(),
                "nM".into(),
            ]);
        }

        if !has_option(&args, "--outReadsUnmapped") {
            args.extend(["--outReadsUnmapped".into(), "None".into()]);
        }

        args.push("--readFilesIn".into());

        let process = if self.launch.paired {
            MapperProcess::spawn_paired_fifo(
                &self.launch.mapper_bin,
                &args,
                Some(self.header.clone()),
            )
            .context("failed to spawn STAR with paired FIFO input")?
        } else {
            MapperProcess::spawn_single_fifo(
                &self.launch.mapper_bin,
                &args,
                Some(self.header.clone()),
            )
            .context("failed to spawn STAR with single FIFO input")?
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
