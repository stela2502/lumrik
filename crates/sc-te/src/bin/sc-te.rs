use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use sc_te::convert_ucsc_rmsk_to_gtf;

#[derive(Debug, Parser)]
#[command(name = "sc-te", about = "Utilities for single-cell transposable-element analysis")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Convert an official UCSC rmsk table to the TE GTF used by Lumrik.
    UcscRmskToGtf {
        /// UCSC rmsk.txt or rmsk.txt.gz file.
        #[arg(long)]
        input: PathBuf,

        /// Output .gtf or .gtf.gz file.
        #[arg(long)]
        output: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Command::UcscRmskToGtf { input, output } => {
            let reader = open_reader(&input)?;
            let writer = open_writer(&output)?;
            let summary = convert_ucsc_rmsk_to_gtf(reader, writer)?;
            eprintln!(
                "converted {} UCSC RepeatMasker rows to {} GTF rows: {}",
                summary.input_rows,
                summary.output_rows,
                output.display()
            );
        }
    }

    Ok(())
}

fn open_reader(path: &Path) -> Result<BufReader<Box<dyn Read>>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader: Box<dyn Read> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::new(reader))
}

fn open_writer(path: &Path) -> Result<BufWriter<Box<dyn Write>>> {
    let file = File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer: Box<dyn Write> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(file)
    };
    Ok(BufWriter::new(writer))
}
