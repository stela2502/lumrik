use anyhow::{Context, Result};
use bam_tide::fastq::{FastqPairReader, FastqRead, SimpleFastqReader};
use clap::Parser;
use rust_htslib::bam;
use sc_mapper::{MappingCall, StreamingMapper, StreamingMapperCli};
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
struct Cli {
    #[command(flatten)]
    mapper: StreamingMapperCli,

    #[arg(long)]
    r1: PathBuf,

    #[arg(long)]
    r2: Option<PathBuf>,

    /// Optional debug table, one line per received mapping call.
    #[arg(long)]
    out_tsv: Option<PathBuf>,

    /// BAM output containing all SAM/BAM records returned by the mapper.
    #[arg(long)]
    out_bam: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let mut mapper = cli.mapper.from_cli()?;

    let mut out_tsv = match &cli.out_tsv {
        Some(path) => {
            let mut out = std::fs::File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?;

            writeln!(out, "read_id\tn_records")?;

            Some(out)
        }

        None => None,
    };

    /*
     * IMPORTANT:
     *
     * We cannot create the BAM writer yet.
     *
     * Its header must come from the mapper output, and that header
     * only becomes available after the mapper stdout reader has
     * initialized its BAM/SAM reader.
     */
    let mut bam_writer: Option<bam::Writer> = None;

    match &cli.r2 {
        Some(r2_path) => {
            run_paired(
                &mut mapper,
                &cli.r1,
                r2_path,
                out_tsv.as_mut(),
                &mut bam_writer,
                &cli.out_bam,
            )?;
        }

        None => {
            run_single(
                &mut mapper,
                &cli.r1,
                out_tsv.as_mut(),
                &mut bam_writer,
                &cli.out_bam,
            )?;
        }
    }

    mapper.finish()?;

    Ok(())
}

fn run_single(
    mapper: &mut StreamingMapper,
    r1_path: &Path,
    mut out_tsv: Option<&mut std::fs::File>,
    bam_writer: &mut Option<bam::Writer>,
    out_bam: &Path,
) -> Result<()> {
    let mut reads = SimpleFastqReader::new(&r1_path.to_path_buf())
        .with_context(|| format!("failed to open R1 FASTQ {}", r1_path.display()))?;

    while let Some(r1) = reads.next_record()? {
        /*
         * Keep this separate from write_mapping_call().
         *
         * mapper.process() needs mutable access to mapper, and
         * write_mapping_call() may also need mutable access when
         * obtaining/caching the header.
         */
        let call = mapper.process(&r1, None)?;

        write_mapping_call(call, out_tsv.as_deref_mut(), bam_writer, out_bam, mapper)?;
    }

    Ok(())
}

fn run_paired(
    mapper: &mut StreamingMapper,
    r1_path: &Path,
    r2_path: &Path,
    mut out_tsv: Option<&mut std::fs::File>,
    bam_writer: &mut Option<bam::Writer>,
    out_bam: &Path,
) -> Result<()> {
    let mut reads = FastqPairReader::from_paths(&r1_path.to_path_buf(), &r2_path.to_path_buf())
        .with_context(|| {
            format!(
                "failed to open paired FASTQs R1={} R2={}",
                r1_path.display(),
                r2_path.display()
            )
        })?;

    while let Some((r1, r2)) = reads.next_pair()? {
        let call = mapper.process(&r1, Some(&r2))?;

        write_mapping_call(call, out_tsv.as_deref_mut(), bam_writer, out_bam, mapper)?;
    }

    Ok(())
}

fn write_mapping_call(
    call: Option<MappingCall>,
    mut out_tsv: Option<&mut std::fs::File>,
    bam_writer: &mut Option<bam::Writer>,
    out_bam: &Path,
    mapper: &mut StreamingMapper,
) -> Result<()> {
    let Some(call) = call else {
        return Ok(());
    };

    if let Some(out) = out_tsv.as_deref_mut() {
        writeln!(out, "{}\t{}", call.read_id, call.records.records.len())?;
    }

    /*
     * Create the BAM writer lazily.
     *
     * If we already received a mapping record, the stdout BAM
     * reader necessarily parsed the mapper header first.
     */
    if bam_writer.is_none() {
        let header = mapper
            .header()
            .context("mapper produced records but no SAM/BAM header was available")?
            .clone();

        *bam_writer = Some(
            bam::Writer::from_path(out_bam, &header, bam::Format::Bam)
                .with_context(|| format!("failed to create BAM {}", out_bam.display()))?,
        );
    }

    let writer = bam_writer
        .as_mut()
        .expect("BAM writer was initialized above");

    for rec in &call.records.records {
        writer.write(&rec.record)?;
    }

    Ok(())
}
