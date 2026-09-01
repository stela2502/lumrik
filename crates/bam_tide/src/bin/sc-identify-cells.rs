// bam_tide/src/bin/sc-identify-cells.rs

use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use int_to_str::IntToStr;
use sc_primer::{PrimerCli, PrimerDetector, PrimerMatch};
use scdata::cell_data::GeneUmiHash;

use bam_tide::fastq::record::FastqRecord;
use bam_tide::fastq::writer::FastqWriter;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Identify one insert signature in single-cell FASTQ reads"
)]
struct Cli {
    #[command(flatten)]
    primer: PrimerCli,

    /// Input FASTQ or FASTQ.GZ .
    #[arg(short, long, required =true, num_args = 1..)]
    fastq: Vec<PathBuf>,

    /// Output TSV.
    #[arg(short, long)]
    out: PathBuf,

    /// Feature id used as GeneUmiHash(feature_id, umi_id).
    /// This is totally user defined and represents the INSERT you search for.
    #[arg(long, default_value_t = 0)]
    insert_id: u64,

    /// Print progress every N reads. Set to 0 to disable.
    #[arg(long, default_value_t = 1_000_000)]
    progress_every: usize,

    /// Optional paired R2 FASTQ or FASTQ.GZ.
    /// Must have the same number/order of records as --fastq.
    #[arg(long, num_args = 1..)]
    r2_fastq: Vec<PathBuf>,

    /// Write matched R1 records.
    #[arg(long)]
    r1_out: Option<PathBuf>,

    /// Output FASTQ containing R2 records whose paired R1 matched.
    #[arg(long)]
    r2_out: Option<PathBuf>,

    /// Gzip level for written fastq files.
    #[arg(long, default_value_t = 1)]
    out_gzip_level: u32,

    /// Gzip written fastq files.
    #[arg(long, default_value_t = false)]
    out_gzip: bool,
}

#[derive(Debug, Default)]
struct CellHitSet {
    reads: usize,
    molecules: HashSet<GeneUmiHash>,
}

#[derive(Debug, Default)]
struct Stats {
    total: usize,
    matched: usize,
    failed: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let detector = cli
        .primer
        .detector()
        .map_err(anyhow::Error::msg)
        .context("failed to build primer detector")?;

    if !cli.r2_fastq.is_empty() && cli.r2_fastq.len() != cli.fastq.len() {
        anyhow::bail!("--r2-fastq must contain the same number of files as --fastq");
    }

    if cli.r2_out.is_some() && cli.r2_fastq.is_empty() {
        anyhow::bail!("--r2-out requires --r2-fastq");
    }

    let mut r1_writer = match &cli.r1_out {
        Some(path) => Some(FastqWriter::new(
            path,
            cli.out_gzip || path.extension().is_some_and(|e| e == "gz"),
            cli.out_gzip_level,
        )?),
        None => None,
    };

    let mut r2_writer = match &cli.r2_out {
        Some(path) => Some(FastqWriter::new(
            path,
            cli.out_gzip || path.extension().is_some_and(|e| e == "gz"),
            cli.out_gzip_level,
        )?),
        None => None,
    };

    let mut stats = Stats::default();
    let mut cells: HashMap<u64, CellHitSet> = HashMap::new();

    for (i, fastq) in cli.fastq.iter().enumerate() {
        eprintln!("reading FASTQ {}", fastq.display());

        let mut reader = open_fastq(fastq)
            .with_context(|| format!("failed to open FASTQ {}", fastq.display()))?;

        let mut r2_reader = if cli.r2_fastq.is_empty() {
            None
        } else {
            let r2_fastq = &cli.r2_fastq[i];

            eprintln!("reading paired R2 FASTQ {}", r2_fastq.display());

            Some(
                open_fastq(r2_fastq)
                    .with_context(|| format!("failed to open R2 FASTQ {}", r2_fastq.display()))?,
            )
        };

        while let Some(record) = read_fastq_record(&mut reader)? {
            let r2_record = if let Some(r2_reader) = r2_reader.as_mut() {
                Some(
                    read_fastq_record(r2_reader.as_mut())?
                        .context("R2 FASTQ ended before R1 FASTQ")?,
                )
            } else {
                None
            };

            stats.total += 1;

            match process_record(&detector, &record, cli.insert_id)? {
                Some((cell_id, gene_umi)) => {
                    stats.matched += 1;

                    let entry = cells.entry(cell_id).or_default();
                    entry.reads += 1;
                    entry.molecules.insert(gene_umi);

                    if let Some(writer) = r1_writer.as_mut() {
                        writer.write(&record)?;
                    }

                    if let (Some(writer), Some(r2)) = (r2_writer.as_mut(), r2_record.as_ref()) {
                        writer.write(r2)?;
                    }
                }
                None => {
                    stats.failed += 1;
                }
            }

            if cli.progress_every > 0 && stats.total % cli.progress_every == 0 {
                let matched_pct = 100.0 * stats.matched as f64 / stats.total as f64;
                let failed_pct = 100.0 * stats.failed as f64 / stats.total as f64;

                eprintln!(
                    "processed={} matched={} ({:.2}%) failed={} ({:.2}%) cells={}",
                    stats.total,
                    stats.matched,
                    matched_pct,
                    stats.failed,
                    failed_pct,
                    cells.len()
                );
            }
        }

        if let Some(r2_reader) = r2_reader.as_mut() {
            if read_fastq_record(r2_reader.as_mut())?.is_some() {
                anyhow::bail!(
                    "R2 FASTQ {} contains more records than R1 FASTQ {}",
                    cli.r2_fastq[i].display(),
                    fastq.display()
                );
            }
        }
    }

    write_output(&cli.out, cli.insert_id, &cells)?;

    eprintln!(
        "done: total={} matched={} ({:.2}%) failed={} ({:.2}%) cells={}",
        stats.total,
        stats.matched,
        100.0 * stats.matched as f64 / stats.total.max(1) as f64,
        stats.failed,
        100.0 * stats.failed as f64 / stats.total.max(1) as f64,
        cells.len()
    );

    if let Some(writer) = r1_writer {
        writer.finish()?;
    }

    if let Some(writer) = r2_writer {
        writer.finish()?;
    }

    Ok(())
}

fn process_record(
    detector: &PrimerDetector,
    record: &FastqRecord,
    insert_id: u64,
) -> Result<Option<(u64, GeneUmiHash)>> {
    let Some(hit) = detector
        .detect_first(&record.seq, &record.qual)
        .map_err(anyhow::Error::msg)?
    else {
        return Ok(None);
    };

    let Some(cell_id) = hit.bd_cell_id else {
        return Ok(None);
    };

    let Some(umi_seq) = first_segment_seq(record, &hit, "UMI") else {
        return Ok(None);
    };

    let umi_id = IntToStr::new(umi_seq).into_u64();

    Ok(Some((cell_id, GeneUmiHash(insert_id, umi_id))))
}

fn first_segment_seq<'a>(
    record: &'a FastqRecord,
    hit: &PrimerMatch,
    name: &str,
) -> Option<&'a [u8]> {
    let segment = hit.segments.iter().find(|segment| segment.name == name)?;
    let range = segment.ranges.first()?;

    if range.start > range.end || range.end > record.seq.len() {
        return None;
    }

    Some(&record.seq[range.start..range.end])
}

fn write_output(path: &Path, insert_id: u64, cells: &HashMap<u64, CellHitSet>) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create output {}", path.display()))?;

    let mut out = BufWriter::new(file);

    writeln!(out, "cell_id\tinsert_id\tumi_count\tread_count")?;

    let mut rows: Vec<_> = cells.iter().collect();
    rows.sort_by_key(|(cell_id, _)| **cell_id);

    for (cell_id, hits) in rows {
        writeln!(
            out,
            "{cell_id}\t{insert_id}\t{}\t{}",
            hits.molecules.len(),
            hits.reads
        )?;
    }

    Ok(())
}

fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;

    let is_gz = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

    if is_gz {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn read_fastq_record(reader: &mut dyn BufRead) -> Result<Option<FastqRecord>> {
    let mut id = String::new();

    if reader.read_line(&mut id)? == 0 {
        return Ok(None);
    }

    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();

    reader.read_line(&mut seq)?;
    reader.read_line(&mut plus)?;
    reader.read_line(&mut qual)?;

    let id = id.trim_end();

    if !id.starts_with('@') {
        anyhow::bail!("invalid FASTQ record: id line does not start with @: {id}");
    }

    if !plus.starts_with('+') {
        anyhow::bail!("invalid FASTQ record for {id}: plus line does not start with +");
    }

    let id = id.trim_start_matches('@').to_string();
    let seq = seq.trim_end().as_bytes().to_vec();

    // Existing FastqRecord stores raw Phred scores, not ASCII FASTQ qualities.
    let qual: Vec<u8> = qual
        .trim_end()
        .as_bytes()
        .iter()
        .map(|q| q.saturating_sub(33))
        .collect();

    Ok(Some(FastqRecord::new(id, &seq, &qual)))
}
