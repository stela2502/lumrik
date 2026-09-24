use anyhow::{Context, Result};
use clap::Parser;
use fast_tag_mapper::{AlignmentStrand, FastTagMapper, FeatureEntry};
use flate2::read::MultiGzDecoder;
use rust_htslib::bam::{self, Read};
use sc_primer::{PrimerCli, PrimerDetector, PrimerMatch};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use bam_tide::fastq::record::FastqRecord;

const TILE_BASES: usize = 48;
const TILE_STEP: usize = 24;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Find reference sequences in FASTQ or BAM reads and retain cell/UMI identity"
)]
struct Cli {
    /// FASTA containing sequences to search for. Long entries are tiled internally.
    #[arg(long)]
    fasta: PathBuf,

    /// Output TSV.
    #[arg(short, long)]
    out: PathBuf,

    /// FASTQ/FASTQ.GZ input. Multiple files are processed consecutively.
    #[arg(long, num_args = 1.., conflicts_with = "bam")]
    fastq: Vec<PathBuf>,

    /// Optional paired R2 FASTQ/FASTQ.GZ. Both mates are searched; cell/UMI
    /// identity is recovered from R1.
    #[arg(long, num_args = 1.., requires = "fastq")]
    r2_fastq: Vec<PathBuf>,

    /// BAM input. Records may be mapped or completely unmapped; raw record
    /// sequence is searched. Existing CB/UB tags are used when present.
    #[arg(long, num_args = 1.., conflicts_with = "fastq")]
    bam: Vec<PathBuf>,

    /// Single-cell chemistry/read structure used to recover cell and UMI from
    /// read sequence when tags are unavailable.
    #[command(flatten)]
    primer: PrimerCli,

    /// Minimum exact 16-mer locators required before verifying a feature hit.
    #[arg(long, default_value_t = 4)]
    min_hits: u32,

    /// Print progress every N reads/records. Set to 0 to disable.
    #[arg(long, default_value_t = 1_000_000)]
    progress_every: usize,
}

#[derive(Debug, Clone)]
struct ReferenceEntry {
    name: String,
    seq: Vec<u8>,
}

#[derive(Debug, Clone)]
struct Tile {
    reference_index: usize,
    start: usize,
    end: usize,
}

#[derive(Debug, Clone)]
struct Hit {
    reference_index: usize,
    start: usize,
    end: usize,
    read_start: usize,
    read_end: usize,
    strand: AlignmentStrand,
}

#[derive(Debug, Default)]
struct FeatureStats {
    reads: usize,
    max_separable_hits_per_read: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.fastq.is_empty() && cli.bam.is_empty() {
        anyhow::bail!("one of --fastq or --bam is required");
    }
    if !cli.r2_fastq.is_empty() && cli.r2_fastq.len() != cli.fastq.len() {
        anyhow::bail!("--r2-fastq must contain the same number of files as --fastq");
    }

    let references = read_fasta(&cli.fasta)?;
    let (mapper, tiles) = build_mapper(&references, cli.min_hits)?;
    let detector = cli
        .primer
        .detector()
        .map_err(anyhow::Error::msg)
        .context("failed to build primer detector")?;

    let file = File::create(&cli.out)
        .with_context(|| format!("failed to create {}", cli.out.display()))?;
    let mut out = BufWriter::new(file);
    writeln!(out, "cell_id\tumi\tfeature\tfeature_start\tfeature_end\tstrand\tmatched_sequence")?;

    let mut total = 0usize;
    let mut matched_reads = 0usize;
    let mut feature_stats: HashMap<usize, FeatureStats> = HashMap::new();

    for (i, path) in cli.fastq.iter().enumerate() {
        let mut r1 = open_fastq(path)?;
        let mut r2 = if cli.r2_fastq.is_empty() {
            None
        } else {
            Some(open_fastq(&cli.r2_fastq[i])?)
        };
        while let Some(record) = read_fastq_record(&mut r1)? {
            let mate = if let Some(reader) = r2.as_mut() {
                Some(read_fastq_record(reader.as_mut())?.context("R2 FASTQ ended before R1")?)
            } else {
                None
            };
            total += 1;
            let identity = primer_identity(&detector, &record);
            let mut hits = find_hits(&mapper, &tiles, &record.seq);
            if let Some(mate) = mate.as_ref() {
                hits.extend(find_hits(&mapper, &tiles, &mate.seq));
            }
            let hits = collapse_hits(hits);
            if !hits.is_empty() {
                matched_reads += 1;
                write_hits(&mut out, identity.as_ref(), &hits, &references)?;
                update_stats(&mut feature_stats, &hits);
            }
            progress(&cli, total, matched_reads);
        }
        if let Some(reader) = r2.as_mut() {
            if read_fastq_record(reader.as_mut())?.is_some() {
                anyhow::bail!("R2 FASTQ contains more records than R1 FASTQ");
            }
        }
    }

    for path in &cli.bam {
        let mut reader = bam::Reader::from_path(path)
            .with_context(|| format!("failed to open BAM {}", path.display()))?;
        for result in reader.records() {
            let record = result?;
            total += 1;
            let seq = record.seq().as_bytes();
            let identity = bam_identity(&record).or_else(|| primer_identity_bytes(&detector, &seq, &record.qual()));
            let hits = collapse_hits(find_hits(&mapper, &tiles, &seq));
            if !hits.is_empty() {
                matched_reads += 1;
                write_hits(&mut out, identity.as_ref(), &hits, &references)?;
                update_stats(&mut feature_stats, &hits);
            }
            progress(&cli, total, matched_reads);
        }
    }

    out.flush()?;
    eprintln!("done: reads={} reads_with_hits={}", total, matched_reads);
    eprintln!("feature\treads_with_feature\tmax_separable_hits_per_read");
    let mut rows: Vec<_> = feature_stats.into_iter().collect();
    rows.sort_by_key(|(idx, _)| *idx);
    for (idx, stats) in rows {
        eprintln!(
            "{}\t{}\t{}",
            references[idx].name, stats.reads, stats.max_separable_hits_per_read
        );
    }
    Ok(())
}

fn build_mapper(references: &[ReferenceEntry], min_hits: u32) -> Result<(FastTagMapper, Vec<Tile>)> {
    let mut mapper = FastTagMapper::new().with_min_hits(min_hits);
    let mut tiles = Vec::new();
    let mut feature_id = 1u64;
    for (reference_index, reference) in references.iter().enumerate() {
        if reference.seq.len() < 16 {
            anyhow::bail!("FASTA entry '{}' is shorter than 16 bp", reference.name);
        }
        let tile_len = reference.seq.len().min(TILE_BASES);
        let mut start = 0usize;
        loop {
            let end = (start + tile_len).min(reference.seq.len());
            mapper.add_feature(
                &reference.seq[start..end],
                FeatureEntry::new(feature_id, reference.name.clone(), "sequence_query"),
            );
            tiles.push(Tile { reference_index, start, end });
            feature_id += 1;
            if end == reference.seq.len() {
                break;
            }
            start = (start + TILE_STEP).min(reference.seq.len() - tile_len);
        }
    }
    Ok((mapper, tiles))
}

fn find_hits(mapper: &FastTagMapper, tiles: &[Tile], seq: &[u8]) -> Vec<Hit> {
    mapper
        .align_all(seq, None)
        .into_iter()
        .filter_map(|alignment| {
            let tile = tiles.get(alignment.feature_index)?;
            let read_start = alignment.query_start;
            let read_end = read_start.checked_add(tile.end - tile.start)?;
            Some(Hit {
                reference_index: tile.reference_index,
                start: tile.start,
                end: tile.end,
                read_start,
                read_end,
                strand: alignment.strand,
            })
        })
        .collect()
}

/// Adjacent overlapping tiles are evidence for one reference occurrence, not
/// independent copies. Merge them before reporting so long references remain
/// readable and repeated ONT occurrences stay separable.
fn collapse_hits(mut hits: Vec<Hit>) -> Vec<Hit> {
    hits.sort_by_key(|h| (h.reference_index, strand_key(h.strand), h.read_start, h.start));
    let mut out: Vec<Hit> = Vec::new();
    for hit in hits {
        if let Some(last) = out.last_mut() {
            if last.reference_index == hit.reference_index
                && last.strand == hit.strand
                && hit.read_start <= last.read_end
                && hit.start <= last.end
            {
                last.end = last.end.max(hit.end);
                last.read_end = last.read_end.max(hit.read_end);
                continue;
            }
        }
        out.push(hit);
    }
    out
}

fn strand_key(strand: AlignmentStrand) -> u8 {
    match strand {
        AlignmentStrand::Forward => 0,
        AlignmentStrand::Reverse => 1,
    }
}

fn write_hits(
    out: &mut impl Write,
    identity: Option<&(String, String)>,
    hits: &[Hit],
    references: &[ReferenceEntry],
) -> Result<()> {
    let (cell, umi) = identity
        .map(|(cell, umi)| (cell.as_str(), umi.as_str()))
        .unwrap_or((".", "."));
    for hit in hits {
        let reference = &references[hit.reference_index];
        let strand = match hit.strand {
            AlignmentStrand::Forward => "+",
            AlignmentStrand::Reverse => "-",
        };
        let matched = String::from_utf8_lossy(&reference.seq[hit.start..hit.end]);
        writeln!(
            out,
            "{cell}\t{umi}\t{}\t{}\t{}\t{strand}\t{matched}",
            reference.name, hit.start, hit.end
        )?;
    }
    Ok(())
}

fn update_stats(stats: &mut HashMap<usize, FeatureStats>, hits: &[Hit]) {
    let mut per_feature: BTreeMap<usize, usize> = BTreeMap::new();
    for hit in hits {
        *per_feature.entry(hit.reference_index).or_default() += 1;
    }
    for (feature, copies) in per_feature {
        let entry = stats.entry(feature).or_default();
        entry.reads += 1;
        entry.max_separable_hits_per_read = entry.max_separable_hits_per_read.max(copies);
    }
}

fn primer_identity(detector: &PrimerDetector, record: &FastqRecord) -> Option<(String, String)> {
    primer_identity_bytes(detector, &record.seq, &record.qual)
}

fn primer_identity_bytes(detector: &PrimerDetector, seq: &[u8], qual: &[u8]) -> Option<(String, String)> {
    let hit = detector.detect_first(seq, qual).ok()??;
    let cell = hit.cell_seq.as_ref().map(|x| String::from_utf8_lossy(x).into_owned())
        .or_else(|| hit.bd_cell_id.map(|x| x.to_string()))?;
    let umi = first_segment_seq(seq, &hit, "UMI")?;
    Some((cell, String::from_utf8_lossy(umi).into_owned()))
}

fn bam_identity(record: &bam::Record) -> Option<(String, String)> {
    use bam::record::Aux;
    let cell = match record.aux(b"CB").ok()? { Aux::String(x) => x.to_string(), _ => return None };
    let umi = match record.aux(b"UB").ok()? { Aux::String(x) => x.to_string(), _ => return None };
    Some((cell, umi))
}

fn first_segment_seq<'a>(seq: &'a [u8], hit: &PrimerMatch, name: &str) -> Option<&'a [u8]> {
    let segment = hit.segments.iter().find(|segment| segment.name == name)?;
    let range = segment.ranges.first()?;
    seq.get(range.start..range.end)
}

fn read_fasta(path: &Path) -> Result<Vec<ReferenceEntry>> {
    let reader = BufReader::new(File::open(path).with_context(|| format!("failed to open FASTA {}", path.display()))?);
    let mut entries = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() { continue; }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(old) = name.take() {
                entries.push(ReferenceEntry { name: old, seq: std::mem::take(&mut seq) });
            }
            name = Some(header.split_whitespace().next().unwrap_or(header).to_string());
        } else {
            seq.extend(line.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(old) = name {
        entries.push(ReferenceEntry { name: old, seq });
    }
    if entries.is_empty() { anyhow::bail!("FASTA contains no entries"); }
    Ok(entries)
}

fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("failed to open FASTQ {}", path.display()))?;
    if path.extension().is_some_and(|x| x == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn read_fastq_record(reader: &mut dyn BufRead) -> Result<Option<FastqRecord>> {
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 { return Ok(None); }
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();
    reader.read_line(&mut seq)?;
    reader.read_line(&mut plus)?;
    reader.read_line(&mut qual)?;
    if !header.starts_with('@') || !plus.starts_with('+') {
        anyhow::bail!("invalid FASTQ record near {}", header.trim_end());
    }
    let seq = seq.trim_end().as_bytes().to_vec();
    let qual: Vec<u8> = qual.trim_end().as_bytes().iter().map(|q| q.saturating_sub(33)).collect();
    if seq.len() != qual.len() { anyhow::bail!("FASTQ sequence/quality length mismatch"); }
    Ok(Some(FastqRecord::new(header.trim_start_matches('@').trim_end().to_string(), &seq, &qual)))
}

fn progress(cli: &Cli, total: usize, matched: usize) {
    if cli.progress_every > 0 && total % cli.progress_every == 0 {
        eprintln!("processed={} reads_with_hits={}", total, matched);
    }
}
