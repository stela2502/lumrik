use clap::Parser;
use flate2::read::MultiGzDecoder;
use sc_primer::{BdCellVersion, RhapsodyWhitelist};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

const LINKER1: &[u8; 4] = b"GTGA";
const LINKER2: &[u8; 4] = b"GACA";

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Survey observed BD Rhapsody v2 linker sequences from R1 FASTQ reads"
)]
struct Cli {
    /// R1 FASTQ or FASTQ.GZ file.
    #[arg(long)]
    r1: PathBuf,

    /// Stop after this many reads. Zero means scan the complete file.
    #[arg(long, default_value_t = 0)]
    max_reads: u64,
}

fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>, String> {
    let file =
        File::open(path).map_err(|e| format!("failed to open R1 '{}': {e}", path.display()))?;

    let is_gz = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

    let reader: Box<dyn Read> = if is_gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    Ok(Box::new(BufReader::new(reader)))
}

fn read_fastq_record<R: BufRead>(
    reader: &mut R,
    record_no: u64,
) -> Result<Option<(Vec<u8>, Vec<u8>)>, String> {
    let mut header = Vec::new();
    let mut seq = Vec::new();
    let mut plus = Vec::new();
    let mut qual = Vec::new();

    if reader
        .read_until(b'\n', &mut header)
        .map_err(|e| format!("failed reading FASTQ record {} header: {e}", record_no + 1))?
        == 0
    {
        return Ok(None);
    }

    for (name, line) in [
        ("sequence", &mut seq),
        ("plus", &mut plus),
        ("quality", &mut qual),
    ] {
        if reader
            .read_until(b'\n', line)
            .map_err(|e| format!("failed reading FASTQ record {} {name}: {e}", record_no + 1))?
            == 0
        {
            return Err(format!(
                "truncated FASTQ at record {} while reading {name}",
                record_no + 1
            ));
        }
    }

    trim_newline(&mut header);
    trim_newline(&mut seq);
    trim_newline(&mut plus);
    trim_newline(&mut qual);

    if !header.starts_with(b"@") {
        return Err(format!(
            "FASTQ record {} header does not start with '@'",
            record_no + 1
        ));
    }
    if !plus.starts_with(b"+") {
        return Err(format!(
            "FASTQ record {} separator does not start with '+'",
            record_no + 1
        ));
    }
    if seq.len() != qual.len() {
        return Err(format!(
            "FASTQ record {} sequence/quality lengths differ: {} vs {}",
            record_no + 1,
            seq.len(),
            qual.len()
        ));
    }

    seq.make_ascii_uppercase();
    Ok(Some((seq, qual)))
}

fn trim_newline(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

fn hamming(observed: &[u8], expected: &[u8]) -> usize {
    observed
        .iter()
        .zip(expected.iter())
        .filter(|(a, b)| a != b)
        .count()
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    let mut reader = open_fastq(&cli.r1)?;
    let whitelist = RhapsodyWhitelist::builtin(BdCellVersion::V2_384);

    let mut total = 0u64;
    let mut called = 0u64;
    let mut counts = HashMap::<[u8; 8], u64>::new();

    loop {
        if cli.max_reads != 0 && total >= cli.max_reads {
            break;
        }

        let Some((seq, qual)) = read_fastq_record(&mut reader, total)? else {
            break;
        };
        total += 1;

        let Some(call) = whitelist.call(&seq, &qual, 0, 0, 4) else {
            continue;
        };

        let linker1 = &seq[call.c1.1..call.c2.0];
        let linker2 = &seq[call.c2.1..call.c3.0];
        if linker1.len() != 4 || linker2.len() != 4 {
            continue;
        }

        let mut signature = [0u8; 8];
        signature[..4].copy_from_slice(linker1);
        signature[4..].copy_from_slice(linker2);
        *counts.entry(signature).or_default() += 1;
        called += 1;
    }

    let mut ranked = counts.into_iter().collect::<Vec<_>>();
    ranked.sort_unstable_by(|(sig_a, count_a), (sig_b, count_b)| {
        count_b.cmp(count_a).then_with(|| sig_a.cmp(sig_b))
    });

    println!("signature\tlinker1\tlinker2\tcount\tfraction\tlinker1_mismatches\tlinker2_mismatches\ttotal_mismatches");
    for (signature, count) in ranked {
        let d1 = hamming(&signature[..4], LINKER1);
        let d2 = hamming(&signature[4..], LINKER2);
        let fraction = if called == 0 {
            0.0
        } else {
            count as f64 / called as f64
        };

        println!(
            "{}\t{}\t{}\t{}\t{:.8}\t{}\t{}\t{}",
            String::from_utf8_lossy(&signature),
            String::from_utf8_lossy(&signature[..4]),
            String::from_utf8_lossy(&signature[4..]),
            count,
            fraction,
            d1,
            d2,
            d1 + d2,
        );
    }

    eprintln!("R1: {}", cli.r1.display());
    eprintln!("reads scanned: {total}");
    eprintln!("BD v2.384 cells resolved by fuzzy whitelist: {called}");
    if total != 0 {
        eprintln!(
            "resolved fraction: {:.2}%",
            called as f64 * 100.0 / total as f64
        );
    }

    Ok(())
}
