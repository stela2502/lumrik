use clap::Parser;
use sc_primer::{PrimerCli, PrimerDetector};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Identify primer/barcode/UMI structure in one DNA sequence or a line-delimited sequence file"
)]
struct Cli {
    #[command(flatten)]
    primer: PrimerCli,

    /// DNA sequence to inspect, or a path to a text file containing one sequence per line.
    #[arg(long)]
    seq: String,

    /// Optional FASTQ quality string. Only valid when --seq is a literal sequence.
    #[arg(long)]
    qual: Option<String>,
}

fn sequence_qual(seq: &[u8], qual: Option<&str>) -> Result<Vec<u8>, String> {
    let qual = qual
        .map(str::as_bytes)
        .map(ToOwned::to_owned)
        .unwrap_or_else(|| vec![b'I'; seq.len()]);

    if seq.len() != qual.len() {
        return Err(format!(
            "sequence and quality have different lengths: seq={} qual={}",
            seq.len(),
            qual.len()
        ));
    }

    Ok(qual)
}

fn print_single(detector: &PrimerDetector, seq: &[u8], qual: &[u8]) -> Result<(), String> {
    let attempts = detector.explain_all(seq, qual)?;

    if attempts.is_empty() {
        println!("no primer attempts");
        return Ok(());
    }

    let mut matched = 0usize;

    for attempt in &attempts {
        if !attempt.ok {
            continue;
        }

        matched += 1;
        let prefix = std::str::from_utf8(&seq[..attempt.offset]).unwrap_or("<non-utf8>");
        let cell_seq = attempt.cell_seq.as_deref().unwrap_or("");
        println!(
            "  prefix: {}bp [0..{}] {}\n  cell_seq: {}",
            prefix, attempt.offset, attempt.offset, cell_seq,
        );
        println!(
            "offset: {} orientation: {:?} status: OK reason: {}",
            attempt.offset, attempt.orientation, attempt.reason
        );

        for segment in &attempt.segments {
            let dna = std::str::from_utf8(&seq[segment.range.start..segment.range.end])
                .unwrap_or("<non-utf8>");

            println!(
                "  {}: {}bp [{}..{}] {} | {} | {}",
                segment.name,
                segment.range.end.saturating_sub(segment.range.start),
                segment.range.start,
                segment.range.end,
                dna,
                if segment.ok { "OK" } else { "FAIL" },
                segment.reason
            );
        }
    }

    if matched == 0 {
        println!("summary: no complete primer match");
        println!("reason: {}", detector.explain_failure(seq, qual)?);
    } else {
        println!("summary: {matched} complete primer match(es)\n");
    }

    Ok(())
}

fn print_file(detector: &PrimerDetector, path: &Path) -> Result<(), String> {
    let file = File::open(path)
        .map_err(|e| format!("failed to open sequence file '{}': {e}", path.display()))?;
    let reader = BufReader::new(file);

    let mut total = 0usize;
    let mut valid = 0usize;
    let mut matches = 0usize;
    let mut errors = BTreeMap::<String, usize>::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| {
            format!(
                "failed to read line {} from '{}': {e}",
                line_no + 1,
                path.display()
            )
        })?;
        let seq = line.trim().as_bytes();
        if seq.is_empty() {
            continue;
        }

        total += 1;
        let qual = vec![b'I'; seq.len()];

        match detector.detect_all(seq, &qual) {
            Ok(hits) if !hits.is_empty() => {
                valid += 1;
                matches += hits.len();
            }
            Ok(_) => {
                let reason = detector.explain_failure(seq, &qual)?;
                *errors.entry(reason).or_default() += 1;
            }
            Err(error) => {
                *errors
                    .entry(format!("detector error: {error}"))
                    .or_default() += 1;
            }
        }
    }

    println!("file: {}", path.display());
    println!("sequences: {total}");
    println!("valid: {valid}");
    println!("invalid: {}", total.saturating_sub(valid));
    println!("complete primer matches: {matches}");
    println!("errors:");
    if errors.is_empty() {
        println!("  0\tnone");
    } else {
        for (reason, count) in errors {
            println!("  {count}\t{reason}");
        }
    }

    Ok(())
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    let detector = cli.primer.detector()?;
    let path = Path::new(&cli.seq);

    if path.is_file() {
        if cli.qual.is_some() {
            return Err("--qual cannot be used when --seq points to a file".to_string());
        }
        print_file(&detector, path)
    } else {
        let seq = cli.seq.as_bytes();
        let qual = sequence_qual(seq, cli.qual.as_deref())?;
        print_single(&detector, seq, &qual)
    }
}
