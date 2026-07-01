use clap::Parser;
use sc_primer::PrimerCli;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Identify primer/barcode/UMI structure in one DNA sequence"
)]
struct Cli {
    #[command(flatten)]
    primer: PrimerCli,

    /// DNA sequence to inspect.
    #[arg(long)]
    seq: String,

    /// Optional FASTQ quality string. If omitted, dummy high-quality scores are used.
    #[arg(long)]
    qual: Option<String>,
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();

    let seq = cli.seq.as_bytes();

    let qual: Vec<u8> = cli
        .qual
        .unwrap_or_else(|| "I".repeat(seq.len()))
        .into_bytes();

    if seq.len() != qual.len() {
        return Err(format!(
            "sequence and quality have different lengths: seq={} qual={}",
            seq.len(),
            qual.len()
        ));
    }

    let detector = cli.primer.detector()?;
    let attempts = detector.explain_all(seq, &qual)?;

    if attempts.is_empty() {
        println!("no primer attempts");
        return Ok(());
    }

    let mut matched = 0usize;

    for attempt in &attempts {
        if attempt.ok {
            let prefix = std::str::from_utf8(&seq[..attempt.offset]).unwrap_or("<non-utf8>");
            let cell_seq = match &attempt.cell_seq {
                Some(seq) => seq,
                None => &String::new(),
            };
            println!(
                "  prefix: {}bp [0..{}] {}\n  cell_seq: {}",
                prefix, attempt.offset, attempt.offset, cell_seq,
            );
        }

        if !attempt.ok {
            continue;
        } else {
            matched += 1;
        }
        println!(
            "offset: {} orientation: {:?} status: {} reason: {}",
            attempt.offset,
            attempt.orientation,
            if attempt.ok { "OK" } else { "FAIL" },
            attempt.reason
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

        //println!();
    }

    if matched == 0 {
        println!("summary: no complete primer match");
    } else {
        println!("summary: {matched} complete primer match(es)\n");
    }

    Ok(())
}
