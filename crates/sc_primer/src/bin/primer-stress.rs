use clap::Parser;
use fast_tag_mapper::{BuiltinTagSet, FastTagMapper};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use int_to_dna::IntToDna;
use onehot_dna::OneHotSequence;
use mapping_info::MappingInfo;
use sc_primer::{Chemistry, Grammar, PrimerDetector};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

#[derive(Debug, Clone, Copy, clap::ValueEnum)]
enum FeatureSet {
    BdSampleMouse,
    BdSampleHuman,
}

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Stress-test primer detection and cumulative downstream barcode/UMI work on real FASTQ reads"
)]
struct Cli {
    #[arg(long)]
    r1: PathBuf,

    /// R2 is required for molecule_identity(), because the hard UMI key includes R2 bases.
    #[arg(long)]
    r2: Option<PathBuf>,

    #[arg(long, value_enum, num_args = 1.., default_values_t = [Chemistry::BdV2_384])]
    chemistry: Vec<Chemistry>,

    /// Custom primer/read structure grammar. Overrides --chemistry.
    #[arg(long)]
    primer_structure: Option<String>,

    #[arg(long, default_value_t = 250_000)]
    max_reads: usize,

    #[arg(long, default_value_t = 3)]
    iterations: usize,

    /// Print the first N read pairs for which primer detection fails.
    #[arg(long, default_value_t = 0)]
    failed_reads: usize,

    #[arg(long, default_value_t = false)]
    no_forward_only: bool,

    /// Built-in feature set used for the FastTagMapper stress stage.
    #[arg(long, value_enum, default_value_t = FeatureSet::BdSampleMouse)]
    feature_set: FeatureSet,

    /// Minimum unique 8-mer votes required by FastTagMapper.
    #[arg(long, default_value_t = 4)]
    feature_min_hits: u32,

    /// Also benchmark the serial FASTQ output path after normalization.
    #[arg(long, default_value_t = false)]
    io_stress: bool,

    /// Gzip level used by the compressed FASTQ output stress stage.
    #[arg(long, default_value_t = 6)]
    io_gzip_level: u32,
}

#[derive(Clone)]
struct ReadPair {
    r1_id: Vec<u8>,
    r1_seq: Vec<u8>,
    r1_qual: Vec<u8>,
    r2_id: Vec<u8>,
    r2_seq: Vec<u8>,
    r2_qual: Vec<u8>,
}

fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(path).map_err(|e| format!("failed to open '{}': {e}", path.display()))?;
    let gz = path
        .extension()
        .and_then(|x| x.to_str())
        .is_some_and(|x| x.eq_ignore_ascii_case("gz"));
    let reader: Box<dyn Read> = if gz {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::with_capacity(1024 * 1024, reader)))
}

fn trim_newline(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

fn next_fastq<R: BufRead>(
    reader: &mut R,
    record_no: usize,
) -> Result<Option<(Vec<u8>, Vec<u8>, Vec<u8>)>, String> {
    let mut header = Vec::new();
    let mut seq = Vec::new();
    let mut plus = Vec::new();
    let mut qual = Vec::new();
    if reader
        .read_until(b'\n', &mut header)
        .map_err(|e| e.to_string())?
        == 0
    {
        return Ok(None);
    }
    if reader
        .read_until(b'\n', &mut seq)
        .map_err(|e| e.to_string())?
        == 0
        || reader
            .read_until(b'\n', &mut plus)
            .map_err(|e| e.to_string())?
            == 0
        || reader
            .read_until(b'\n', &mut qual)
            .map_err(|e| e.to_string())?
            == 0
    {
        return Err(format!("truncated FASTQ at record {}", record_no + 1));
    }
    trim_newline(&mut header);
    trim_newline(&mut seq);
    trim_newline(&mut qual);
    if seq.len() != qual.len() {
        return Err(format!(
            "sequence/quality length mismatch at record {}",
            record_no + 1
        ));
    }
    seq.make_ascii_uppercase();
    if header.first() == Some(&b'@') {
        header.remove(0);
    }
    Ok(Some((header, seq, qual)))
}

fn load_reads(r1: &Path, r2: Option<&Path>, max_reads: usize) -> Result<Vec<ReadPair>, String> {
    let mut r1_reader = open_fastq(r1)?;
    let mut r2_reader = match r2 {
        Some(path) => Some(open_fastq(path)?),
        None => None,
    };
    let mut reads = Vec::with_capacity(if max_reads == 0 { 250_000 } else { max_reads });
    while max_reads == 0 || reads.len() < max_reads {
        let Some((r1_id, r1_seq, r1_qual)) = next_fastq(&mut r1_reader, reads.len())? else {
            break;
        };
        let (r2_id, r2_seq, r2_qual) = if let Some(reader) = r2_reader.as_mut() {
            let Some((id, seq, qual)) = next_fastq(reader, reads.len())? else {
                return Err("R2 ended before R1".to_string());
            };
            (id, seq, qual)
        } else {
            (Vec::new(), Vec::new(), Vec::new())
        };
        reads.push(ReadPair {
            r1_id,
            r1_seq,
            r1_qual,
            r2_id,
            r2_seq,
            r2_qual,
        });
    }
    Ok(reads)
}

fn report(
    label: &str,
    n: usize,
    called: usize,
    best: Duration,
    total: Duration,
    iterations: usize,
) {
    let n = n as f64;
    let best_s = best.as_secs_f64();
    let mean_s = total.as_secs_f64() / iterations as f64;
    eprintln!("{label}");
    eprintln!(
        "  calls: {called}/{} ({:.2}%)",
        n as usize,
        100.0 * called as f64 / n
    );
    eprintln!(
        "  best: {:.3}s = {:.0} reads/s = {:.3} us/read",
        best_s,
        n / best_s,
        best_s * 1e6 / n
    );
    eprintln!(
        "  mean: {:.3}s = {:.0} reads/s = {:.3} us/read",
        mean_s,
        n / mean_s,
        mean_s * 1e6 / n
    );
}

fn bench<F>(label: &str, reads: &[ReadPair], iterations: usize, mut f: F) -> Result<(), String>
where
    F: FnMut(&ReadPair) -> Result<bool, String>,
{
    let mut best = Duration::MAX;
    let mut total = Duration::ZERO;
    let mut called = 0usize;
    for iteration in 0..iterations {
        let start = Instant::now();
        let mut this_called = 0usize;
        for read in reads {
            if f(read)? {
                this_called += 1;
            }
        }
        let elapsed = start.elapsed();
        best = best.min(elapsed);
        total += elapsed;
        if iteration == 0 {
            called = this_called;
        }
        if called != this_called {
            return Err(format!("non-deterministic call count in {label}"));
        }
        std::hint::black_box(this_called);
    }
    report(label, reads.len(), called, best, total, iterations);
    Ok(())
}

fn write_fastq_record<W: Write>(writer: &mut W, read: &ReadPair) -> Result<(), String> {
    writer.write_all(b"@").map_err(|e| e.to_string())?;
    writer.write_all(&read.r2_id).map_err(|e| e.to_string())?;
    writer.write_all(b"\n").map_err(|e| e.to_string())?;
    writer.write_all(&read.r2_seq).map_err(|e| e.to_string())?;
    writer.write_all(b"\n+\n").map_err(|e| e.to_string())?;
    writer.write_all(&read.r2_qual).map_err(|e| e.to_string())?;
    writer.write_all(b"\n").map_err(|e| e.to_string())?;
    Ok(())
}

fn bench_io<F>(
    label: &str,
    n_input: usize,
    n_output: usize,
    iterations: usize,
    mut f: F,
) -> Result<(), String>
where
    F: FnMut() -> Result<u64, String>,
{
    let mut best = Duration::MAX;
    let mut total = Duration::ZERO;
    let mut expected = None;
    for _ in 0..iterations {
        let start = Instant::now();
        let bytes = f()?;
        let elapsed = start.elapsed();
        best = best.min(elapsed);
        total += elapsed;
        if let Some(want) = expected {
            if bytes != want {
                return Err(format!("non-deterministic byte count in {label}"));
            }
        } else {
            expected = Some(bytes);
        }
        std::hint::black_box(bytes);
    }
    let best_s = best.as_secs_f64();
    let mean_s = total.as_secs_f64() / iterations as f64;
    let bytes = expected.unwrap_or(0);
    eprintln!("{label}");
    eprintln!("  output: {n_output}/{n_input} reads; {bytes} uncompressed FASTQ bytes/pass");
    eprintln!(
        "  best: {:.3}s = {:.0} input reads/s = {:.0} output reads/s",
        best_s,
        n_input as f64 / best_s,
        n_output as f64 / best_s
    );
    eprintln!(
        "  mean: {:.3}s = {:.0} input reads/s = {:.0} output reads/s",
        mean_s,
        n_input as f64 / mean_s,
        n_output as f64 / mean_s
    );
    Ok(())
}


fn trim_r2_primer_readthrough(
    r1_seq: &[u8], r2_seq: &mut Vec<u8>, r2_qual: &mut Vec<u8>, primer_match: &sc_primer::PrimerMatch,
) -> bool {
    use sc_primer::Orientation;
    if primer_match.orientation != Orientation::Forward || primer_match.primer_start >= primer_match.insert_start || primer_match.insert_start > r1_seq.len() { return false; }
    let r1_primer = &r1_seq[primer_match.primer_start..primer_match.insert_start];
    if r1_primer.len() < 2 || r2_seq.len() < 2 { return false; }
    let expected_seq = PrimerDetector::reverse_complement(r1_primer);
    let expected = OneHotSequence::from_iupac_bytes(&expected_seq);
    let observed = OneHotSequence::from_iupac_bytes(r2_seq);
    let lookup = expected.mask_at(0).unwrap();
    let external: Vec<_> = (1..expected.len()).map(|i| expected.mask_at(i).unwrap()).collect();
    let mut hit = observed.find_with_external_after(lookup, &external);
    if hit.is_none() {
        let max_proof = expected.len().min(r2_seq.len());
        for proof in (2..=max_proof).rev() {
            let start = r2_seq.len() - proof;
            if observed.find_next_with_external_after(lookup, &external[..proof - 1], start) == Some(start) { hit = Some(start); break; }
        }
    }
    let Some(start) = hit else { return false; };
    r2_seq.truncate(start); r2_qual.truncate(start); true
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    if cli.iterations == 0 {
        return Err("--iterations must be >= 1".to_string());
    }
    eprintln!(
        "loading {}{} ...",
        cli.r1.display(),
        cli.r2
            .as_ref()
            .map(|p| format!(" + {}", p.display()))
            .unwrap_or_default()
    );
    let start = Instant::now();
    let reads = load_reads(&cli.r1, cli.r2.as_deref(), cli.max_reads)?;
    let elapsed = start.elapsed();
    if reads.is_empty() {
        return Err("FASTQ contained no reads".to_string());
    }
    eprintln!(
        "loaded {} pairs in {:.3}s ({:.0} pairs/s; I/O + decompression excluded from benchmarks)",
        reads.len(),
        elapsed.as_secs_f64(),
        reads.len() as f64 / elapsed.as_secs_f64()
    );
    eprintln!(
        "chemistry: {}",
        cli.chemistry
            .iter()
            .map(|c| c.name())
            .collect::<Vec<_>>()
            .join(", ")
    );
    if let Some(structure) = cli.primer_structure.as_deref() {
        eprintln!("primer structure override: {structure}");
    }
    eprintln!("iterations: {}", cli.iterations);

    let detector = if let Some(structure) = cli.primer_structure.as_deref() {
        PrimerDetector::from_grammar(Grammar::parse("custom", structure)?)
    } else {
        PrimerDetector::from_chemistries(cli.chemistry.iter().copied())
    }
    .map_err(|e| e.to_string())?;

    if cli.failed_reads > 0 {
        let mut failed = 0usize;
        eprintln!("failed primer detections (first {}):", cli.failed_reads);
        for (idx, r) in reads.iter().enumerate() {
            if detector
                .detect_first(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?
                .is_some()
            {
                continue;
            }
            failed += 1;
            eprintln!("--- failed read {} (input record {}) ---", failed, idx + 1);
            eprintln!("R1 id: {}", String::from_utf8_lossy(&r.r1_id));
            eprintln!("R1: {}", String::from_utf8_lossy(&r.r1_seq));
            if cli.r2.is_some() {
                eprintln!("R2 id: {}", String::from_utf8_lossy(&r.r2_id));
                eprintln!("R2: {}", String::from_utf8_lossy(&r.r2_seq));
            }
            if failed == cli.failed_reads {
                break;
            }
        }
        eprintln!("printed {failed} failed primer detections");
    }

    let mut feature_mapper = FastTagMapper::new();
    let _ = match cli.feature_set {
        FeatureSet::BdSampleMouse => feature_mapper.add_builtin(BuiltinTagSet::Mouse),
        FeatureSet::BdSampleHuman => feature_mapper.add_builtin(BuiltinTagSet::Human),
    };

    let feature_mapper = feature_mapper.with_min_hits(cli.feature_min_hits);
    eprintln!(
        "feature mapper: {:?}, {} features, min_hits={}",
        cli.feature_set,
        feature_mapper.feature_count(),
        feature_mapper.min_hits()
    );

    bench("1. detect_first", &reads, cli.iterations, |r| {
        Ok(detector
            .detect_first(&r.r1_seq, &r.r1_qual)
            .map_err(|e| e.to_string())?
            .is_some())
    })?;

    bench(
        "2. detect + get_cell + get_umi",
        &reads,
        cli.iterations,
        |r| {
            let Some(m) = detector
                .detect_first(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?
            else {
                return Ok(false);
            };
            if detector.grammar_for_match(&m).is_unbarcoded() {
                return Ok(true);
            }
            let cell = m
                .get_cell(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?;
            let umi = m
                .get_umi(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?;
            std::hint::black_box((
                cell.seq.len(),
                cell.qual.len(),
                umi.seq.len(),
                umi.qual.len(),
            ));
            Ok(true)
        },
    )?;

    bench(
        "3. detect + slices + normalized-cell clone",
        &reads,
        cli.iterations,
        |r| {
            let Some(m) = detector
                .detect_first(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?
            else {
                return Ok(false);
            };
            if detector.grammar_for_match(&m).is_unbarcoded() {
                return Ok(true);
            }
            let cell = m
                .get_cell(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?;
            let umi = m
                .get_umi(&r.r1_seq, &r.r1_qual)
                .map_err(|e| e.to_string())?;
            let normalized_cell = m.cell_seq.clone().unwrap_or_else(|| cell.seq.clone());
            std::hint::black_box((normalized_cell, umi));
            Ok(true)
        },
    )?;

    if cli.r2.is_some() {
        bench(
            "3b. detect + current R2 primer read-through trimming",
            &reads,
            cli.iterations,
            |r| {
                let Some(m) = detector
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                else {
                    return Ok(false);
                };
                let mut r2_seq = r.r2_seq.clone();
                let mut r2_qual = r.r2_qual.clone();
                let trimmed = trim_r2_primer_readthrough(
                    &r.r1_seq,
                    &mut r2_seq,
                    &mut r2_qual,
                    &m,
                );
                std::hint::black_box((trimmed, r2_seq.len()));
                Ok(true)
            },
        )?;

        bench(
            "4. detect + slices + clone + molecule_identity",
            &reads,
            cli.iterations,
            |r| {
                let Some(m) = detector
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                else {
                    return Ok(false);
                };
                if detector.grammar_for_match(&m).is_unbarcoded() {
                    let Some(id) = detector
                        .grammar_for_match(&m)
                        .molecule_identity_if_exact(None, None, &r.r1_seq, &r.r2_seq)
                        .map_err(|e| e.to_string())?
                    else {
                        return Ok(false);
                    };
                    std::hint::black_box(id);
                    return Ok(true);
                }
                let cell = m
                    .get_cell(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?;
                let umi = m
                    .get_umi(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?;
                let normalized_cell = m.cell_seq.clone().unwrap_or_else(|| cell.seq.clone());
                let Some(id) = detector
                    .grammar_for_match(&m)
                    .molecule_identity_if_exact(Some(&normalized_cell), Some(&umi.seq), &r.r1_seq, &r.r2_seq)
                    .map_err(|e| e.to_string())?
                else {
                    return Ok(false);
                };
                std::hint::black_box(id);
                Ok(true)
            },
        )?;

        bench(
            "5. current identity path + duplicate IntToDna cell/UMI encoding",
            &reads,
            cli.iterations,
            |r| {
                let Some(m) = detector
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                else {
                    return Ok(false);
                };
                if detector.grammar_for_match(&m).is_unbarcoded() {
                    let Some(id) = detector
                        .grammar_for_match(&m)
                        .molecule_identity_if_exact(None, None, &r.r1_seq, &r.r2_seq)
                        .map_err(|e| e.to_string())?
                    else {
                        return Ok(false);
                    };
                    std::hint::black_box(id);
                    return Ok(true);
                }
                let cell = m
                    .get_cell(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?;
                let umi = m
                    .get_umi(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?;
                let normalized_cell = m.cell_seq.clone().unwrap_or_else(|| cell.seq.clone());
                let Some(identity) = detector
                    .grammar_for_match(&m)
                    .molecule_identity_if_exact(Some(&normalized_cell), Some(&umi.seq), &r.r1_seq, &r.r2_seq)
                    .map_err(|e| e.to_string())?
                else {
                    return Ok(false);
                };
                let cell_id = IntToDna::new(&normalized_cell).into_u64();
                let umi_id = IntToDna::new(&umi.seq).into_u64();
                std::hint::black_box((identity, cell_id, umi_id));
                Ok(true)
            },
        )?;
    } else {
        eprintln!("4/5 skipped: provide --r2 to benchmark molecule_identity and duplicate cell/UMI encoding");
    }

    {
        let mut best = Duration::MAX;
        let mut total = Duration::ZERO;
        let mut hits = 0usize;
        for iteration in 0..cli.iterations {
            let mut stats = MappingInfo::new(None, 0.0, reads.len());
            let start = Instant::now();
            let mut this_hits = 0usize;
            for r in &reads {
                let hit = feature_mapper.map_feature_id(&r.r1_seq, &mut stats);
                std::hint::black_box(hit);
                if hit.is_some() {
                    this_hits += 1;
                }
            }
            let elapsed = start.elapsed();
            best = best.min(elapsed);
            total += elapsed;
            if iteration == 0 {
                hits = this_hits;
            }
            if hits != this_hits {
                return Err("non-deterministic hit count in FastTagMapper".to_string());
            }
        }
        report(
            "6a. FastTagMapper only on R1",
            reads.len(),
            hits,
            best,
            total,
            cli.iterations,
        );
    }

    // Production-shaped cumulative benchmark, except that R1 is intentionally
    // fed to FastTagMapper. We are measuring mapper cost, not biological calls.
    // MappingInfo is kept for the full pass, matching normal usage.
    {
        let mut best = Duration::MAX;
        let mut total = Duration::ZERO;
        let mut called = 0usize;
        for iteration in 0..cli.iterations {
            let mut stats = MappingInfo::new(None, 0.0, reads.len());
            let start = Instant::now();
            let mut this_called = 0usize;
            for r in &reads {
                let Some(m) = detector
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                else {
                    continue;
                };
                if detector.grammar_for_match(&m).is_unbarcoded() {
                    if cli.r2.is_some() {
                        let Some(identity) = detector
                            .grammar_for_match(&m)
                            .molecule_identity_if_exact(None, None, &r.r1_seq, &r.r2_seq)
                            .map_err(|e| e.to_string())?
                        else {
                            continue;
                        };
                        std::hint::black_box(identity);
                    }
                } else {
                    let cell = m
                        .get_cell(&r.r1_seq, &r.r1_qual)
                        .map_err(|e| e.to_string())?;
                    let umi = m
                        .get_umi(&r.r1_seq, &r.r1_qual)
                        .map_err(|e| e.to_string())?;
                    let normalized_cell = m.cell_seq.clone().unwrap_or_else(|| cell.seq.clone());
                    if cli.r2.is_some() {
                        let Some(identity) = detector
                            .grammar_for_match(&m)
                            .molecule_identity_if_exact(
                                Some(&normalized_cell),
                                Some(&umi.seq),
                                &r.r1_seq,
                                &r.r2_seq,
                            )
                            .map_err(|e| e.to_string())?
                        else {
                            continue;
                        };
                        let cell_id = IntToDna::new(&normalized_cell).into_u64();
                        let umi_id = IntToDna::new(&umi.seq).into_u64();
                        std::hint::black_box((identity, cell_id, umi_id));
                    }
                }
                let hit = feature_mapper.map_feature_id(&r.r1_seq, &mut stats);
                std::hint::black_box(hit);
                this_called += 1;
            }
            let elapsed = start.elapsed();
            best = best.min(elapsed);
            total += elapsed;
            if iteration == 0 {
                called = this_called;
            }
            if called != this_called {
                return Err("non-deterministic call count in stage 6b".to_string());
            }
        }
        report(
            "6b. current identity path + FastTagMapper(R1)",
            reads.len(),
            called,
            best,
            total,
            cli.iterations,
        );
    }


    if cli.r2.is_some() {
        let mut best = Duration::MAX;
        let mut total = Duration::ZERO;
        let mut called = 0usize;
        let mut trimmed = 0usize;
        let mut feature_hits = 0usize;
        for iteration in 0..cli.iterations {
            let mut stats = MappingInfo::new(None, 0.0, reads.len());
            let start = Instant::now();
            let mut this_called = 0usize;
            let mut this_trimmed = 0usize;
            let mut this_feature_hits = 0usize;
            for r in &reads {
                let Some(m) = detector
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                else {
                    continue;
                };

                if detector.grammar_for_match(&m).is_unbarcoded() {
                    let Some(identity) = detector
                        .grammar_for_match(&m)
                        .molecule_identity_if_exact(None, None, &r.r1_seq, &r.r2_seq)
                        .map_err(|e| e.to_string())?
                    else {
                        continue;
                    };
                    std::hint::black_box(identity);
                } else {
                    let cell = m
                        .get_cell(&r.r1_seq, &r.r1_qual)
                        .map_err(|e| e.to_string())?;
                    let umi = m
                        .get_umi(&r.r1_seq, &r.r1_qual)
                        .map_err(|e| e.to_string())?;
                    let normalized_cell = m.cell_seq.clone().unwrap_or_else(|| cell.seq.clone());
                    let Some(identity) = detector
                        .grammar_for_match(&m)
                        .molecule_identity_if_exact(
                            Some(&normalized_cell),
                            Some(&umi.seq),
                            &r.r1_seq,
                            &r.r2_seq,
                        )
                        .map_err(|e| e.to_string())?
                    else {
                        continue;
                    };
                    let cell_id = IntToDna::new(&normalized_cell).into_u64();
                    let umi_id = IntToDna::new(&umi.seq).into_u64();
                    std::hint::black_box((identity, cell_id, umi_id));
                }

                // Match IlluminaNormalizer ordering: feature-tag reads leave before
                // genomic R2 cleanup, while genomic candidates pay the read-through scan.
                if feature_mapper.map_feature_id(&r.r2_seq, &mut stats).is_some() {
                    this_feature_hits += 1;
                    this_called += 1;
                    continue;
                }

                let mut r2_seq = r.r2_seq.clone();
                let mut r2_qual = r.r2_qual.clone();
                if trim_r2_primer_readthrough(&r.r1_seq, &mut r2_seq, &mut r2_qual, &m) {
                    this_trimmed += 1;
                }
                std::hint::black_box((r2_seq.len(), r2_qual.len()));
                this_called += 1;
            }
            let elapsed = start.elapsed();
            best = best.min(elapsed);
            total += elapsed;
            if iteration == 0 {
                called = this_called;
                trimmed = this_trimmed;
                feature_hits = this_feature_hits;
            }
            if called != this_called || trimmed != this_trimmed || feature_hits != this_feature_hits {
                return Err("non-deterministic counts in stage 6c".to_string());
            }
        }
        report(
            "6c. production-shaped identity + FastTagMapper(R2) + R2 read-through",
            reads.len(),
            called,
            best,
            total,
            cli.iterations,
        );
        eprintln!("  feature hits: {feature_hits}; R2 primer read-through trims: {trimmed}");
    } else {
        eprintln!("6c skipped: provide --r2 for production-shaped R2 processing");
    }

    if cli.io_stress {
        if cli.r2.is_none() {
            eprintln!("7a/7b skipped: --io-stress requires --r2");
        } else {
            // Build the output membership once, outside the timed I/O passes.
            // This approximates the serial R2 output population after primer
            // resolution and feature-tag removal. Deduplication is deliberately
            // not included here: these stages isolate serialization/compression.
            let mut output_indices = Vec::new();
            let mut stats = MappingInfo::new(None, 0.0, reads.len());
            for (idx, r) in reads.iter().enumerate() {
                let Some(_m) = detector
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                else {
                    continue;
                };
                if feature_mapper
                    .map_feature_id(&r.r2_seq, &mut stats)
                    .is_none()
                {
                    output_indices.push(idx);
                }
            }
            eprintln!(
                "I/O stress population: {} R2 reads from {} input pairs ({:.2}%)",
                output_indices.len(),
                reads.len(),
                100.0 * output_indices.len() as f64 / reads.len() as f64
            );

            bench_io(
                "7a. serial FASTQ serialization -> sink (no compression)",
                reads.len(),
                output_indices.len(),
                cli.iterations,
                || {
                    let mut sink = std::io::sink();
                    let mut bytes = 0u64;
                    for &idx in &output_indices {
                        let r = &reads[idx];
                        write_fastq_record(&mut sink, r)?;
                        bytes += (r.r2_id.len() + r.r2_seq.len() + r.r2_qual.len() + 7) as u64;
                    }
                    Ok(bytes)
                },
            )?;

            bench_io(
                "7b. serial FASTQ serialization + gzip -> sink",
                reads.len(),
                output_indices.len(),
                cli.iterations,
                || {
                    let sink = std::io::sink();
                    let mut gz = GzEncoder::new(sink, Compression::new(cli.io_gzip_level));
                    let mut bytes = 0u64;
                    for &idx in &output_indices {
                        let r = &reads[idx];
                        write_fastq_record(&mut gz, r)?;
                        bytes += (r.r2_id.len() + r.r2_seq.len() + r.r2_qual.len() + 7) as u64;
                    }
                    gz.finish().map_err(|e| e.to_string())?;
                    Ok(bytes)
                },
            )?;
        }
    }

    if !cli.no_forward_only {
        let forward = if let Some(structure) = cli.primer_structure.as_deref() {
            PrimerDetector::from_grammar(Grammar::parse("custom", structure)?)
                .map_err(|e| e.to_string())?
        } else {
            PrimerDetector::from_chemistries(cli.chemistry.iter().copied())
                .map_err(|e| e.to_string())?
        }
        .with_reverse_complement_detection(false);
        bench(
            "reference: forward-only detect_first",
            &reads,
            cli.iterations,
            |r| {
                Ok(forward
                    .detect_first(&r.r1_seq, &r.r1_qual)
                    .map_err(|e| e.to_string())?
                    .is_some())
            },
        )?;
    }
    Ok(())
}
