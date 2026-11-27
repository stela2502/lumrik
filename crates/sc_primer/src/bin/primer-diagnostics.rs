use lumrik_clap::{Parser, ReadableFile};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use sc_primer::{BdCellVersion, Chemistry, Grammar, Orientation, PrimerDetector, RhapsodyWhitelist};
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const LINKER1: &[u8; 4] = b"GTGA";
const LINKER2: &[u8; 4] = b"GACA";

#[derive(Debug, Parser)]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(
    author,
    version,
    about = "Inspect single-cell primer recovery and diagnose failed barcode calls from R1 FASTQ"
)]
struct Cli {
    /// R1 FASTQ or FASTQ.GZ file.
    #[arg(long)]
    r1: ReadableFile,

    /// Write resolved R1 reads with cell IDs to a new FASTQ. A .gz suffix enables gzip.
    #[arg(long)]
    out_fastq: Option<PathBuf>,

    /// Write rejected R1 reads with the failure cause in the read ID. A .gz suffix enables gzip.
    #[arg(long)]
    failed_fastq: Option<PathBuf>,

    /// Primer chemistry to diagnose.
    #[arg(long, value_enum, default_value_t = Chemistry::BdV2_384)]
    chemistry: Chemistry,

    /// Custom sc_primer structure, e.g. 'FIXED:ACGT+CELL:8+UMI:6'.
    #[arg(long, conflicts_with = "chemistry")]
    grammar: Option<String>,

    /// Stop after this many reads. Zero means scan the complete file.
    #[arg(long, default_value_t = 0)]
    max_reads: u64,

    /// First BD SEARCH shift to inspect.
    #[arg(long, default_value_t = 0)]
    shift_start: usize,

    /// Last BD SEARCH shift to inspect, inclusive.
    #[arg(long, default_value_t = 4)]
    shift_end: usize,

    /// Run expensive second-pass forensic profiling for reads rejected by the
    /// production primer detector. Off by default: the normal diagnostic path
    /// performs exactly one production detection pass per read.
    #[arg(long, default_value_t = false)]
    deep_failures: bool,
}

fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>, String> {
    let file =
        File::open(path).map_err(|e| format!("failed to open R1 '{}': {e}", path.display()))?;

    let reader: Box<dyn Read> = if is_gzip_path(path) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    Ok(Box::new(BufReader::new(reader)))
}

struct FastqRecord {
    header: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl FastqRecord {
    fn uppercase_sequence(&self) -> Cow<'_, [u8]> {
        if self.seq.iter().any(u8::is_ascii_lowercase) {
            Cow::Owned(self.seq.to_ascii_uppercase())
        } else {
            Cow::Borrowed(&self.seq)
        }
    }

    fn write_with_cell_id(&self, out: &mut dyn Write, cell_id: u64) -> std::io::Result<()> {
        self.write_annotated(out, format_args!("#cell_id:{cell_id}"))
    }

    fn write_failure(&self, out: &mut dyn Write, reason: &str) -> std::io::Result<()> {
        let label: String = reason.chars().map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                ch
            } else {
                '_'
            }
        }).collect();
        self.write_annotated(out, format_args!("#failure:{label}"))
    }

    fn write_annotated(
        &self,
        out: &mut dyn Write,
        suffix: std::fmt::Arguments<'_>,
    ) -> std::io::Result<()> {
        let id_end = self.header
            .iter()
            .position(u8::is_ascii_whitespace)
            .unwrap_or(self.header.len());
        out.write_all(&self.header[..id_end])?;
        out.write_fmt(suffix)?;
        out.write_all(&self.header[id_end..])?;
        out.write_all(b"\n")?;
        out.write_all(&self.seq)?;
        out.write_all(b"\n+\n")?;
        out.write_all(&self.qual)?;
        out.write_all(b"\n")
    }
}

fn export_read(
    out: &mut dyn Write,
    record: &FastqRecord,
    cell_id: u64,
    record_no: u64,
) -> Result<(), String> {
    record
        .write_with_cell_id(out, cell_id)
        .map_err(|e| format!("failed to write FASTQ record {record_no}: {e}"))
}

fn read_fastq_record<R: BufRead>(
    reader: &mut R,
    record_no: u64,
) -> Result<Option<FastqRecord>, String> {
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

    Ok(Some(FastqRecord { header, seq, qual }))
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

fn bd_version(chemistry: Chemistry) -> Option<BdCellVersion> {
    match chemistry {
        Chemistry::BdV1 => Some(BdCellVersion::V1),
        Chemistry::BdV2_96 => Some(BdCellVersion::V2_96),
        Chemistry::BdV2_384 => Some(BdCellVersion::V2_384),
        Chemistry::BdV2_384Vdj => Some(BdCellVersion::V2_384Vdj),
        _ => None,
    }
}

fn run_bd(
    cli: &Cli,
    version: BdCellVersion,
    mut out: Option<&mut dyn Write>,
    mut rejected: Option<&mut dyn Write>,
) -> Result<(), String> {
    let mut reader = open_fastq(&cli.r1)?;
    let whitelist = RhapsodyWhitelist::builtin(version);

    let mut total = 0u64;
    let mut called = 0u64;
    let mut linker_counts = HashMap::<[u8; 8], u64>::new();
    let mut failed_profiles = HashMap::<(u32, u32, u32, usize), u64>::new();
    let mut failed_without_profile = 0u64;
    let mut failed_linker_distances = HashMap::<usize, u64>::new();
    let mut failed_linker_signatures = HashMap::<([u8; 8], usize, usize), u64>::new();
    let mut failed_without_linker_profile = 0u64;
    let mut failed_linker_cell_profiles = HashMap::<(usize, u32, u32, u32), u64>::new();
    let mut failed_linker_examples =
        HashMap::<usize, Vec<(u64, usize, [u8; 8], u32, u32, u32, Vec<u8>)>>::new();
    let mut forward_calls = 0u64;
    let mut reverse_calls = 0u64;

    loop {
        if cli.max_reads != 0 && total >= cli.max_reads {
            break;
        }

        let Some(record) = read_fastq_record(&mut reader, total)? else {
            break;
        };
        total += 1;
        let seq = record.uppercase_sequence();
        let qual = &record.qual;

        // BD cell cassettes are positional. SEARCH supplies only the small
        // tolerated cassette shift; diagnostics must not scan transcript
        // sequence for linker-like motifs.
        let mut call = whitelist.call(&seq, qual, 0, cli.shift_start, cli.shift_end);
        let mut reverse = false;

        if call.is_none() {
            let rc_seq = PrimerDetector::reverse_complement(&seq);
            let mut rc_qual = qual.to_vec();
            rc_qual.reverse();
            call = whitelist.call(&rc_seq, &rc_qual, 0, cli.shift_start, cli.shift_end);
            reverse = call.is_some();
        }

        if let Some(call) = call {
            if let Some(writer) = out.as_deref_mut() {
                export_read(writer, &record, call.cell_id, total)?;
            }
            called += 1;
            if reverse {
                reverse_calls += 1;
            } else {
                forward_calls += 1;
            }

            if let Some(bd) = call.diagnostics() {
                *linker_counts.entry(bd.linker_signature).or_default() += 1;
            }
            continue;
        }

        if let Some(writer) = rejected.as_deref_mut() {
            let reason = whitelist.explain_call_failure(&seq, 0, cli.shift_start, cli.shift_end);
            record.write_failure(writer, &reason)
                .map_err(|e| format!("failed to write rejected FASTQ record {total}: {e}"))?;
        }

        if cli.deep_failures {
            if let Some(profile) =
                whitelist.best_mismatch_profile(&seq, 0, cli.shift_start, cli.shift_end)
            {
                *failed_profiles
                    .entry((
                        profile.c1_mismatches,
                        profile.c2_mismatches,
                        profile.c3_mismatches,
                        profile.shift,
                    ))
                    .or_default() += 1;
            } else {
                failed_without_profile += 1;
            }

            let mut best_linker = None::<([u8; 8], usize, usize)>;
            for shift in cli.shift_start..=cli.shift_end {
                let Some(coords) = whitelist.coords(shift) else {
                    continue;
                };
                if seq.len() < coords.c3.0 {
                    continue;
                }
                let linker1 = &seq[coords.c1.1..coords.c2.0];
                let linker2 = &seq[coords.c2.1..coords.c3.0];
                if linker1.len() != 4 || linker2.len() != 4 {
                    continue;
                }

                let mut signature = [0u8; 8];
                signature[..4].copy_from_slice(linker1);
                signature[4..].copy_from_slice(linker2);
                let distance = hamming(linker1, LINKER1) + hamming(linker2, LINKER2);

                if best_linker.is_none_or(|(_, current_distance, current_shift)| {
                    distance < current_distance
                        || (distance == current_distance && shift < current_shift)
                }) {
                    best_linker = Some((signature, distance, shift));
                }
            }

            if let Some((signature, distance, shift)) = best_linker {
                *failed_linker_distances.entry(distance).or_default() += 1;
                *failed_linker_signatures
                    .entry((signature, distance, shift))
                    .or_default() += 1;

                // Evaluate cell-barcode distance at the SAME shift chosen by the
                // linker, so the two diagnostics describe one candidate cassette.
                if let Some(profile) = whitelist.best_mismatch_profile(&seq, 0, shift, shift) {
                    *failed_linker_cell_profiles
                        .entry((
                            distance,
                            profile.c1_mismatches,
                            profile.c2_mismatches,
                            profile.c3_mismatches,
                        ))
                        .or_default() += 1;

                    if matches!(distance, 4 | 5) {
                        let examples = failed_linker_examples.entry(distance).or_default();
                        if examples.len() < 10 {
                            examples.push((
                                total,
                                shift,
                                signature,
                                profile.c1_mismatches,
                                profile.c2_mismatches,
                                profile.c3_mismatches,
                                seq.as_ref().to_vec(),
                            ));
                        }
                    }
                }
            } else {
                failed_without_linker_profile += 1;
            }
        }
    }

    let failed = total.saturating_sub(called);

    println!("# summary");
    println!("metric\tvalue");
    println!("chemistry\t{}", cli.chemistry.name());
    println!("reads_scanned\t{total}");
    println!("reads_resolved\t{called}");
    println!("reads_failed\t{failed}");
    println!("reads_forward\t{forward_calls}");
    println!("reads_reverse_complement\t{reverse_calls}");
    if total != 0 {
        println!("resolved_fraction\t{:.8}", called as f64 / total as f64);
    }
    if cli.deep_failures {
        println!("failed_without_mismatch_profile\t{failed_without_profile}");
        println!("failed_without_linker_profile\t{failed_without_linker_profile}");
    }

    if !linker_counts.is_empty() {
        let mut ranked = linker_counts.into_iter().collect::<Vec<_>>();
        ranked.sort_unstable_by(|(sig_a, count_a), (sig_b, count_b)| {
            count_b.cmp(count_a).then_with(|| sig_a.cmp(sig_b))
        });

        println!();
        println!("# linker_signatures_resolved_reads");
        println!("signature\tlinker1\tlinker2\tcount\tfraction_resolved\tlinker1_mismatches\tlinker2_mismatches\ttotal_mismatches");
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
    }

    if !failed_linker_distances.is_empty() {
        let mut distances = failed_linker_distances.into_iter().collect::<Vec<_>>();
        distances.sort_unstable_by_key(|(distance, _)| *distance);

        println!();
        println!("# failed_best_linker_distance");
        println!("linker_mismatches\tcount\tfraction_failed");
        for (distance, count) in distances {
            let fraction = if failed == 0 {
                0.0
            } else {
                count as f64 / failed as f64
            };
            println!("{distance}\t{count}\t{fraction:.8}");
        }

        let mut signatures = failed_linker_signatures.into_iter().collect::<Vec<_>>();
        signatures.sort_unstable_by(
            |((sig_a, dist_a, shift_a), count_a), ((sig_b, dist_b, shift_b), count_b)| {
                count_b
                    .cmp(count_a)
                    .then_with(|| dist_a.cmp(dist_b))
                    .then_with(|| shift_a.cmp(shift_b))
                    .then_with(|| sig_a.cmp(sig_b))
            },
        );

        println!();
        println!("# failed_best_linker_signatures");
        println!(
            "signature\tlinker1\tlinker2\tlinker_mismatches\tbest_shift\tcount\tfraction_failed"
        );
        for ((signature, distance, shift), count) in signatures {
            let fraction = if failed == 0 {
                0.0
            } else {
                count as f64 / failed as f64
            };
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.8}",
                String::from_utf8_lossy(&signature),
                String::from_utf8_lossy(&signature[..4]),
                String::from_utf8_lossy(&signature[4..]),
                distance,
                shift,
                count,
                fraction,
            );
        }
    }

    if !failed_linker_cell_profiles.is_empty() {
        let mut rows = failed_linker_cell_profiles.into_iter().collect::<Vec<_>>();
        rows.sort_unstable_by(|((al, a1, a2, a3), ac), ((bl, b1, b2, b3), bc)| {
            al.cmp(bl)
                .then_with(|| (a1 + a2 + a3).cmp(&(b1 + b2 + b3)))
                .then_with(|| bc.cmp(ac))
                .then_with(|| a1.cmp(b1))
                .then_with(|| a2.cmp(b2))
                .then_with(|| a3.cmp(b3))
        });
        let mut totals = HashMap::<usize, u64>::new();
        for ((linker_distance, _, _, _), count) in &rows {
            *totals.entry(*linker_distance).or_default() += *count;
        }
        println!();
        println!("# failed_cell_mismatches_by_linker_distance");
        println!("linker_mismatches\tc1_mismatches\tc2_mismatches\tc3_mismatches\ttotal_cell_mismatches\tcount\tfraction_with_linker_distance");
        for ((linker_distance, c1, c2, c3), count) in rows {
            let denominator = totals.get(&linker_distance).copied().unwrap_or(0);
            let fraction = if denominator == 0 {
                0.0
            } else {
                count as f64 / denominator as f64
            };
            println!(
                "{linker_distance}\t{c1}\t{c2}\t{c3}\t{}\t{count}\t{fraction:.8}",
                c1 + c2 + c3
            );
        }
    }

    if !failed_linker_examples.is_empty() {
        println!();
        println!("# failed_linker_distance_4_5_first_10_reads");
        println!("linker_mismatches\tread_number\tbest_shift\tsignature\tc1_mismatches\tc2_mismatches\tc3_mismatches\ttotal_cell_mismatches\tr1_sequence");
        for linker_distance in [4usize, 5usize] {
            if let Some(examples) = failed_linker_examples.get(&linker_distance) {
                for (read_number, shift, signature, c1, c2, c3, seq) in examples {
                    println!(
                        "{linker_distance}\t{read_number}\t{shift}\t{}\t{c1}\t{c2}\t{c3}\t{}\t{}",
                        String::from_utf8_lossy(signature),
                        c1 + c2 + c3,
                        String::from_utf8_lossy(seq)
                    );
                }
            }
        }
    }

    if !failed_profiles.is_empty() {
        let mut profiles = failed_profiles.into_iter().collect::<Vec<_>>();
        profiles.sort_unstable_by(
            |((a1, a2, a3, ashift), acount), ((b1, b2, b3, bshift), bcount)| {
                let atotal = a1 + a2 + a3;
                let btotal = b1 + b2 + b3;
                atotal
                    .cmp(&btotal)
                    .then_with(|| bcount.cmp(acount))
                    .then_with(|| a1.cmp(b1))
                    .then_with(|| a2.cmp(b2))
                    .then_with(|| a3.cmp(b3))
                    .then_with(|| ashift.cmp(bshift))
            },
        );

        let mut total_distribution = HashMap::<u32, u64>::new();
        for ((c1, c2, c3, _), count) in &profiles {
            *total_distribution.entry(c1 + c2 + c3).or_default() += *count;
        }
        let mut totals = total_distribution.into_iter().collect::<Vec<_>>();
        totals.sort_unstable_by_key(|(distance, _)| *distance);

        println!();
        println!("# failed_total_barcode_mismatches");
        println!("total_mismatches\tcount\tfraction_failed");
        for (distance, count) in totals {
            let fraction = if failed == 0 {
                0.0
            } else {
                count as f64 / failed as f64
            };
            println!("{distance}\t{count}\t{fraction:.8}");
        }

        println!();
        println!("# failed_best_block_profile");
        println!("c1_mismatches\tc2_mismatches\tc3_mismatches\ttotal_mismatches\tbest_shift\tcount\tfraction_failed");
        for ((c1, c2, c3, shift), count) in profiles {
            let fraction = if failed == 0 {
                0.0
            } else {
                count as f64 / failed as f64
            };
            println!(
                "{c1}\t{c2}\t{c3}\t{}\t{shift}\t{count}\t{fraction:.8}",
                c1 + c2 + c3
            );
        }
    }

    eprintln!("R1: {}", cli.r1.display());
    eprintln!("chemistry: {}", cli.chemistry.name());
    eprintln!("reads scanned: {total}");
    eprintln!("reads resolved: {called}");
    if total != 0 {
        eprintln!(
            "resolved fraction: {:.2}%",
            called as f64 * 100.0 / total as f64
        );
    }

    Ok(())
}

fn run_generic(
    cli: &Cli,
    detector: PrimerDetector,
    mut out: Option<&mut dyn Write>,
    mut rejected: Option<&mut dyn Write>,
) -> Result<(), String> {
    let mut reader = open_fastq(&cli.r1)?;
    let chemistry_name = &detector.grammar().name;
    let mut exported = 0u64;
    let mut total = 0u64;
    let mut called = 0u64;
    let mut failures = HashMap::<String, u64>::new();

    loop {
        if cli.max_reads != 0 && total >= cli.max_reads {
            break;
        }
        let Some(record) = read_fastq_record(&mut reader, total)? else {
            break;
        };
        total += 1;
        let seq = record.uppercase_sequence();
        let qual = &record.qual;

        match detector
            .detect_first(&seq, qual)
            .map_err(|e| e.to_string())?
        {
            Some(hit) => {
                if let Some(writer) = out.as_deref_mut() {
                    let cell_id = hit.bd_cell_id.or_else(|| {
                        let cell = if let Some(corrected) = &hit.cell_seq {
                            Cow::Borrowed(corrected.as_slice())
                        } else {
                            let observed = hit.get_cell(&seq, qual).ok()?.seq;
                            Cow::Owned(match hit.orientation {
                                Orientation::Forward => observed,
                                Orientation::ReverseComplement => {
                                    PrimerDetector::reverse_complement(&observed)
                                }
                            })
                        };
                        Grammar::encode_cell_id(&cell).ok()
                    });
                    if let Some(cell_id) = cell_id {
                        export_read(writer, &record, cell_id, total)?;
                        exported += 1;
                    }
                }
                called += 1;
            }
            None => {
                let reason = detector
                    .explain_failure(&seq, qual)
                    .unwrap_or_else(|e| format!("diagnostic error: {e}"));
                if let Some(writer) = rejected.as_deref_mut() {
                    let header_reason = if cli.grammar.is_some() || reason == "no complete primer match" {
                        detector.explain_grammar_failure(&seq, qual)?
                    } else {
                        reason.clone()
                    };
                    record.write_failure(writer, &header_reason)
                        .map_err(|e| format!("failed to write rejected FASTQ record {total}: {e}"))?;
                }
                *failures.entry(reason).or_default() += 1;
            }
        }
    }

    println!("# summary");
    println!("metric\tvalue");
    println!("chemistry\t{chemistry_name}");
    println!("reads_scanned\t{total}");
    println!("reads_resolved\t{called}");
    println!("reads_failed\t{}", total.saturating_sub(called));
    if total != 0 {
        println!("resolved_fraction\t{:.8}", called as f64 / total as f64);
    }

    if out.is_some() {
        println!("reads_exported\t{exported}");
        println!("resolved_without_cell_id\t{}", called - exported);
    }

    if !failures.is_empty() {
        let mut ranked = failures.into_iter().collect::<Vec<_>>();
        ranked.sort_unstable_by(|(reason_a, count_a), (reason_b, count_b)| {
            count_b.cmp(count_a).then_with(|| reason_a.cmp(reason_b))
        });
        println!();
        println!("# failure_reasons");
        println!("reason\tcount\tfraction_failed");
        let failed = total.saturating_sub(called);
        for (reason, count) in ranked {
            let fraction = if failed == 0 {
                0.0
            } else {
                count as f64 / failed as f64
            };
            println!("{reason}\t{count}\t{fraction:.8}");
        }
    }

    eprintln!("R1: {}", cli.r1.display());
    eprintln!("chemistry: {chemistry_name}");
    eprintln!("reads scanned: {total}");
    eprintln!("reads resolved: {called}");
    if total != 0 {
        eprintln!(
            "resolved fraction: {:.2}%",
            called as f64 * 100.0 / total as f64
        );
    }
    Ok(())
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension().is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

enum FastqOutput {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl FastqOutput {
    fn create(path: &Path) -> Result<Self, String> {
        let file = File::create_new(path)
            .map_err(|e| format!("failed to create FASTQ '{}': {e}", path.display()))?;
        let writer = BufWriter::new(file);
        Ok(if is_gzip_path(path) {
            Self::Gzip(GzEncoder::new(writer, Compression::default()))
        } else {
            Self::Plain(writer)
        })
    }

    fn writer(&mut self) -> &mut dyn Write {
        match self {
            Self::Plain(writer) => writer,
            Self::Gzip(writer) => writer,
        }
    }

    fn finish(self) -> std::io::Result<()> {
        match self {
            Self::Plain(mut writer) => writer.flush(),
            Self::Gzip(writer) => writer.finish()?.flush(),
        }
    }
}

fn main() -> Result<(), String> {
    let cli = Cli::parse();
    if cli.shift_start > cli.shift_end {
        return Err(format!(
            "--shift-start ({}) must be <= --shift-end ({})",
            cli.shift_start, cli.shift_end
        ));
    }

    let custom_grammar = cli.grammar.as_deref()
        .map(|structure| Grammar::parse("custom", structure))
        .transpose()?;
    let version = if custom_grammar.is_some() { None } else { bd_version(cli.chemistry) };
    let detector = match (custom_grammar, version) {
        (Some(grammar), _) => Some(PrimerDetector::from_grammar(grammar)?),
        (None, Some(BdCellVersion::V2_96 | BdCellVersion::V2_384 | BdCellVersion::V2_384Vdj)) => None,
        (None, _) => Some(PrimerDetector::from_chemistry(cli.chemistry)?),
    };
    if cli.out_fastq.is_some() {
        if let Some(detector) = &detector {
            let cell_len = detector.grammar().cell_len();
            if cell_len == 0 || cell_len > 32 {
                return Err("--out-fastq requires a cell barcode of 1 to 32 bases".to_string());
            }
        }
    }
    let mut out = cli.out_fastq.as_deref().map(FastqOutput::create).transpose()?;
    let mut rejected = cli.failed_fastq.as_deref().map(FastqOutput::create).transpose()?;

    let accepted_writer = out.as_mut().map(FastqOutput::writer);
    let rejected_writer = rejected.as_mut().map(FastqOutput::writer);
    if let Some(detector) = detector {
        run_generic(&cli, detector, accepted_writer, rejected_writer)?;
    } else {
        run_bd(
            &cli,
            version.expect("BD V2 diagnostics require a BD version"),
            accepted_writer,
            rejected_writer,
        )?;
    }
    for (writer, path) in [(out, &cli.out_fastq), (rejected, &cli.failed_fastq)] {
        if let (Some(writer), Some(path)) = (writer, path) {
            writer.finish()
                .map_err(|e| format!("failed to finish FASTQ '{}': {e}", path.display()))?;
        }
    }
    Ok(())
}
