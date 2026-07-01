
use clap::Parser;
use rust_htslib::bam::{self, Read, record::Cigar, record::Record};
use std::collections::VecDeque;
use rust_htslib::bam::record::Aux;

/// CLI tool to report the positions of a mutation in a bam read (MD tag)
/// and report that as normalized position or as a histogram of occurances.
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input BAM file
    #[arg(short, long)]
    input: String,

    /// Output histogram of positions
    #[arg(short='H', long)]
    histogram: bool,

    /// Set bin count
    #[arg(short, long, default_value_t=100)]
    bins: usize,
}

fn main() {
    let args = Args::parse();

    let mut bam = bam::Reader::from_path(&args.input).expect("Failed to open BAM file");

    let mut histogram_bins = vec![0u32; args.bins ];

    for r in bam.records() {
        let record = r.expect("Failed to read record");

        if record.is_unmapped() {
            continue;
        }

        let read_len = record.seq().len();
        if read_len == 0 {
            continue;
        }

        let md_tag = match record.aux(b"MD") {
            Ok(Aux::String(md_bytes)) => md_bytes.to_string(),
            _ => continue, // Skip if MD tag missing or not a string
        };

        let mismatches = extract_mismatches(&record, &md_tag);

        for mismatch_read_pos in mismatches {
            let norm = mismatch_read_pos as f64 / read_len as f64;

            if args.histogram {
                let bin = (norm * args.bins as f64 ).floor() as usize;
                if bin < args.bins {
                    histogram_bins[bin] += 1;
                }
            } else {
                println!("{:.4}", norm);
            }
        }
    }

    if args.histogram {
        for (i, count) in histogram_bins.iter().enumerate() {
            println!("{:>3}: {}", i, count);
        }
    }
}

fn extract_mismatches(record: &Record, md_tag: &str) -> Vec<u32> {
    let mut mismatches = Vec::new();
    let cigar = record.cigar();
    
    #[allow(unused_variables)]
    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut ref_to_read = VecDeque::new();

    for c in cigar.iter() {
        match c {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                for _ in 0..*l {
                    ref_to_read.push_back(Some(read_pos));
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(l) => {
                read_pos += *l;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                for _ in 0..*l {
                    ref_to_read.push_back(None);
                    ref_pos += 1;
                }
            }
            Cigar::SoftClip(l) => {
                read_pos += *l;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    let mut chars = md_tag.chars().peekable();
    let mut ref_pos_in_tag = 0;

    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            let mut num = String::new();
            while let Some(&d) = chars.peek() {
                if d.is_ascii_digit() {
                    num.push(d);
                    chars.next();
                } else {
                    break;
                }
            }
            ref_pos_in_tag += num.parse::<usize>().unwrap();
        } else if c == '^' {
            chars.next(); // skip ^
            while let Some(&d) = chars.peek() {
                if d.is_ascii_alphabetic() {
                    ref_pos_in_tag += 1;
                    chars.next();
                } else {
                    break;
                }
            }
        } else if c.is_ascii_alphabetic() {
            if let Some(Some(read_pos)) = ref_to_read.get(ref_pos_in_tag) {
                mismatches.push(*read_pos);
            }
            ref_pos_in_tag += 1;
            chars.next();
        } else {
            chars.next(); // unexpected char
        }
    }

    mismatches
}
