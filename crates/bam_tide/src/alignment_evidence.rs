use anyhow::{Context, Result};
use rust_htslib::bam::{
    self, Read,
    record::{Aux, Cigar},
};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct BamAlignmentEvidence {
    pub query: String,
    pub target: String,
    pub start: u64,
    pub end: u64,
    pub reverse: bool,
    pub cigar: String,
    pub mapq: u8,
    pub edit_distance: Option<u32>,
    pub query_start: u32,
    pub query_end: u32,
    pub query_len: u32,
    pub secondary: bool,
    pub supplementary: bool,
}

pub fn read_alignment_evidence(path: impl AsRef<Path>) -> Result<Vec<BamAlignmentEvidence>> {
    let path = path.as_ref();
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let targets: Vec<String> = reader
        .header()
        .target_names()
        .iter()
        .map(|x| String::from_utf8_lossy(x).into_owned())
        .collect();
    let mut out = Vec::new();
    for rec in reader.records() {
        let rec = rec.with_context(|| format!("reading BAM {}", path.display()))?;
        if rec.is_unmapped() || rec.tid() < 0 {
            continue;
        }
        let query_len = rec.seq_len() as u32;
        let mut left_clip = 0u32;
        let mut right_clip = 0u32;
        if let Some(c) = rec.cigar().iter().next() {
            if let Cigar::SoftClip(n) | Cigar::HardClip(n) = c {
                left_clip = *n;
            }
        }
        if let Some(c) = rec.cigar().iter().last() {
            if let Cigar::SoftClip(n) | Cigar::HardClip(n) = c {
                right_clip = *n;
            }
        }
        let full_query_len = query_len.saturating_add(
            rec.cigar()
                .iter()
                .filter_map(|c| {
                    if let Cigar::HardClip(n) = c {
                        Some(*n)
                    } else {
                        None
                    }
                })
                .sum::<u32>(),
        );
        let edit_distance = match rec.aux(b"NM") {
            Ok(Aux::U8(v)) => Some(v as u32),
            Ok(Aux::U16(v)) => Some(v as u32),
            Ok(Aux::U32(v)) => Some(v),
            Ok(Aux::I8(v)) if v >= 0 => Some(v as u32),
            Ok(Aux::I16(v)) if v >= 0 => Some(v as u32),
            Ok(Aux::I32(v)) if v >= 0 => Some(v as u32),
            _ => None,
        };
        let tid = rec.tid() as usize;
        out.push(BamAlignmentEvidence {
            query: String::from_utf8_lossy(rec.qname()).into_owned(),
            target: targets
                .get(tid)
                .cloned()
                .unwrap_or_else(|| format!("tid:{}", tid)),
            start: rec.pos().max(0) as u64,
            end: rec.cigar().end_pos().max(0) as u64,
            reverse: rec.is_reverse(),
            cigar: rec.cigar().to_string(),
            mapq: rec.mapq(),
            edit_distance,
            query_start: left_clip,
            query_end: full_query_len.saturating_sub(right_clip),
            query_len: full_query_len,
            secondary: rec.is_secondary(),
            supplementary: rec.is_supplementary(),
        });
    }
    Ok(out)
}
