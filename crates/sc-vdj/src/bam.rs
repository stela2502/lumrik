use crate::posterior::BamReadEvidence;
use anyhow::{Context, Result};
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use std::path::Path;

pub trait BamIdentityResolver {
    fn resolve(&self, record: &bam::Record) -> Option<(String, String)>;
}

#[derive(Debug, Clone)]
pub struct QnameIdentityResolver {
    pub separator: u8,
    pub cell_field: usize,
    pub umi_field: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct AuxTagIdentityResolver {
    pub cell_tag: [u8; 2],
    pub umi_tag: [u8; 2],
}

impl Default for AuxTagIdentityResolver {
    fn default() -> Self {
        Self {
            cell_tag: *b"CB",
            umi_tag: *b"UB",
        }
    }
}

impl BamIdentityResolver for AuxTagIdentityResolver {
    fn resolve(&self, record: &bam::Record) -> Option<(String, String)> {
        let cell = match record.aux(&self.cell_tag).ok()? {
            Aux::String(s) => s.to_string(),
            _ => return None,
        };
        let umi = match record.aux(&self.umi_tag).ok()? {
            Aux::String(s) => s.to_string(),
            _ => return None,
        };
        Some((cell, umi))
    }
}

impl BamIdentityResolver for QnameIdentityResolver {
    fn resolve(&self, record: &bam::Record) -> Option<(String, String)> {
        let q = record.qname();
        let fields: Vec<&[u8]> = q.split(|b| *b == self.separator).collect();
        let cell = String::from_utf8(fields.get(self.cell_field)?.to_vec()).ok()?;
        let umi = String::from_utf8(fields.get(self.umi_field)?.to_vec()).ok()?;
        Some((cell, umi))
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct NelruneIdentityResolver;

impl BamIdentityResolver for NelruneIdentityResolver {
    fn resolve(&self, record: &bam::Record) -> Option<(String, String)> {
        if let Some(identity) = AuxTagIdentityResolver::default().resolve(record) {
            return Some(identity);
        }

        let fields: Vec<&[u8]> = record.qname().split(|b| *b == b'|').collect();
        let cell = decode_hex_ascii(*fields.get(1)?)?;
        let umi = decode_hex_ascii(*fields.get(3)?)?;
        Some((cell, umi))
    }
}

fn decode_hex_ascii(input: &[u8]) -> Option<String> {
    if input.is_empty() || input.len() % 2 != 0 {
        return None;
    }
    let mut decoded = Vec::with_capacity(input.len() / 2);
    for pair in input.chunks_exact(2) {
        let hi = hex_nibble(pair[0])?;
        let lo = hex_nibble(pair[1])?;
        decoded.push((hi << 4) | lo);
    }
    String::from_utf8(decoded).ok()
}

fn hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

pub fn read_bam<P: AsRef<Path>, R: BamIdentityResolver>(
    path: P,
    resolver: &R,
) -> Result<Vec<BamReadEvidence>> {
    read_bam_filtered(path, resolver, |_| true)
}

pub fn read_bam_filtered<P, R, F>(
    path: P,
    resolver: &R,
    keep_cell: F,
) -> Result<Vec<BamReadEvidence>>
where
    P: AsRef<Path>,
    R: BamIdentityResolver,
    F: Fn(&str) -> bool,
{
    let path = path.as_ref();
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let header = reader.header().to_owned();
    let mut record = bam::Record::new();
    let mut out = Vec::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let Some((cell, umi)) = resolver.resolve(&record) else {
            continue;
        };
        if !keep_cell(&cell) {
            continue;
        }
        let seq = record.seq().as_bytes();
        let read_name = String::from_utf8_lossy(record.qname()).to_string();
        let tid = record.tid();
        let (chr, start, end, ref_blocks) = if tid >= 0 {
            let chr = String::from_utf8_lossy(header.tid2name(tid as u32)).to_string();
            let s = record.pos().max(0) as u32;
            let blocks = aligned_blocks(&record);
            let e = blocks.last().map(|x| x.1).unwrap_or(s);
            (Some(chr), Some(s), Some(e), blocks)
        } else {
            (None, None, None, Vec::new())
        };
        out.push(BamReadEvidence {
            cell,
            umi,
            read_name,
            sequence: seq,
            chr,
            ref_start: start,
            ref_end: end,
            ref_blocks,
            mapq: record.mapq(),
            is_reverse: record.is_reverse(),
            is_secondary: record.is_secondary(),
            is_supplementary: record.is_supplementary(),
        });
    }
    Ok(out)
}

fn aligned_blocks(record: &bam::Record) -> Vec<(u32, u32)> {
    let mut pos = record.pos().max(0) as u32;
    let mut block_start: Option<u32> = None;
    let mut out = Vec::new();
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                if block_start.is_none() {
                    block_start = Some(pos)
                }
                pos = pos.saturating_add(n);
            }
            Cigar::Del(n) => {
                if block_start.is_none() {
                    block_start = Some(pos)
                }
                pos = pos.saturating_add(n);
            }
            Cigar::RefSkip(n) => {
                if let Some(s) = block_start.take() {
                    if pos > s {
                        out.push((s, pos));
                    }
                }
                pos = pos.saturating_add(n);
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if let Some(s) = block_start {
        if pos > s {
            out.push((s, pos));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_nelrune_mapper_qname_hex_fields() {
        let mut record = bam::Record::new();
        record.set_qname(
            b"read|435454415447544343474347434341544154544143545447544347|2828|434147434743|2828|",
        );
        let (cell, umi) = NelruneIdentityResolver.resolve(&record).expect("identity");
        assert_eq!(cell, "CTTATGTCCGCGCCATATTACTTGTCG");
        assert_eq!(umi, "CAGCGC");
    }
}
