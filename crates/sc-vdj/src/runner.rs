use crate::cellrep::{
    AlignmentGeometry, BamFeatureEvidence, BamFeatureSequenceParts, CellEvidenceVdj, EvidenceId,
    MapperEvidence, SequencePart,
};
use crate::index::VdjIndex;
use crate::recombination::{process_chain_work, ChainWork, Recombination};
use anyhow::{Context, Result};
use int_to_str::IntToStr;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::path::Path;

pub trait BamIdentityResolver {
    fn cell(&self, record: &bam::Record) -> Option<String>;
}
#[derive(Debug, Clone, Copy, Default)]
pub struct NelruneIdentityResolver;
impl BamIdentityResolver for NelruneIdentityResolver {
    fn cell(&self, record: &bam::Record) -> Option<String> {
        if let Ok(Aux::String(s)) = record.aux(b"CB") {
            return Some(s.to_string());
        }
        let fields: Vec<_> = record.qname().split(|b| *b == b'|').collect();
        decode_hex_ascii(fields.get(1).copied()?)
    }
}

#[derive(Debug, Clone)]
pub struct VdjRunnerConfig {
    pub min_sequence_overlap: usize,
}
impl Default for VdjRunnerConfig {
    fn default() -> Self {
        Self {
            min_sequence_overlap: 12,
        }
    }
}

#[derive(Debug)]
pub struct VdjRunner {
    pub index: VdjIndex,
    pub evidence: CellEvidenceVdj,
    pub cell_names: HashMap<u64, String>,
    pub config: VdjRunnerConfig,
    flush_id: u32,
}
impl VdjRunner {
    pub fn new(index: VdjIndex, config: VdjRunnerConfig) -> Self {
        Self {
            index,
            evidence: CellEvidenceVdj::new(),
            cell_names: HashMap::new(),
            config,
            flush_id: 0,
        }
    }
    pub fn read_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &mut self,
        path: P,
        resolver: &R,
    ) -> Result<usize> {
        let mut reader = bam::Reader::from_path(path.as_ref())
            .with_context(|| format!("opening {}", path.as_ref().display()))?;
        let header = reader.header().to_owned();
        let mut n = 0usize;
        let mut entry = 0u32;
        for rec in reader.records() {
            let rec = rec?;
            if rec.is_unmapped() {
                continue;
            }
            let Some(cell) = resolver.cell(&rec) else {
                continue;
            };
            let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
            let tid = rec.tid();
            if tid < 0 {
                continue;
            }
            let chr = String::from_utf8_lossy(header.tid2name(tid as u32)).into_owned();
            let blocks = aligned_blocks(&rec);
            let segment_ids = self.index.overlapping(&chr, &blocks);
            if segment_ids.is_empty() {
                continue;
            }
            let geometry = AlignmentGeometry {
                tid,
                start: blocks.first().map_or(0, |x| x.0),
                end: blocks.last().map_or(0, |x| x.1),
                is_reverse: rec.is_reverse(),
                is_secondary: rec.is_secondary(),
                is_supplementary: rec.is_supplementary(),
                mapq: rec.mapq(),
                ref_blocks: blocks,
            };
            let mappings = segment_ids
                .into_iter()
                .map(|segment_id| MapperEvidence {
                    segment_id,
                    alignment: geometry.clone(),
                })
                .collect();
            let part = SequencePart {
                bases: rec.seq().as_bytes(),
                qualities: rec.qual().to_vec(),
            };
            let mut sequence = BamFeatureSequenceParts::default();
            if rec.is_last_in_template() {
                sequence.r2 = Some(part)
            } else {
                sequence.r1 = Some(part)
            };
            self.evidence.push(
                cell_id,
                BamFeatureEvidence {
                    id: EvidenceId {
                        flush: self.flush_id,
                        entry,
                    },
                    sequence,
                    mappings,
                },
            );
            self.cell_names.entry(cell_id).or_insert(cell);
            entry = entry.wrapping_add(1);
            n += 1;
        }
        self.flush_id = self.flush_id.wrapping_add(1);
        Ok(n)
    }
    pub fn identify(&self) -> Vec<(u64, Vec<Recombination>)> {
        // Materialize independent cell/locus work units first.  The execution
        // line below is intentionally serial today; this Vec is the future
        // Rayon boundary (into_par_iter) without changing assembly internals.
        let work: Vec<_> = self
            .evidence
            .cells()
            .flat_map(|(cell_id, cell)| {
                cell.chains(&self.index).into_iter().map(move |chain| {
                    (
                        cell_id,
                        ChainWork {
                            cell,
                            index: &self.index,
                            chain,
                            min_overlap: self.config.min_sequence_overlap,
                        },
                    )
                })
            })
            .collect();

        let mut by_cell = HashMap::<u64, Vec<Recombination>>::new();
        for (cell_id, chain_work) in work {
            by_cell
                .entry(cell_id)
                .or_default()
                .extend(process_chain_work(chain_work));
        }

        let mut out: Vec<_> = by_cell.into_iter().collect();
        out.sort_by_key(|x| x.0);
        out
    }
}

fn aligned_blocks(record: &bam::Record) -> Vec<(u32, u32)> {
    let mut ref_pos = record.pos().max(0) as u32;
    let mut out = Vec::new();
    let mut block_start = None;
    for c in record.cigar().iter() {
        match *c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                if block_start.is_none() {
                    block_start = Some(ref_pos)
                }
                ref_pos = ref_pos.saturating_add(n)
            }
            Cigar::Del(n) => {
                if block_start.is_none() {
                    block_start = Some(ref_pos)
                }
                ref_pos = ref_pos.saturating_add(n)
            }
            Cigar::RefSkip(n) => {
                if let Some(a) = block_start.take() {
                    if a < ref_pos {
                        out.push((a, ref_pos))
                    }
                }
                ref_pos = ref_pos.saturating_add(n)
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if let Some(a) = block_start {
        if a < ref_pos {
            out.push((a, ref_pos))
        }
    }
    out
}
fn decode_hex_ascii(input: &[u8]) -> Option<String> {
    if input.is_empty() || input.len() % 2 != 0 {
        return None;
    }
    let mut o = Vec::with_capacity(input.len() / 2);
    for x in input.chunks_exact(2) {
        o.push((hex(x[0])? << 4) | hex(x[1])?)
    }
    String::from_utf8(o).ok()
}
fn hex(x: u8) -> Option<u8> {
    match x {
        b'0'..=b'9' => Some(x - b'0'),
        b'a'..=b'f' => Some(x - b'a' + 10),
        b'A'..=b'F' => Some(x - b'A' + 10),
        _ => None,
    }
}
