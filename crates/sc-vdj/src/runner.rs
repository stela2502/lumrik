use crate::cellrep::{
    AlignmentGeometry, BamFeatureEvidence, BamFeatureSequenceParts, CellEvidenceVdj, CellEvidence, EvidenceId,
    MapperEvidence, SequencePart,
};
use crate::index::VdjIndex;
use crate::recombination::{
    process_chain_work, rescue_missing_constants_from_bam,
    rescue_missing_constants_from_bam_with_report, ChainWork, Recombination,
    RecombinationEvidenceRescanReport,
};
use anyhow::{Context, Result};
use int_to_str::IntToStr;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

const EVIDENCE_BATCH_SIZE: usize = 200_000;

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
    threads: usize,
}
impl VdjRunner {
    pub fn new(index: VdjIndex, config: VdjRunnerConfig) -> Self {
        Self {
            index,
            evidence: CellEvidenceVdj::new(),
            cell_names: HashMap::new(),
            config,
            flush_id: 0,
            threads: 1,
        }
    }
    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads.max(1);
    }

    pub fn read_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &mut self,
        path: P,
        resolver: &R,
    ) -> Result<usize> {
        let mut reader = bam::Reader::from_path(path.as_ref())
            .with_context(|| format!("opening {}", path.as_ref().display()))?;
        if self.threads > 1 {
            reader
                .set_threads(self.threads)
                .context("configuring multithreaded BAM decoding")?;
        }
        let header = reader.header().to_owned();
        let mut n = 0usize;
        let mut entry = 0u32;
        let mut batch = Vec::<(u64, BamFeatureEvidence)>::with_capacity(EVIDENCE_BATCH_SIZE);

        // Mapper BAMs emit the records belonging to one physical query together.
        // Keep only the immediately preceding query key so paired/supplementary
        // records share an EvidenceId without retaining a BAM-sized QNAME map.
        // A full batch is flushed only when a new query begins, so one physical
        // fragment is never split merely because it crossed the 20k boundary.
        let mut last_query: Option<(u64, Vec<u8>)> = None;
        let mut current_evidence_id: Option<EvidenceId> = None;

        for rec in reader.records() {
            let rec = rec?;
            if rec.is_unmapped() {
                continue;
            }
            let Some(cell) = resolver.cell(&rec) else {
                continue;
            };
            let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
            let query_key = (cell_id, rec.qname().to_vec());

            if last_query.as_ref() != Some(&query_key) {
                if batch.len() >= EVIDENCE_BATCH_SIZE {
                    let full = std::mem::replace(
                        &mut batch,
                        Vec::with_capacity(EVIDENCE_BATCH_SIZE),
                    );
                    self.evidence.consume_batch(
                        full,
                        &self.index,
                        self.config.min_sequence_overlap,
                    );
                }
                last_query = Some(query_key);
                current_evidence_id = None;
            }

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

            let id = *current_evidence_id.get_or_insert_with(|| {
                let id = EvidenceId {
                    flush: self.flush_id,
                    entry,
                };
                entry = entry.wrapping_add(1);
                id
            });

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

            batch.push((
                cell_id,
                BamFeatureEvidence {
                    id,
                    sequence,
                    mappings,
                },
            ));
            self.cell_names.entry(cell_id).or_insert(cell);
            n += 1;
        }

        if !batch.is_empty() {
            self.evidence.consume_batch(
                batch,
                &self.index,
                self.config.min_sequence_overlap,
            );
        }
        self.flush_id = self.flush_id.wrapping_add(1);
        Ok(n)
    }
    /// Re-scan the retained BAM with cell-specific reconstructed receptor baits.
    ///
    /// This independently rediscovers receptor reads and records direct
    /// receptor-to-constant linkage support. Missing constant calls may be
    /// filled when one constant segment has uniquely stronger direct support.
    pub fn rediscover_receptor_linkage_from_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<usize> {
        rescue_missing_constants_from_bam(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
        )
    }

    pub fn rediscover_receptor_linkage_from_bam_with_report<
        P: AsRef<Path>,
        R: BamIdentityResolver,
    >(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<RecombinationEvidenceRescanReport> {
        rescue_missing_constants_from_bam_with_report(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
        )
    }

    pub fn rescue_missing_constants_from_bam<P: AsRef<Path>, R: BamIdentityResolver>(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<usize> {
        rescue_missing_constants_from_bam(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
        )
    }

    pub fn rescue_missing_constants_from_bam_with_report<
        P: AsRef<Path>,
        R: BamIdentityResolver,
    >(
        &self,
        path: P,
        resolver: &R,
        calls: &mut [(u64, Vec<Recombination>)],
    ) -> Result<RecombinationEvidenceRescanReport> {
        rescue_missing_constants_from_bam_with_report(
            path,
            resolver,
            &self.index,
            calls,
            self.threads,
        )
    }
    pub fn identify(&self) -> Vec<(u64, Vec<Recombination>)> {
        // Cells are the natural production-scale Rayon boundary.  Keep the
        // chains within one cell serial: there are normally many more cells
        // than worker threads, this preserves cell-local cache locality, and
        // it avoids scheduling tiny nested Rayon jobs.
        //
        // `cell` and `&self.index` are shared immutable references.  Neither
        // CellEvidence nor VdjIndex is cloned into a worker.
        let cells: Vec<_> = self.evidence.cells().collect();
        let process_cell = |(cell_id, cell): (u64, &CellEvidence)| {
            let recombinations = cell
                .chains(&self.index)
                .into_iter()
                .flat_map(|chain| {
                    process_chain_work(ChainWork {
                        cell,
                        index: &self.index,
                        chain,
                        min_overlap: self.config.min_sequence_overlap,
                    })
                })
                .collect::<Vec<_>>();

            (cell_id, recombinations)
        };

        // Use a private Rayon pool for real parallel work, but make
        // `--threads 1` genuinely serial rather than falling back to the
        // process-global Rayon pool.
        let mut out = if self.threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .expect("building sc-vdj Rayon pool")
                .install(|| {
                    cells
                        .into_par_iter()
                        .map(process_cell)
                        .collect::<Vec<_>>()
                })
        } else {
            cells.into_iter().map(process_cell).collect::<Vec<_>>()
        };

        // Each worker returns one complete result vector for one cell, so
        // there is no shared mutable result map and no cross-thread merge.
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
