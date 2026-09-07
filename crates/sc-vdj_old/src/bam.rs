use crate::posterior::BamReadEvidence;
use crate::reference::VdjReference;
use crate::cell_evidence::{absorb_provisional, BamFeatureSequenceParts, CellEvidence, CellRecombinationResults, InitialMapperResults, MappedSegment};
use rayon::prelude::*;
use int_to_str::IntToStr;
use scdata::CellHash;
use anyhow::{Context, Result};
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::path::Path;


#[derive(Debug, Clone, Default)]
pub struct BamEvidenceStats {
    pub receptor_chromosomes: usize,
    pub indexed_segment_intervals: usize,
    pub total_records: usize,
    pub receptor_chromosome_records: usize,
    pub segment_overlap_records: usize,
    pub called_cell_records: usize,
    pub routed_evidence_records: usize,
    pub batch_flushes: usize,
    pub compacted_sequences: usize,
    pub redundant_segments: usize,
    pub non_receptor_chromosome_records: usize,
    pub non_segment_overlap_records: usize,
    pub locus_records: [usize; 7],
}

pub const DEFAULT_ROUTED_EVIDENCE_BATCH_SIZE: usize = 1_000_000;

/// Cell-keyed receptor evidence. `scdata::CellHash` owns the u64 cell key and
/// the 256-way bucket partition; the value is the complete provisional/final
/// V(D)J evidence for that cell.
#[derive(Debug)]
pub struct RoutedBamEvidence {
    cells: CellHash<CellEvidence>,
}

impl RoutedBamEvidence {
    pub fn new() -> Self {
        Self {
            cells: CellHash::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.cells.cell_count()
    }

    pub fn is_empty(&self) -> bool {
        self.cells.cells_are_empty()
    }

    pub(crate) fn sequence_count(&self) -> usize {
        self.cells
            .buckets()
            .iter()
            .flat_map(|bucket| bucket.values())
            .flat_map(|cell| cell.final_results.values())
            .map(Vec::len)
            .sum()
    }

    pub(crate) fn ensure_cell(&mut self, cell: String) {
        let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
        self.cells
            .entry_cell(cell_id)
            .or_insert_with(|| CellEvidence {
                cell,
                ..CellEvidence::default()
            });
    }

    pub(crate) fn into_cells(
        self,
    ) -> Vec<(
        String,
        HashMap<crate::types::Chain, Vec<CellRecombinationResults>>,
    )> {
        self.cells
            .into_iter()
            .flat_map(|bucket| bucket.into_values())
            .map(|cell| (cell.cell, cell.final_results))
            .collect()
    }

    /// Collapse all currently provisional BAM feature sequences into the
    /// persistent per-cell recombination representation. The cell buckets are
    /// independent, so this is the only Rayon boundary needed by the collector.
    fn finalize_provisional(&mut self) -> (usize, usize) {
        self.cells
            .buckets_mut()
            .par_iter_mut()
            .map(|bucket| {
                bucket.values_mut().fold((0usize, 0usize), |(p, m), cell| {
                    let (processed, merged) = absorb_provisional(cell);
                    (p + processed, m + merged)
                })
            })
            .reduce(|| (0usize, 0usize), |(pa, ma), (pb, mb)| (pa + pb, ma + mb))
    }
}

#[derive(Debug, Clone, Copy)]
struct SegmentInterval {
    start: u32,
    end: u32,
    segment_index: usize,
    chain: crate::types::Chain,
    kind: crate::types::SegmentKind,
}

/// Small BAM-TID routing table derived directly from the serialized VDJ index.
/// Chromosomes without indexed V/D/J/C segments have an empty interval vector
/// and are rejected before cell/UMI parsing or sequence extraction.
struct BamLocusRouter {
    by_tid: Vec<Vec<SegmentInterval>>,
}

impl BamLocusRouter {
    fn new(header: &bam::HeaderView, reference: &VdjReference) -> Self {
        let mut tid_by_name = HashMap::<String, usize>::new();
        for tid in 0..header.target_count() {
            tid_by_name.insert(
                String::from_utf8_lossy(header.tid2name(tid)).into_owned(),
                tid as usize,
            );
        }
        let mut by_tid = vec![Vec::new(); header.target_count() as usize];
        for (segment_index, segment) in reference.segments.iter().enumerate() {
            let Some(&tid) = tid_by_name.get(&segment.chr) else {
                continue;
            };
            if segment.end <= segment.start {
                continue;
            }
            by_tid[tid].push(SegmentInterval {
                start: segment.start,
                end: segment.end,
                segment_index,
                chain: segment.chain,
                kind: segment.kind,
            });
        }
        for intervals in &mut by_tid {
            intervals.sort_by_key(|interval| (interval.start, interval.end, chain_index(interval.chain), interval.segment_index));
        }
        Self { by_tid }
    }

    fn has_receptor_chromosome(&self, tid: i32) -> bool {
        tid >= 0
            && self
                .by_tid
                .get(tid as usize)
                .is_some_and(|intervals| !intervals.is_empty())
    }

    fn matching_segments(&self, record: &bam::Record) -> Vec<MappedSegment> {
        self.matching_segments_for_blocks(record.tid(), &aligned_blocks(record))
    }

    fn matching_segments_for_blocks(
        &self,
        tid: i32,
        blocks: &[(u32, u32)],
    ) -> Vec<MappedSegment> {
        if tid < 0 || blocks.is_empty() {
            return Vec::new();
        }
        let Some(intervals) = self.by_tid.get(tid as usize) else {
            return Vec::new();
        };
        if intervals.is_empty() {
            return Vec::new();
        }

        let mut out = Vec::new();
        for &(block_start, block_end) in blocks {
            for interval in intervals {
                if interval.start >= block_end {
                    break;
                }
                if interval.end > block_start
                    && !out.iter().any(|hit: &MappedSegment| hit.segment_index == interval.segment_index)
                {
                    out.push(MappedSegment {
                        segment_index: interval.segment_index,
                        chain: interval.chain,
                        kind: interval.kind,
                    });
                }
            }
        }
        out.sort_by_key(|hit| (chain_index(hit.chain), hit.kind, hit.segment_index));
        out
    }

    fn matching_chains(&self, record: &bam::Record) -> Vec<crate::types::Chain> {
        let mut chains = Vec::new();
        for hit in self.matching_segments(record) {
            if !chains.contains(&hit.chain) {
                chains.push(hit.chain);
            }
        }
        chains
    }

}

fn chain_index(chain: crate::types::Chain) -> usize {
    match chain {
        crate::types::Chain::Igh => 0,
        crate::types::Chain::Igk => 1,
        crate::types::Chain::Igl => 2,
        crate::types::Chain::Tra => 3,
        crate::types::Chain::Trb => 4,
        crate::types::Chain::Trg => 5,
        crate::types::Chain::Trd => 6,
    }
}

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
        let qualities = record.qual().to_vec();
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
            qualities,
            chr,
            ref_start: start,
            ref_end: end,
            ref_blocks,
            mapped_segment_indices: Vec::new(),
            mapq: record.mapq(),
            is_reverse: record.is_reverse(),
            is_secondary: record.is_secondary(),
            is_supplementary: record.is_supplementary(),
        });
    }
    Ok(out)
}

/// Stream a whole-transcriptome BAM exactly once and retain only called-cell
/// alignments whose reference blocks overlap an actual indexed V/D/J/C segment.
/// Evidence is routed immediately into independent cell/locus buckets. No
/// temporary BAMs or shard files are created.
pub fn read_bam_receptor_evidence<P, R, F>(
    path: P,
    resolver: &R,
    keep_cell: F,
    reference: &VdjReference,
) -> Result<(RoutedBamEvidence, BamEvidenceStats)>
where
    P: AsRef<Path>,
    R: BamIdentityResolver,
    F: Fn(&str) -> bool,
{
    read_bam_receptor_evidence_with_progress(path, resolver, keep_cell, reference, |_| {})
}

/// Progress-enabled form of [`read_bam_receptor_evidence`]. The callback is
/// invoked every 50,000 BAM records and once at EOF.
pub fn read_bam_receptor_evidence_with_progress<P, R, F, G>(
    path: P,
    resolver: &R,
    keep_cell: F,
    reference: &VdjReference,
    progress: G,
) -> Result<(RoutedBamEvidence, BamEvidenceStats)>
where
    P: AsRef<Path>,
    R: BamIdentityResolver,
    F: Fn(&str) -> bool,
    G: FnMut(&BamEvidenceStats),
{
    read_bam_receptor_evidence_with_progress_and_batch_size(
        path,
        resolver,
        keep_cell,
        reference,
        DEFAULT_ROUTED_EVIDENCE_BATCH_SIZE,
        progress,
    )
}

/// Same routing pipeline as [`read_bam_receptor_evidence_with_progress`], but
/// with an explicit staging threshold. This is primarily useful for integration
/// tests that need to exercise repeated collect -> compact -> collect cycles
/// without manufacturing a million BAM records.
pub fn read_bam_receptor_evidence_with_progress_and_batch_size<P, R, F, G>(
    path: P,
    resolver: &R,
    keep_cell: F,
    reference: &VdjReference,
    batch_size: usize,
    mut progress: G,
) -> Result<(RoutedBamEvidence, BamEvidenceStats)>
where
    P: AsRef<Path>,
    R: BamIdentityResolver,
    F: Fn(&str) -> bool,
    G: FnMut(&BamEvidenceStats),
{
    if batch_size == 0 {
        anyhow::bail!("candidate batch size must be greater than zero");
    }
    let path = path.as_ref();
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let header = reader.header().to_owned();
    let router = BamLocusRouter::new(&header, reference);
    let mut stats = BamEvidenceStats {
        receptor_chromosomes: router.by_tid.iter().filter(|intervals| !intervals.is_empty()).count(),
        indexed_segment_intervals: router.by_tid.iter().map(|intervals| intervals.len()).sum(),
        ..BamEvidenceStats::default()
    };
    let mut evidence = RoutedBamEvidence::new();
    let mut provisional_count = 0usize;
    let mut next_sequence_id = 0u64;
    let mut provisional_by_qname = HashMap::<Box<[u8]>, (u64, usize)>::new();
    let mut record = bam::Record::new();

    while let Some(result) = reader.read(&mut record) {
        result?;
        stats.total_records += 1;
        if stats.total_records % 50_000 == 0 {
            progress(&stats);
        }
        if record.is_secondary() {
            continue;
        }

        if !router.has_receptor_chromosome(record.tid()) {
            stats.non_receptor_chromosome_records += 1;
            continue;
        }
        stats.receptor_chromosome_records += 1;

        let mapped_segments = router.matching_segments(&record);
        if mapped_segments.is_empty() {
            stats.non_segment_overlap_records += 1;
            continue;
        }
        stats.segment_overlap_records += 1;

        let Some((cell, _umi)) = resolver.resolve(&record) else {
            continue;
        };
        if !keep_cell(&cell) {
            continue;
        }
        stats.called_cell_records += 1;

        let cell_id = IntToStr::new(cell.as_bytes()).into_u64();
        let read2 = record.is_last_in_template();
        let existing_idx = provisional_by_qname
            .get(record.qname())
            .and_then(|&(stored_cell_id, idx)| (stored_cell_id == cell_id).then_some(idx));

        let cell_evidence = evidence
            .cells
            .entry_cell(cell_id)
            .or_insert_with(|| CellEvidence {
                cell,
                ..CellEvidence::default()
            });

        let provisional_idx = if let Some(idx) = existing_idx {
            idx
        } else {
            let idx = cell_evidence.provisional.len();
            let sequence_id = next_sequence_id;
            next_sequence_id = next_sequence_id.wrapping_add(1);
            cell_evidence.provisional.push((
                BamFeatureSequenceParts {
                    sequence_id,
                    ..BamFeatureSequenceParts::default()
                },
                Vec::new(),
            ));
            provisional_by_qname.insert(
                record.qname().to_vec().into_boxed_slice(),
                (cell_id, idx),
            );
            idx
        };

        let (parts, mappings) = &mut cell_evidence.provisional[provisional_idx];
        let sequence_slot = if read2 { &mut parts.r2 } else { &mut parts.r1 };
        if sequence_slot.is_none() {
            *sequence_slot = Some((record.seq().as_bytes(), record.qual().to_vec()));
        }
        mappings.push(InitialMapperResults {
            segments: mapped_segments.clone(),
            read2,
            ref_blocks: aligned_blocks(&record),
            mapq: record.mapq(),
            is_supplementary: record.is_supplementary(),
        });

        let mut chains = Vec::new();
        for segment in &mapped_segments {
            if !chains.contains(&segment.chain) {
                chains.push(segment.chain);
            }
        }
        for chain in chains {
            provisional_count += 1;
            stats.locus_records[chain_index(chain)] += 1;
            stats.routed_evidence_records += 1;
        }

        if provisional_count >= batch_size {
            let (_processed, merged) = evidence.finalize_provisional();
            provisional_count = 0;
            provisional_by_qname.clear();
            stats.batch_flushes += 1;
            stats.redundant_segments += merged;
        }
    }

    if provisional_count > 0 {
        let (_processed, merged) = evidence.finalize_provisional();
        stats.batch_flushes += 1;
        stats.redundant_segments += merged;
    }
    stats.compacted_sequences = evidence.sequence_count();
    progress(&stats);
    Ok((evidence, stats))
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
    fn routes_only_actual_segment_overlaps() {
        let router = BamLocusRouter {
            by_tid: vec![vec![
                SegmentInterval { start: 100, end: 150, segment_index: 0, chain: crate::types::Chain::Igh, kind: crate::types::SegmentKind::V },
                SegmentInterval { start: 300, end: 350, segment_index: 1, chain: crate::types::Chain::Igk, kind: crate::types::SegmentKind::C },
            ]],
        };
        assert_eq!(router.matching_segments_for_blocks(0, &[(120, 130)]).iter().map(|hit| (hit.segment_index, hit.chain, hit.kind)).collect::<Vec<_>>(), vec![(0, crate::types::Chain::Igh, crate::types::SegmentKind::V)]);
        assert!(router.matching_segments_for_blocks(0, &[(200, 250)]).is_empty());
        assert_eq!(router.matching_segments_for_blocks(0, &[(320, 330)]).iter().map(|hit| (hit.segment_index, hit.chain, hit.kind)).collect::<Vec<_>>(), vec![(1, crate::types::Chain::Igk, crate::types::SegmentKind::C)]);
    }

    #[test]
    fn nested_tra_trd_segments_can_route_to_both_loci() {
        let router = BamLocusRouter {
            by_tid: vec![vec![
                SegmentInterval { start: 100, end: 160, segment_index: 0, chain: crate::types::Chain::Tra, kind: crate::types::SegmentKind::C },
                SegmentInterval { start: 130, end: 180, segment_index: 1, chain: crate::types::Chain::Trd, kind: crate::types::SegmentKind::C },
            ]],
        };
        assert_eq!(
            router.matching_segments_for_blocks(0, &[(140, 150)]).iter().map(|hit| (hit.chain, hit.kind)).collect::<Vec<_>>(),
            vec![(crate::types::Chain::Tra, crate::types::SegmentKind::C), (crate::types::Chain::Trd, crate::types::SegmentKind::C)]
        );
    }

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
