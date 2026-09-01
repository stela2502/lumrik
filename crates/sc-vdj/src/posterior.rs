use crate::align::{best_oriented_local_alignment, best_oriented_local_alignment_with_reverse, OrientedLocalAlignment};
use crate::gex::ExpressionMatrix;
use crate::mapper::VdjMapper;
use crate::reference::VdjReference;
use crate::score::{score_development, DevelopmentEvidence};
use crate::sterile::{SterileAccumulator, SterileProfile};
use crate::types::{Chain, SegmentKind};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct BamReadEvidence {
    pub cell: String,
    pub umi: String,
    pub read_name: String,
    pub sequence: Vec<u8>,
    pub chr: Option<String>,
    pub ref_start: Option<u32>,
    pub ref_end: Option<u32>,
    /// Exact aligned reference blocks; RefSkip introns are not counted as transcriptional coverage.
    pub ref_blocks: Vec<(u32, u32)>,
    pub mapq: u8,
    pub is_reverse: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
}

#[derive(Debug, Clone)]
pub struct GermlineSegmentSupport {
    pub segment_index: usize,
    pub id: String,
    pub kind: SegmentKind,
    pub local_alignment_score: i32,
    pub supporting_umis: usize,
    pub supporting_reads: usize,
    pub locus_fraction: f64,
    pub distance_to_recombination_center: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecombinationStage {
    None,
    Dj,
    Vj,
    Vdj,
}

#[derive(Debug, Clone)]
pub struct ReadSegmentAlignment {
    pub score: i32,
    /// Coordinates on the read in transcript orientation. When `reverse_complement`
    /// is true these coordinates refer to the reverse-complemented BAM sequence.
    pub query_start: usize,
    pub query_end: usize,
    pub reference_start: usize,
    pub reference_end: usize,
    pub reverse_complement: bool,
}

#[derive(Debug, Clone)]
pub struct RearrangementSupportingRead {
    pub umi: String,
    pub read_name: String,
    /// Sequence exactly as stored in BAM. It is not silently reverse-complemented.
    pub sequence: Vec<u8>,
    pub bam_is_reverse: bool,
    pub is_supplementary: bool,
    pub supports_v: bool,
    pub supports_j: bool,
    pub supports_d: bool,
    pub supports_c: bool,
    pub v_alignment: Option<ReadSegmentAlignment>,
    pub j_alignment: Option<ReadSegmentAlignment>,
    pub d_alignment: Option<ReadSegmentAlignment>,
    pub c_alignment: Option<ReadSegmentAlignment>,
}

#[derive(Debug, Clone)]
pub struct RearrangementCall {
    pub chain: Chain,
    pub stage: RecombinationStage,
    pub v: Option<GermlineSegmentSupport>,
    pub d: Option<GermlineSegmentSupport>,
    pub j: Option<GermlineSegmentSupport>,
    pub c: Option<GermlineSegmentSupport>,
    /// UMIs with coherent support for both the selected V and J segments.
    pub total_supporting_umis: usize,
    /// Original BAM reads from those coherent UMIs that support at least one
    /// segment in the final call. Used only for auditable sequence diagnostics.
    pub supporting_reads: Vec<RearrangementSupportingRead>,
    pub notation: String,
}

#[derive(Debug, Clone)]
pub struct CellVdjSummary {
    pub cell: String,
    pub rearrangements: Vec<RearrangementCall>,
    pub sterile: Vec<SterileProfile>,
    pub development: DevelopmentEvidence,
}

#[derive(Debug, Clone)]
pub struct PosteriorConfig {
    pub sterile_bins: usize,
    pub min_seed_hits: u32,
    pub candidate_segments_per_kind: usize,
    pub min_local_score: i32,
    /// Reject C candidates implausibly far from the chain's J recombination center.
    pub max_constant_distance_bp: u64,
    /// Maximum V/J overlap on one oriented read. Large overlap usually means
    /// the two germlines are explaining the same conserved sequence.
    pub max_vj_alignment_overlap: usize,
    /// V evidence must reach this close to the 3' end of the germline V.
    pub max_v_end_distance: usize,
    /// J evidence must begin this close to the 5' end of the germline J.
    pub max_j_start_distance: usize,
    /// Minimum observed sequence overlap when separate V- and J-anchored reads
    /// from the same UMI are used to establish a continuous V->J path.
    pub min_vj_read_overlap: usize,
}
impl Default for PosteriorConfig {
    fn default() -> Self {
        Self {
            sterile_bins: 64,
            min_seed_hits: 2,
            candidate_segments_per_kind: 8,
            min_local_score: 18,
            max_constant_distance_bp: 1_000_000,
            max_vj_alignment_overlap: 12,
            max_v_end_distance: 35,
            max_j_start_distance: 35,
            min_vj_read_overlap: 20,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct BestReadMatch {
    segment_index: usize,
    seed_hits: u32,
    local_score: i32,
    alignment: OrientedLocalAlignment,
}

struct SeededRead<'a> {
    read: &'a BamReadEvidence,
    reverse_sequence: Vec<u8>,
    seed_scores: HashMap<usize, u32>,
    best_by_kind: HashMap<SegmentKind, BestReadMatch>,
    alignments: HashMap<usize, OrientedLocalAlignment>,
}

struct SelectedSegment {
    support: GermlineSegmentSupport,
    umis: HashSet<String>,
}

pub struct PosteriorAnalyzer<'a> {
    reference: &'a VdjReference,
    mapper: &'a VdjMapper,
    config: PosteriorConfig,
}
impl<'a> PosteriorAnalyzer<'a> {
    pub fn new(
        reference: &'a VdjReference,
        mapper: &'a VdjMapper,
        config: PosteriorConfig,
    ) -> Self {
        Self {
            reference,
            mapper,
            config,
        }
    }

    pub fn analyze<E: ExpressionMatrix + Sync>(
        &self,
        reads: impl IntoIterator<Item = BamReadEvidence>,
        gex: &E,
    ) -> Vec<CellVdjSummary> {
        let mut cells: HashMap<String, Vec<BamReadEvidence>> = HashMap::new();
        for r in reads {
            if !r.is_secondary {
                cells.entry(r.cell.clone()).or_default().push(r);
            }
        }
        let mut out: Vec<_> = cells
            .into_par_iter()
            .map(|(cell, reads)| self.analyze_cell(cell, reads, gex))
            .collect();
        out.sort_by(|a, b| a.cell.cmp(&b.cell));
        out
    }

    fn analyze_cell<E: ExpressionMatrix>(
        &self,
        cell: String,
        reads: Vec<BamReadEvidence>,
        gex: &E,
    ) -> CellVdjSummary {
        let seeded: Vec<SeededRead<'_>> = reads
            .par_iter()
            .map(|read| {
                let reverse_sequence = crate::sequence::reverse_complement(&read.sequence);
                let ranked = self
                    .mapper
                    .segment_seed_scores_ranked_oriented(&read.sequence, &reverse_sequence);
                let seed_scores: HashMap<usize, u32> = ranked
                    .iter()
                    .map(|&(idx, hits, _)| (idx, hits))
                    .collect();
                let (best_by_kind, alignments) =
                    self.best_read_matches(read, &reverse_sequence, &ranked);
                SeededRead {
                    read,
                    reverse_sequence,
                    seed_scores,
                    best_by_kind,
                    alignments,
                }
            })
            .collect();

        let mut sterile_acc: HashMap<Chain, SterileAccumulator> = Chain::ALL
            .into_iter()
            .filter_map(|c| {
                SterileAccumulator::new(self.reference, c, self.config.sterile_bins).map(|a| (c, a))
            })
            .collect();
        let mut chain_reads: HashMap<Chain, Vec<&SeededRead<'_>>> = HashMap::new();
        for evidence in &seeded {
            let mut seen_chain = HashSet::new();
            for best in evidence.best_by_kind.values() {
                let chain = self.reference.segments[best.segment_index].chain;
                if seen_chain.insert(chain) {
                    chain_reads.entry(chain).or_default().push(evidence);
                }
            }
            if self.mapper.map(&evidence.read.sequence).is_none() {
                if let Some(chr) = &evidence.read.chr {
                    for &(s, e) in &evidence.read.ref_blocks {
                        for (chain, acc) in sterile_acc.iter_mut() {
                            if let Some((lchr, ls, le)) = self.reference.locus_bounds(*chain) {
                                if chr == lchr && e > ls && s < le {
                                    acc.observe(chr, s, e, &evidence.read.umi);
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut rearrangements: Vec<_> = Chain::ALL
            .par_iter()
            .filter_map(|&chain| {
                chain_reads
                    .get(&chain)
                    .and_then(|rs| self.call_chain(chain, rs))
            })
            .collect();
        rearrangements.sort_by_key(|call| call.chain);

        let sterile: Vec<_> = Chain::ALL
            .into_iter()
            .filter_map(|c| sterile_acc.remove(&c).map(|a| a.finish()))
            .collect();
        let b_locus = sterile
            .iter()
            .filter(|p| p.chain.is_bcr())
            .map(locus_signal)
            .fold(0.0, f64::max);
        let t_locus = sterile
            .iter()
            .filter(|p| p.chain.is_tcr())
            .map(locus_signal)
            .fold(0.0, f64::max);
        let development = score_development(gex, &cell, b_locus, t_locus);
        CellVdjSummary {
            cell,
            rearrangements,
            sterile,
            development,
        }
    }

    fn best_read_matches(
        &self,
        read: &BamReadEvidence,
        reverse_sequence: &[u8],
        ranked: &[(usize, u32, u32)],
    ) -> (HashMap<SegmentKind, BestReadMatch>, HashMap<usize, OrientedLocalAlignment>) {
        let mut candidates: HashMap<SegmentKind, Vec<(usize, u32)>> = HashMap::new();
        for &(segment_index, seed_hits, _) in ranked {
            if seed_hits < self.config.min_seed_hits {
                continue;
            }
            let kind = self.reference.segments[segment_index].kind;
            let bucket = candidates.entry(kind).or_default();
            if bucket.len() < self.config.candidate_segments_per_kind {
                bucket.push((segment_index, seed_hits));
            }
        }

        let candidate_list: Vec<_> = candidates.into_values().flatten().collect();
        let aligned: Vec<_> = candidate_list
            .into_iter()
            .map(|(segment_index, seed_hits)| {
                let segment = &self.reference.segments[segment_index];
                let alignment = best_oriented_local_alignment_with_reverse(
                    &read.sequence,
                    reverse_sequence,
                    &segment.sequence,
                );
                (segment_index, seed_hits, alignment)
            })
            .collect();

        let mut best = HashMap::new();
        let mut alignments = HashMap::with_capacity(aligned.len());
        for (segment_index, seed_hits, alignment) in aligned {
            alignments.insert(segment_index, alignment);
            let local = alignment.alignment.score;
            if local < self.config.min_local_score {
                continue;
            }
            let segment = &self.reference.segments[segment_index];
            let candidate = BestReadMatch { segment_index, seed_hits, local_score: local, alignment };
            let replace = best.get(&segment.kind).map_or(true, |old: &BestReadMatch| {
                candidate.local_score > old.local_score
                    || (candidate.local_score == old.local_score
                        && (candidate.seed_hits > old.seed_hits
                            || (candidate.seed_hits == old.seed_hits
                                && candidate.segment_index < old.segment_index)))
            });
            if replace { best.insert(segment.kind, candidate); }
        }
        (best, alignments)
    }

    fn call_chain(&self, chain: Chain, reads: &[&SeededRead<'_>]) -> Option<RearrangementCall> {
        let (v, j) = rayon::join(
            || self.best_segment(chain, SegmentKind::V, reads, None),
            || self.best_segment(chain, SegmentKind::J, reads, None),
        );
        let (Some(v), Some(j)) = (v, j) else {
            return None;
        };

        // A rearrangement requires molecule-level coherence: the same UMI must
        // support the selected V and J.  This deliberately replaces the old
        // cell/chain-wide UMI count that produced support ~= every UMI in a cell.
        let candidate_vj = umi_intersection(&v.umis, &j.umis);
        let coherent_vj = self.coherent_vj_umis(
            reads,
            &candidate_vj,
            v.support.segment_index,
            j.support.segment_index,
        );
        if coherent_vj.is_empty() {
            return None;
        }

        // D is short and often seeds poorly. Search it only in UMIs already
        // shown to coherently support this V/J pair; do not require D to define
        // the chain support count.
        let d = if chain.has_d() {
            self.best_segment(chain, SegmentKind::D, reads, Some(&coherent_vj))
        } else {
            None
        };

        // C is supporting context, not evidence that creates a VJ call.  It must
        // be physically plausible and share at least one coherent VJ UMI.
        let c = self
            .best_segment(chain, SegmentKind::C, reads, None)
            .filter(|candidate| {
                candidate.support.distance_to_recombination_center
                    <= self.config.max_constant_distance_bp
                    && candidate.umis.iter().any(|umi| coherent_vj.contains(umi))
            });

        let stage = if chain.has_d() {
            RecombinationStage::Vdj
        } else {
            RecombinationStage::Vj
        };
        let notation = format_call(
            chain,
            Some(&v.support),
            d.as_ref().map(|x| &x.support),
            Some(&j.support),
            c.as_ref().map(|x| &x.support),
        );
        let supporting_reads = self.collect_supporting_reads(
            reads,
            &coherent_vj,
            v.support.segment_index,
            j.support.segment_index,
            d.as_ref().map(|x| x.support.segment_index),
            c.as_ref().map(|x| x.support.segment_index),
        );
        Some(RearrangementCall {
            chain,
            stage,
            v: Some(v.support),
            d: d.map(|x| x.support),
            j: Some(j.support),
            c: c.map(|x| x.support),
            total_supporting_umis: coherent_vj.len(),
            supporting_reads,
            notation,
        })
    }

    fn coherent_vj_umis(
        &self,
        reads: &[&SeededRead<'_>],
        candidate_umis: &HashSet<String>,
        v_index: usize,
        j_index: usize,
    ) -> HashSet<String> {
        let mut by_umi: HashMap<&str, Vec<&SeededRead<'_>>> = HashMap::new();
        for evidence in reads {
            if candidate_umis.contains(&evidence.read.umi) {
                by_umi
                    .entry(evidence.read.umi.as_str())
                    .or_default()
                    .push(*evidence);
            }
        }
        by_umi
            .into_iter()
            .filter_map(|(umi, evidence)| {
                self.umi_has_observed_vj_path(&evidence, v_index, j_index)
                    .then(|| umi.to_string())
            })
            .collect()
    }

    fn umi_has_observed_vj_path(
        &self,
        reads: &[&SeededRead<'_>],
        v_index: usize,
        j_index: usize,
    ) -> bool {
        let v_ref = &self.reference.segments[v_index].sequence;
        let j_ref = &self.reference.segments[j_index].sequence;
        let mut v_reads: Vec<(usize, Vec<u8>)> = Vec::new();
        let mut j_reads: Vec<(usize, Vec<u8>)> = Vec::new();

        for (read_idx, evidence) in reads.iter().enumerate() {
            let v_match = evidence
                .best_by_kind
                .get(&SegmentKind::V)
                .filter(|x| x.segment_index == v_index);
            let j_match = evidence
                .best_by_kind
                .get(&SegmentKind::J)
                .filter(|x| x.segment_index == j_index);

            if let (Some(v), Some(j)) = (v_match, j_match) {
                if self.valid_vj_geometry(v, j, v_ref.len()) {
                    return true;
                }
            }
            if let Some(v) = v_match {
                if self.v_reaches_junction(v, v_ref.len()) {
                    v_reads.push((read_idx, oriented_read(&evidence.read.sequence, v.alignment.reverse_complement)));
                }
            }
            if let Some(j) = j_match {
                if self.j_starts_at_junction(j) {
                    j_reads.push((read_idx, oriented_read(&evidence.read.sequence, j.alignment.reverse_complement)));
                }
            }
        }

        for (v_idx, v_read) in &v_reads {
            for (j_idx, j_read) in &j_reads {
                if v_idx == j_idx {
                    continue;
                }
                let Some(merged) = overlap_reads(v_read, j_read, self.config.min_vj_read_overlap)
                else {
                    continue;
                };
                let v = best_oriented_local_alignment(&merged, v_ref);
                let j = best_oriented_local_alignment(&merged, j_ref);
                if !v.reverse_complement
                    && !j.reverse_complement
                    && self.valid_oriented_vj_geometry(v, j, v_ref.len())
                {
                    return true;
                }
            }
        }
        false
    }

    fn valid_vj_geometry(&self, v: &BestReadMatch, j: &BestReadMatch, v_len: usize) -> bool {
        if v.alignment.reverse_complement != j.alignment.reverse_complement {
            return false;
        }
        self.valid_oriented_vj_geometry(v.alignment, j.alignment, v_len)
    }

    fn valid_oriented_vj_geometry(
        &self,
        v: OrientedLocalAlignment,
        j: OrientedLocalAlignment,
        v_len: usize,
    ) -> bool {
        if !self.v_alignment_reaches_junction(v, v_len) || !self.j_alignment_starts_at_junction(j) {
            return false;
        }
        let va = v.alignment;
        let ja = j.alignment;
        if va.query_start > ja.query_start {
            return false;
        }
        va.query_end.saturating_sub(ja.query_start) <= self.config.max_vj_alignment_overlap
    }

    fn v_reaches_junction(&self, v: &BestReadMatch, v_len: usize) -> bool {
        self.v_alignment_reaches_junction(v.alignment, v_len)
    }

    fn v_alignment_reaches_junction(&self, v: OrientedLocalAlignment, v_len: usize) -> bool {
        v_len.saturating_sub(v.alignment.reference_end) <= self.config.max_v_end_distance
    }

    fn j_starts_at_junction(&self, j: &BestReadMatch) -> bool {
        self.j_alignment_starts_at_junction(j.alignment)
    }

    fn j_alignment_starts_at_junction(&self, j: OrientedLocalAlignment) -> bool {
        j.alignment.reference_start <= self.config.max_j_start_distance
    }

    fn collect_supporting_reads(
        &self,
        reads: &[&SeededRead<'_>],
        coherent_umis: &HashSet<String>,
        v_index: usize,
        j_index: usize,
        d_index: Option<usize>,
        c_index: Option<usize>,
    ) -> Vec<RearrangementSupportingRead> {
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        for evidence in reads {
            if !coherent_umis.contains(&evidence.read.umi) {
                continue;
            }
            let v_match = evidence
                .best_by_kind
                .get(&SegmentKind::V)
                .filter(|x| x.segment_index == v_index);
            let j_match = evidence
                .best_by_kind
                .get(&SegmentKind::J)
                .filter(|x| x.segment_index == j_index);
            let c_match = c_index.and_then(|idx| {
                evidence
                    .best_by_kind
                    .get(&SegmentKind::C)
                    .filter(|x| x.segment_index == idx)
            });
            let d_match = d_index.and_then(|idx| {
                let alignment = evidence
                    .alignments
                    .get(&idx)
                    .copied()
                    .unwrap_or_else(|| {
                        best_oriented_local_alignment_with_reverse(
                            &evidence.read.sequence,
                            &evidence.reverse_sequence,
                            &self.reference.segments[idx].sequence,
                        )
                    });
                (alignment.alignment.score >= self.config.min_local_score).then_some(alignment)
            });
            let supports_v = v_match.is_some();
            let supports_j = j_match.is_some();
            let supports_c = c_match.is_some();
            let supports_d = d_match.is_some();
            if !(supports_v || supports_j || supports_d || supports_c) {
                continue;
            }
            let key = (
                evidence.read.read_name.clone(),
                evidence.read.sequence.clone(),
            );
            if !seen.insert(key) {
                continue;
            }
            out.push(RearrangementSupportingRead {
                umi: evidence.read.umi.clone(),
                read_name: evidence.read.read_name.clone(),
                sequence: evidence.read.sequence.clone(),
                bam_is_reverse: evidence.read.is_reverse,
                is_supplementary: evidence.read.is_supplementary,
                supports_v,
                supports_j,
                supports_d,
                supports_c,
                v_alignment: v_match.map(|x| read_segment_alignment(x.alignment)),
                j_alignment: j_match.map(|x| read_segment_alignment(x.alignment)),
                d_alignment: d_match.map(read_segment_alignment),
                c_alignment: c_match.map(|x| read_segment_alignment(x.alignment)),
            });
        }
        out.sort_by(|a, b| {
            a.umi
                .cmp(&b.umi)
                .then_with(|| a.read_name.cmp(&b.read_name))
                .then_with(|| a.sequence.cmp(&b.sequence))
        });
        out
    }

    fn best_segment(
        &self,
        chain: Chain,
        kind: SegmentKind,
        reads: &[&SeededRead<'_>],
        restrict_umis: Option<&HashSet<String>>,
    ) -> Option<SelectedSegment> {
        let mut seed_totals: HashMap<usize, u32> = HashMap::new();
        for evidence in reads {
            if restrict_umis.is_some_and(|allowed| !allowed.contains(&evidence.read.umi)) {
                continue;
            }
            if kind == SegmentKind::D {
                for (&idx, &score) in &evidence.seed_scores {
                    let s = &self.reference.segments[idx];
                    if s.chain == chain && s.kind == kind {
                        *seed_totals.entry(idx).or_insert(0) += score;
                    }
                }
            } else if let Some(best) = evidence.best_by_kind.get(&kind) {
                let s = &self.reference.segments[best.segment_index];
                if s.chain == chain {
                    *seed_totals.entry(best.segment_index).or_insert(0) += best.seed_hits;
                }
            }
        }
        let mut cand: Vec<_> = seed_totals.into_iter().collect();
        cand.sort_by_key(|x| std::cmp::Reverse(x.1));
        cand.truncate(self.config.candidate_segments_per_kind);
        // D segments are short and can seed poorly; fall back to the complete,
        // usually small D set, but only on coherent V/J UMIs.
        if cand.is_empty() && kind == SegmentKind::D {
            cand = self
                .reference
                .segments_for(chain, kind)
                .map(|(i, _)| (i, 0))
                .collect();
        }

        cand.into_par_iter()
            .filter_map(|(idx, _)| {
                let seg = &self.reference.segments[idx];
                if kind == SegmentKind::C
                    && seg.distance_to_recombination_center > self.config.max_constant_distance_bp
                {
                    return None;
                }

                let mut total = 0i32;
                let mut umis = HashSet::new();
                let mut supporting_reads = 0usize;
                for evidence in reads {
                    if restrict_umis.is_some_and(|allowed| !allowed.contains(&evidence.read.umi)) {
                        continue;
                    }

                    // V/J/C evidence is assigned to the best local-alignment match
                    // for this read across all receptor chains. Shared k-mers therefore
                    // nominate candidates but cannot make the same read support every
                    // homologous IG/TR segment. D remains a deliberate exception because
                    // its short germline sequence can be seed-poor after junctional editing.
                    let score = if kind == SegmentKind::D {
                        evidence
                            .alignments
                            .get(&idx)
                            .copied()
                            .unwrap_or_else(|| {
                                best_oriented_local_alignment_with_reverse(
                                    &evidence.read.sequence,
                                    &evidence.reverse_sequence,
                                    &seg.sequence,
                                )
                            })
                            .alignment
                            .score
                    } else {
                        let Some(best) = evidence.best_by_kind.get(&kind) else {
                            continue;
                        };
                        if best.segment_index != idx {
                            continue;
                        }
                        best.local_score
                    };
                    if score >= self.config.min_local_score {
                        total += score;
                        umis.insert(evidence.read.umi.clone());
                        supporting_reads += 1;
                    }
                }
                if supporting_reads == 0 {
                    return None;
                }
                Some(SelectedSegment {
                    support: GermlineSegmentSupport {
                        segment_index: idx,
                        id: seg.name.clone(),
                        kind,
                        local_alignment_score: total,
                        supporting_umis: umis.len(),
                        supporting_reads,
                        locus_fraction: seg.locus_fraction,
                        distance_to_recombination_center: seg.distance_to_recombination_center,
                    },
                    umis,
                })
            })
            .max_by(|a, b| {
                a.support
                    .local_alignment_score
                    .cmp(&b.support.local_alignment_score)
                    .then_with(|| a.support.supporting_umis.cmp(&b.support.supporting_umis))
                    // Reverse the index comparison so the lower stable reference
                    // index wins exact ties independent of Rayon scheduling.
                    .then_with(|| b.support.segment_index.cmp(&a.support.segment_index))
            })
    }
}

fn read_segment_alignment(hit: OrientedLocalAlignment) -> ReadSegmentAlignment {
    ReadSegmentAlignment {
        score: hit.alignment.score,
        query_start: hit.alignment.query_start,
        query_end: hit.alignment.query_end,
        reference_start: hit.alignment.reference_start,
        reference_end: hit.alignment.reference_end,
        reverse_complement: hit.reverse_complement,
    }
}

fn oriented_read(sequence: &[u8], reverse: bool) -> Vec<u8> {
    if reverse {
        crate::sequence::reverse_complement(sequence)
    } else {
        sequence.to_vec()
    }
}

fn overlap_reads(left: &[u8], right: &[u8], min_overlap: usize) -> Option<Vec<u8>> {
    let max_overlap = left.len().min(right.len());
    if max_overlap < min_overlap {
        return None;
    }
    for overlap in (min_overlap..=max_overlap).rev() {
        let a = &left[left.len() - overlap..];
        let b = &right[..overlap];
        let mismatches = a.iter().zip(b).filter(|(x, y)| x != y).count();
        if mismatches > (overlap / 20).max(1) {
            continue;
        }
        let mut merged = Vec::with_capacity(left.len() + right.len() - overlap);
        merged.extend_from_slice(&left[..left.len() - overlap]);
        for (&x, &y) in a.iter().zip(b) {
            merged.push(if x == y { x } else { b'N' });
        }
        merged.extend_from_slice(&right[overlap..]);
        return Some(merged);
    }
    None
}

fn umi_intersection(a: &HashSet<String>, b: &HashSet<String>) -> HashSet<String> {
    if a.len() <= b.len() {
        a.iter().filter(|umi| b.contains(*umi)).cloned().collect()
    } else {
        b.iter().filter(|umi| a.contains(*umi)).cloned().collect()
    }
}

fn locus_signal(p: &SterileProfile) -> f64 {
    let umi = (p.total_unique_umis as f64 / 4.0).min(1.0);
    (0.50 * umi + 0.35 * p.breadth + 0.15 * p.distal_fraction).clamp(0.0, 1.0)
}
fn support_id(support: Option<&GermlineSegmentSupport>) -> &str {
    support.map_or("?", |support| support.id.as_str())
}

fn format_call(
    chain: Chain,
    v: Option<&GermlineSegmentSupport>,
    d: Option<&GermlineSegmentSupport>,
    j: Option<&GermlineSegmentSupport>,
    c: Option<&GermlineSegmentSupport>,
) -> String {
    if chain.has_d() {
        format!(
            "{}:{}-({})-{}-{}",
            chain,
            support_id(v),
            support_id(d),
            support_id(j),
            support_id(c)
        )
    } else {
        format!(
            "{}:{}-{}-{}",
            chain,
            support_id(v),
            support_id(j),
            support_id(c)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mapper::VdjMapperConfig;
    use crate::types::VdjSegment;

    fn segment(name: &str, chain: Chain, sequence: &[u8]) -> VdjSegment {
        VdjSegment {
            name: name.to_string(),
            transcript_id: format!("{name}_tx"),
            gene_id: name.to_string(),
            chain,
            kind: SegmentKind::V,
            chr: "chr1".to_string(),
            start: 0,
            end: sequence.len() as u32,
            strand_minus: false,
            locus_rank: 0,
            locus_fraction: 0.0,
            distance_to_recombination_center: 0,
            sequence: sequence.to_vec(),
        }
    }

    #[test]
    fn umi_intersection_is_not_union() {
        let a = HashSet::from(["u1".to_string(), "u2".to_string()]);
        let b = HashSet::from(["u2".to_string(), "u3".to_string()]);
        assert_eq!(umi_intersection(&a, &b), HashSet::from(["u2".to_string()]));
    }

    #[test]
    fn shared_seed_read_supports_only_best_cross_chain_segment() {
        let reference = VdjReference {
            segments: vec![
                segment("IGHV_best", Chain::Igh, b"AACCGGTTAAAA"),
                segment("TRAV_shared", Chain::Tra, b"AACCGGTTCCCC"),
            ],
        };
        let mapper = VdjMapper::new(
            reference.clone(),
            VdjMapperConfig {
                k: 5,
                min_v_hits: 1,
                min_j_hits: 1,
                min_c_hits: 1,
            },
        );
        let analyzer = PosteriorAnalyzer::new(
            &reference,
            &mapper,
            PosteriorConfig {
                min_seed_hits: 1,
                min_local_score: 4,
                ..PosteriorConfig::default()
            },
        );
        let read = BamReadEvidence {
            cell: "cell".to_string(),
            umi: "umi".to_string(),
            read_name: "read".to_string(),
            sequence: b"AACCGGTTAAAA".to_vec(),
            chr: None,
            ref_start: None,
            ref_end: None,
            ref_blocks: Vec::new(),
            mapq: 0,
            is_reverse: false,
            is_secondary: false,
            is_supplementary: false,
        };
        let reverse = crate::sequence::reverse_complement(&read.sequence);
        let ranked = mapper.segment_seed_scores_ranked_oriented(&read.sequence, &reverse);
        let seed_scores: HashMap<usize, u32> = ranked
            .iter()
            .map(|&(idx, hits, _)| (idx, hits))
            .collect();
        assert!(seed_scores.contains_key(&0));
        assert!(seed_scores.contains_key(&1));

        let (best, _) = analyzer.best_read_matches(&read, &reverse, &ranked);
        assert_eq!(best[&SegmentKind::V].segment_index, 0);
    }
}
