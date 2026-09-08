use crate::cellrep::{CellEvidence, ReceptorSequenceEvidence};
use crate::index::{
    reference_base_matches, reverse_complement, Chain, SegmentId, SegmentKind, VdjIndex,
};
use anyhow::{bail, Result};

mod constant_region_linkage;
mod identifier;
mod recombination_evidence_rescan;
pub use identifier::{ReceptorRole, RecombinationId};
pub(crate) use recombination_evidence_rescan::{
    rescue_missing_constants_from_bam, rescue_missing_constants_from_bam_with_report,
    rescue_missing_constants_from_bam_with_report_and_progress,
};
pub use recombination_evidence_rescan::{
    RecombinationEvidenceRescanProgress, RecombinationEvidenceRescanReport,
    RecombinationRescanCall, RecombinationRescanCandidate,
};

use constant_region_linkage::rescue_missing_constant_regions;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalAlignment {
    pub score: i32,
    pub query_start: usize,
    pub query_end: usize,
    pub reference_start: usize,
    pub reference_end: usize,
}

pub fn local_alignment(query: &[u8], reference: &[u8]) -> LocalAlignment {
    #[derive(Clone, Copy, Default)]
    struct Cell {
        score: i32,
        qs: usize,
        rs: usize,
    }
    if query.is_empty() || reference.is_empty() {
        return LocalAlignment {
            score: 0,
            query_start: 0,
            query_end: 0,
            reference_start: 0,
            reference_end: 0,
        };
    }
    let mut prev = vec![Cell::default(); reference.len() + 1];
    let mut curr = vec![Cell::default(); reference.len() + 1];
    let mut best = LocalAlignment {
        score: 0,
        query_start: 0,
        query_end: 0,
        reference_start: 0,
        reference_end: 0,
    };
    for (i, &q) in query.iter().enumerate() {
        curr[0] = Cell::default();
        for (j, &r) in reference.iter().enumerate() {
            let ms = if reference_base_matches(r, q) { 2 } else { -2 };
            let diag = prev[j].score + ms;
            let up = prev[j + 1].score - 3;
            let left = curr[j].score - 3;
            let mut c = Cell::default();
            if diag > 0 && diag >= up && diag >= left {
                c.score = diag;
                if prev[j].score > 0 {
                    c.qs = prev[j].qs;
                    c.rs = prev[j].rs
                } else {
                    c.qs = i;
                    c.rs = j
                }
            } else if up > 0 && up >= left {
                c.score = up;
                c.qs = prev[j + 1].qs;
                c.rs = prev[j + 1].rs
            } else if left > 0 {
                c.score = left;
                c.qs = curr[j].qs;
                c.rs = curr[j].rs
            }
            curr[j + 1] = c;
            let cand = LocalAlignment {
                score: c.score,
                query_start: c.qs,
                query_end: i + 1,
                reference_start: c.rs,
                reference_end: j + 1,
            };
            if cand.score > best.score
                || (cand.score == best.score
                    && cand.score > 0
                    && (
                        cand.query_start,
                        cand.query_end,
                        cand.reference_start,
                        cand.reference_end,
                    ) < (
                        best.query_start,
                        best.query_end,
                        best.reference_start,
                        best.reference_end,
                    ))
            {
                best = cand
            }
        }
        std::mem::swap(&mut prev, &mut curr);
        curr.fill(Cell::default());
    }
    best
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JunctionStructure {
    pub v_del_3: u16,
    pub p_v3: Vec<u8>,
    pub n1: Vec<u8>,
    pub p_d5: Vec<u8>,
    pub d_del_5: Option<u16>,
    pub d_retained: Vec<u8>,
    pub d_del_3: Option<u16>,
    pub p_d3: Vec<u8>,
    pub n2: Vec<u8>,
    pub p_j5: Vec<u8>,
    pub j_del_5: u16,
    pub pn_alternative: bool,
}
impl JunctionStructure {
    pub fn p_v3_len(&self) -> u16 {
        self.p_v3.len() as u16
    }
    pub fn n1_len(&self) -> u16 {
        self.n1.len() as u16
    }
    pub fn p_d5_len(&self) -> u16 {
        self.p_d5.len() as u16
    }
    pub fn d_retained_len(&self) -> u16 {
        self.d_retained.len() as u16
    }
    pub fn p_d3_len(&self) -> u16 {
        self.p_d3.len() as u16
    }
    pub fn n2_len(&self) -> u16 {
        self.n2.len() as u16
    }
    pub fn p_j5_len(&self) -> u16 {
        self.p_j5.len() as u16
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstantRegionEvidence {
    pub segment: SegmentId,
    pub supporting_features: u32,
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReceptorLinkageSupport {
    /// BAM records that independently rediscover this cell-specific CDR3/J bait.
    pub rediscovery_reads: u32,
    /// Rescan reads with coordinate-safe overlap of the AIRR junction.
    pub junction_support_reads: u32,
    /// Rescan reads whose coordinate-safe alignment spans the complete AIRR junction.
    pub junction_spanning_reads: u32,
    /// Junction-supporting reads that disagree with the pre-rescan junction at >=1 base.
    pub junction_conflicting_reads: u32,
    /// Junction bases changed by the Stage-3 read consensus.
    pub junction_refined_bases: u16,
    /// Physical fragments that link the receptor bait to a constant region.
    pub constant_link_fragments: u32,
    /// BAM records that contain both receptor-bait and constant-region support.
    pub constant_spanning_reads: u32,
    /// Constant segment supported by the linkage evidence, if any.
    pub constant_segment: Option<SegmentId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductivityStatus {
    Productive,
    StopCodon,
    MissingVCodingStart,
    MissingVAnchor,
    MissingJAnchor,
    InvalidAnchorOrder,
    TooShort,
}

impl ProductivityStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Productive => "productive",
            Self::StopCodon => "unproductive_stop_codon",
            Self::MissingVCodingStart => "unknown_no_v_cds",
            Self::MissingVAnchor => "unproductive_missing_v_anchor",
            Self::MissingJAnchor => "unproductive_missing_j_anchor",
            Self::InvalidAnchorOrder => "unproductive_invalid_anchor_order",
            Self::TooShort => "unproductive_sequence_too_short",
        }
    }

    pub fn is_unknown(self) -> bool {
        matches!(self, Self::MissingVCodingStart)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Recombination {
    pub chain: Chain,
    pub v: SegmentId,
    pub d: Option<SegmentId>,
    pub j: SegmentId,
    pub junction: JunctionStructure,
    pub constant: Option<ConstantRegionEvidence>,
    pub observed_rearrangement: Vec<u8>,
    pub naive_recombination: Vec<u8>,
    pub observed_receptor_sequence: Vec<u8>,
    pub productive: bool,
    pub in_frame: bool,
    pub stop_codon: bool,
    /// Why productivity could or could not be determined from the observed receptor.
    pub productivity_status: ProductivityStatus,
    /// AIRR junction: conserved V cysteine through conserved J F/W, inclusive.
    pub airr_junction: Vec<u8>,
    pub airr_junction_aa: Vec<u8>,
    /// AIRR CDR3 excludes the conserved V cysteine and J F/W anchors.
    pub cdr3: Vec<u8>,
    pub cdr3_aa: Vec<u8>,
    pub supporting_features: u32,
    pub receptor_linkage: ReceptorLinkageSupport,
    pub stable_id: RecombinationId,
}

pub fn identify_cell_recombinations(
    cell: &CellEvidence,
    index: &VdjIndex,
    min_overlap: usize,
) -> Vec<Recombination> {
    let work: Vec<_> = cell
        .chains(index)
        .into_iter()
        .map(|chain| ChainWork {
            cell,
            index,
            chain,
            min_overlap,
        })
        .collect();
    work.into_iter().flat_map(process_chain_work).collect()
}

/// One independently executable cell/locus unit.  This is deliberately kept
/// free of shared mutable state so callers can later collect `ChainWork` items
/// and process them with Rayon without touching the assembly implementation.
#[derive(Clone, Copy)]
pub struct ChainWork<'a> {
    pub cell: &'a CellEvidence,
    pub index: &'a VdjIndex,
    pub chain: Chain,
    pub min_overlap: usize,
}

/// Perform all summary collapsing and recombination identification for exactly
/// one cell/locus work unit.
pub fn process_chain_work(work: ChainWork<'_>) -> Vec<Recombination> {
    const MAX_FINAL_RECOMBINATIONS_PER_LOCUS: usize = 2;

    let summaries = work.cell.summaries_for_chain(work.chain);

    let mut calls: Vec<_> = summaries
        .iter()
        .filter_map(|summary| identify_summary(summary, work.index, work.chain))
        .collect();

    // The transient summary count may exceed two because unlinked fragments are
    // intentionally kept apart.  Biology constrains the final resolved calls.
    calls.sort_by_key(|r| std::cmp::Reverse(r.supporting_features));
    let mut seen_sequences = std::collections::HashSet::<Vec<u8>>::new();
    calls.retain(|r| seen_sequences.insert(r.observed_rearrangement.clone()));
    calls.truncate(MAX_FINAL_RECOMBINATIONS_PER_LOCUS);

    rescue_missing_constant_regions(&mut calls, work.cell, work.index, work.chain);
    calls
}

fn identify_summary(
    assembled: &ReceptorSequenceEvidence,
    index: &VdjIndex,
    chain: Chain,
) -> Option<Recombination> {
    let raw = assembled.consensus(index);
    if raw.is_empty() {
        return None;
    }

    // Evaluate both transcript orientations as complete hypotheses.  The
    // summary itself is normally already transcript-oriented from mapper
    // geometry, but retaining this check makes the caller robust to ambiguous
    // or synthetic evidence.
    let rc = reverse_complement(&raw);
    let f = best_vj(&raw, index, chain, assembled)?;
    let r = best_vj(&rc, index, chain, assembled)?;
    let (observed, v_id, v_aln, j_id, j_aln) = if f.4 >= r.4 {
        (raw, f.0, f.1, f.2, f.3)
    } else {
        (rc, r.0, r.1, r.2, r.3)
    };
    if v_aln.query_end > j_aln.query_start {
        return None;
    }

    let (d_id, d_aln) = if chain.has_d() {
        best_d_bounded(
            &observed,
            index,
            chain,
            v_aln.query_end,
            j_aln.query_start,
            assembled,
        )
    } else {
        (None, None)
    };
    if chain.has_d() && d_id.is_none() {
        return None;
    }

    let v = &index.segment(v_id)?.sequence;
    let j = &index.segment(j_id)?.sequence;
    let d = d_id
        .and_then(|id| index.segment(id))
        .map(|s| s.sequence.as_slice());
    let measured = measure_junction(&observed, v, v_aln, d, d_aln, j, j_aln)?;
    let constant = best_constant(&observed, index, chain, j_aln.query_end, assembled);
    let rearr_end = j_aln.query_end.min(observed.len());
    let observed_rearrangement = observed[..rearr_end].to_vec();
    let observed_receptor_sequence = observed.clone();
    let productivity = assess_productivity(
        &observed,
        v_aln,
        j_aln,
        index.segment(v_id).and_then(|s| s.coding_start()),
    );
    let in_frame = productivity.in_frame;
    let stop = productivity.stop_codon;
    let productive = productivity.productive;
    let junction = measured.junction;
    let mut recomb = Recombination {
        chain,
        v: v_id,
        d: d_id,
        j: j_id,
        junction,
        constant,
        observed_rearrangement,
        naive_recombination: measured.naive,
        observed_receptor_sequence,
        productive,
        in_frame,
        stop_codon: stop,
        productivity_status: productivity.status,
        airr_junction: productivity.junction.clone(),
        airr_junction_aa: productivity.junction_aa.clone(),
        cdr3: productivity.cdr3.clone(),
        cdr3_aa: productivity.cdr3_aa.clone(),
        supporting_features: assembled.support_features,
        receptor_linkage: ReceptorLinkageSupport::default(),
        stable_id: RecombinationId::placeholder(chain),
    };
    recomb.stable_id = RecombinationId::from_recombination(&recomb, index).ok()?;
    Some(recomb)
}

pub(crate) fn refresh_recombination_from_observed(
    recomb: &mut Recombination,
    index: &VdjIndex,
) -> bool {
    let observed = recomb.observed_rearrangement.clone();
    let Some(v_seg) = index.segment(recomb.v) else { return false; };
    let Some(j_seg) = index.segment(recomb.j) else { return false; };
    let v_aln = local_alignment(&observed, &v_seg.sequence);
    let j_aln = local_alignment(&observed, &j_seg.sequence);
    if v_aln.score <= 0 || j_aln.score <= 0 || v_aln.query_end > j_aln.query_start {
        return false;
    }
    let (d_seq, d_aln) = if let Some(d_id) = recomb.d {
        let Some(d_seg) = index.segment(d_id) else { return false; };
        let aln = local_alignment(
            &observed[v_aln.query_end.min(observed.len())..j_aln.query_start.min(observed.len())],
            &d_seg.sequence,
        );
        if aln.score <= 0 {
            return false;
        }
        let mut shifted = aln;
        shifted.query_start += v_aln.query_end;
        shifted.query_end += v_aln.query_end;
        (Some(d_seg.sequence.as_slice()), Some(shifted))
    } else {
        (None, None)
    };
    let Some(measured) = measure_junction(
        &observed,
        &v_seg.sequence,
        v_aln,
        d_seq,
        d_aln,
        &j_seg.sequence,
        j_aln,
    ) else {
        return false;
    };
    let productivity = assess_productivity(
        &observed,
        v_aln,
        j_aln,
        v_seg.coding_start(),
    );
    recomb.junction = measured.junction;
    recomb.naive_recombination = measured.naive;
    recomb.productive = productivity.productive;
    recomb.in_frame = productivity.in_frame;
    recomb.stop_codon = productivity.stop_codon;
    recomb.productivity_status = productivity.status;
    recomb.airr_junction = productivity.junction;
    recomb.airr_junction_aa = productivity.junction_aa;
    recomb.cdr3 = productivity.cdr3;
    recomb.cdr3_aa = productivity.cdr3_aa;
    if let Ok(id) = RecombinationId::from_recombination(recomb, index) {
        recomb.stable_id = id;
    }
    true
}

fn best_vj(
    observed: &[u8],
    index: &VdjIndex,
    chain: Chain,
    summary: &ReceptorSequenceEvidence,
) -> Option<(SegmentId, LocalAlignment, SegmentId, LocalAlignment, i32)> {
    let (v, va) = best_segment(observed, index, chain, SegmentKind::V, None, summary)?;
    let (j, ja) = best_segment(observed, index, chain, SegmentKind::J, None, summary)?;
    Some((v, va, j, ja, va.score + ja.score))
}

fn best_segment(
    observed: &[u8],
    index: &VdjIndex,
    chain: Chain,
    kind: SegmentKind,
    window: Option<(usize, usize)>,
    summary: &ReceptorSequenceEvidence,
) -> Option<(SegmentId, LocalAlignment)> {
    let (a, b) = window.unwrap_or((0, observed.len()));
    if a >= b || b > observed.len() {
        return None;
    }
    let q = &observed[a..b];

    // Mapper-nominated segments are a strong, cheap prior.  Only fall back to
    // all reference segments if this summary lacks a direct nomination of the
    // requested kind.
    let nominated: Vec<_> = summary
        .segment_ids()
        .filter(|id| {
            index
                .segment(*id)
                .is_some_and(|s| s.chain == chain && s.kind == kind)
        })
        .collect();

    let score_one = |id: SegmentId| {
        let s = index.segment(id)?;
        let mut x = local_alignment(q, &s.sequence);
        x.query_start += a;
        x.query_end += a;
        let support_bonus = summary.segment_support(id).min((i32::MAX / 4) as u32) as i32 * 4;
        let score = x.score.saturating_add(support_bonus);
        Some((id, x, score))
    };

    if !nominated.is_empty() {
        nominated
            .into_iter()
            .filter_map(score_one)
            .max_by_key(|(id, _, score)| (*score, std::cmp::Reverse(*id)))
            .map(|(id, aln, _)| (id, aln))
            .filter(|(_, a)| a.score > 0)
    } else {
        index
            .segments_for(chain, kind)
            .filter_map(|s| score_one(s.id))
            .max_by_key(|(id, _, score)| (*score, std::cmp::Reverse(*id)))
            .map(|(id, aln, _)| (id, aln))
            .filter(|(_, a)| a.score > 0)
    }
}

fn best_d_bounded(
    observed: &[u8],
    index: &VdjIndex,
    chain: Chain,
    a: usize,
    b: usize,
    summary: &ReceptorSequenceEvidence,
) -> (Option<SegmentId>, Option<LocalAlignment>) {
    if a >= b || b > observed.len() {
        return (None, None);
    }
    let nominated: Vec<_> = summary
        .segment_ids()
        .filter(|id| {
            index
                .segment(*id)
                .is_some_and(|s| s.chain == chain && s.kind == SegmentKind::D)
        })
        .collect();

    let eval = |id: SegmentId| {
        let s = index.segment(id)?;
        let mut x = local_alignment(&observed[a..b], &s.sequence);
        x.query_start += a;
        x.query_end += a;
        let bonus = summary.segment_support(id).min((i32::MAX / 4) as u32) as i32 * 4;
        let score = x.score.saturating_add(bonus);
        Some((id, x, score))
    };

    let best = if !nominated.is_empty() {
        nominated
            .into_iter()
            .filter_map(|id| eval(id))
            .max_by_key(|(id, _, score)| (*score, std::cmp::Reverse(*id)))
    } else {
        index
            .segments_for(chain, SegmentKind::D)
            .filter_map(|s| eval(s.id))
            .max_by_key(|(id, _, score)| (*score, std::cmp::Reverse(*id)))
    };

    match best {
        Some((id, a, score)) if score >= 4 => (Some(id), Some(a)),
        _ => (None, None),
    }
}

fn best_constant(
    observed: &[u8],
    index: &VdjIndex,
    chain: Chain,
    j_end: usize,
    summary: &ReceptorSequenceEvidence,
) -> Option<ConstantRegionEvidence> {
    let nominated: Vec<_> = summary
        .segment_ids()
        .filter(|id| {
            index
                .segment(*id)
                .is_some_and(|s| s.chain == chain && s.kind == SegmentKind::C)
        })
        .collect();

    let eval = |id: SegmentId| {
        let s = index.segment(id)?;
        let aln = local_alignment(&observed[j_end.min(observed.len())..], &s.sequence);
        let support = summary.segment_support(id);
        let bonus = support.min((i32::MAX / 8) as u32) as i32 * 8;
        Some((id, aln.score.saturating_add(bonus), support))
    };

    let best = if !nominated.is_empty() {
        nominated
            .into_iter()
            .filter_map(|id| eval(id))
            .max_by_key(|(id, score, _)| (*score, std::cmp::Reverse(*id)))
    } else {
        index
            .segments_for(chain, SegmentKind::C)
            .filter_map(|s| eval(s.id))
            .max_by_key(|(id, score, _)| (*score, std::cmp::Reverse(*id)))
    }?;

    if best.1 <= 0 {
        return None;
    }
    let s = index.segment(best.0)?;
    Some(ConstantRegionEvidence {
        segment: best.0,
        supporting_features: best.2,
        sequence: s.sequence.clone(),
    })
}

struct Measured {
    junction: JunctionStructure,
    naive: Vec<u8>,
}
fn measure_junction(
    observed: &[u8],
    v: &[u8],
    va: LocalAlignment,
    d: Option<&[u8]>,
    da: Option<LocalAlignment>,
    j: &[u8],
    ja: LocalAlignment,
) -> Option<Measured> {
    if va.reference_end > v.len()
        || ja.reference_end > j.len()
        || va.query_end > observed.len()
        || ja.query_start > observed.len()
    {
        return None;
    }
    let v_del_3 = u16::try_from(v.len().checked_sub(va.reference_end)?).ok()?;
    let j_del_5 = u16::try_from(ja.reference_start).ok()?;
    let (p_v3, n1, p_d5, d_del_5, d_retained, d_del_3, p_d3, n2, p_j5, ambiguous) =
        if let (Some(dref), Some(da)) = (d, da) {
            if va.query_end > da.query_start
                || da.query_end > ja.query_start
                || da.reference_end > dref.len()
            {
                return None;
            }
            let ls = split_junction(
                &observed[va.query_end..da.query_start],
                &v[..va.reference_end],
                &dref[da.reference_start..da.reference_end],
            );
            let rs = split_junction(
                &observed[da.query_end..ja.query_start],
                &dref[da.reference_start..da.reference_end],
                &j[ja.reference_start..],
            );
            (
                ls.left_p,
                ls.n,
                ls.right_p,
                Some(da.reference_start as u16),
                dref[da.reference_start..da.reference_end].to_vec(),
                Some((dref.len() - da.reference_end) as u16),
                rs.left_p,
                rs.n,
                rs.right_p,
                ls.ambiguous || rs.ambiguous,
            )
        } else {
            if va.query_end > ja.query_start {
                return None;
            }
            let s = split_junction(
                &observed[va.query_end..ja.query_start],
                &v[..va.reference_end],
                &j[ja.reference_start..],
            );
            (
                s.left_p,
                s.n,
                Vec::new(),
                None,
                Vec::new(),
                None,
                Vec::new(),
                Vec::new(),
                s.right_p,
                s.ambiguous,
            )
        };
    let mut naive = Vec::new();
    naive.extend_from_slice(&v[..va.reference_end]);
    naive.extend_from_slice(&p_v3);
    naive.extend_from_slice(&n1);
    naive.extend_from_slice(&p_d5);
    naive.extend_from_slice(&d_retained);
    naive.extend_from_slice(&p_d3);
    naive.extend_from_slice(&n2);
    naive.extend_from_slice(&p_j5);
    naive.extend_from_slice(&j[ja.reference_start..]);
    Some(Measured {
        junction: JunctionStructure {
            v_del_3,
            p_v3,
            n1,
            p_d5,
            d_del_5,
            d_retained,
            d_del_3,
            p_d3,
            n2,
            p_j5,
            j_del_5,
            pn_alternative: ambiguous,
        },
        naive,
    })
}
const MAX_P: usize = 12;
struct Split {
    left_p: Vec<u8>,
    n: Vec<u8>,
    right_p: Vec<u8>,
    ambiguous: bool,
}
fn split_junction(junction: &[u8], left: &[u8], right: &[u8]) -> Split {
    let lm = MAX_P.min(left.len()).min(junction.len());
    let rm = MAX_P.min(right.len()).min(junction.len());
    let mut c = Vec::new();
    for lp in 0..=lm {
        let el = reverse_complement(&left[left.len() - lp..]);
        if junction.get(..lp) != Some(el.as_slice()) {
            continue;
        }
        for rp in 0..=rm.min(junction.len() - lp) {
            let er = reverse_complement(&right[..rp]);
            if junction.get(junction.len() - rp..) != Some(er.as_slice()) {
                continue;
            }
            c.push((lp, rp))
        }
    }
    c.sort_by(|a, b| {
        (b.0 + b.1)
            .cmp(&(a.0 + a.1))
            .then_with(|| b.0.cmp(&a.0))
            .then_with(|| b.1.cmp(&a.1))
    });
    let (lp, rp) = c.first().copied().unwrap_or((0, 0));
    Split {
        left_p: junction[..lp].to_vec(),
        n: junction[lp..junction.len() - rp].to_vec(),
        right_p: junction[junction.len() - rp..].to_vec(),
        ambiguous: c.len() > 1,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProductivityAssessment {
    productive: bool,
    in_frame: bool,
    stop_codon: bool,
    status: ProductivityStatus,
    junction: Vec<u8>,
    junction_aa: Vec<u8>,
    cdr3: Vec<u8>,
    cdr3_aa: Vec<u8>,
}

fn empty_productivity(status: ProductivityStatus) -> ProductivityAssessment {
    ProductivityAssessment {
        productive: false,
        in_frame: false,
        stop_codon: false,
        status,
        junction: Vec::new(),
        junction_aa: Vec::new(),
        cdr3: Vec::new(),
        cdr3_aa: Vec::new(),
    }
}

/// Determine receptor productivity in the coding frame fixed by the annotated
/// V CDS start.  We no longer choose whichever of three frames happens to make
/// conserved receptor motifs look plausible.
fn assess_productivity(
    observed: &[u8],
    v_aln: LocalAlignment,
    j_aln: LocalAlignment,
    v_coding_start: Option<usize>,
) -> ProductivityAssessment {
    const V_ANCHOR_WINDOW: usize = 90;
    const J_ANCHOR_WINDOW: usize = 120;

    let Some(cds_start) = v_coding_start else {
        return empty_productivity(ProductivityStatus::MissingVCodingStart);
    };
    let rearr_end = j_aln.query_end.min(observed.len());
    if rearr_end < 3 {
        return empty_productivity(ProductivityStatus::TooShort);
    }

    // Project the V reference coding origin through the local alignment. Only
    // the residue modulo three matters, so this remains valid when the CDS
    // start lies outside the aligned V fragment.
    let delta = v_aln.reference_start as isize - cds_start as isize;
    let frame = (v_aln.query_start as isize - delta).rem_euclid(3) as usize;

    let v_start = v_aln.query_start.max(v_aln.query_end.saturating_sub(V_ANCHOR_WINDOW));
    let v_end = v_aln.query_end.min(rearr_end);
    let j_start = j_aln.query_start.min(rearr_end);
    let j_end = j_aln.query_end.min(j_start.saturating_add(J_ANCHOR_WINDOW)).min(rearr_end);

    let mut v_c = None;
    let mut q = frame;
    while q + 3 <= v_end {
        if q >= v_start && codon(observed[q], observed[q + 1], observed[q + 2]) == b'C' {
            v_c = Some(q);
        }
        q += 3;
    }
    let Some(v_c) = v_c else {
        return empty_productivity(ProductivityStatus::MissingVAnchor);
    };

    let mut j_anchor = None;
    let mut q = frame;
    while q + 12 <= j_end {
        if q >= j_start {
            let a0 = codon(observed[q], observed[q + 1], observed[q + 2]);
            let a1 = codon(observed[q + 3], observed[q + 4], observed[q + 5]);
            let a3 = codon(observed[q + 9], observed[q + 10], observed[q + 11]);
            if matches!(a0, b'F' | b'W') && a1 == b'G' && a3 == b'G' {
                j_anchor = Some(q);
                break;
            }
        }
        q += 3;
    }
    let Some(j_anchor) = j_anchor else {
        return empty_productivity(ProductivityStatus::MissingJAnchor);
    };
    if v_c >= j_anchor {
        return empty_productivity(ProductivityStatus::InvalidAnchorOrder);
    }

    let stop_codon = (v_c..rearr_end.saturating_sub(2))
        .step_by(3)
        .any(|q| codon(observed[q], observed[q + 1], observed[q + 2]) == b'*');

    let junction_end = (j_anchor + 3).min(observed.len());
    let junction = observed[v_c..junction_end].to_vec();
    let junction_aa = translate(&junction);
    let cdr3 = if junction.len() >= 6 { junction[3..junction.len() - 3].to_vec() } else { Vec::new() };
    let cdr3_aa = translate(&cdr3);

    ProductivityAssessment {
        productive: !stop_codon,
        in_frame: true,
        stop_codon,
        status: if stop_codon {
            ProductivityStatus::StopCodon
        } else {
            ProductivityStatus::Productive
        },
        junction,
        junction_aa,
        cdr3,
        cdr3_aa,
    }
}

fn translate(seq: &[u8]) -> Vec<u8> {
    seq.chunks_exact(3)
        .map(|x| codon(x[0], x[1], x[2]))
        .collect()
}

fn codon(a: u8, b: u8, c: u8) -> u8 {
    match (a, b, c) {
        (b'T', b'T', b'T' | b'C') => b'F',
        (b'T', b'T', _) => b'L',
        (b'T', b'C', _) => b'S',
        (b'T', b'A', b'T' | b'C') => b'Y',
        (b'T', b'A', _) => b'*',
        (b'T', b'G', b'A') => b'*',
        (b'T', b'G', b'G') => b'W',
        (b'T', b'G', _) => b'C',
        (b'C', b'T', _) => b'L',
        (b'C', b'C', _) => b'P',
        (b'C', b'A', b'T' | b'C') => b'H',
        (b'C', b'A', _) => b'Q',
        (b'C', b'G', _) => b'R',
        (b'A', b'T', b'T' | b'C' | b'A') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'A', b'C', _) => b'T',
        (b'A', b'A', b'T' | b'C') => b'N',
        (b'A', b'A', _) => b'K',
        (b'A', b'G', b'T' | b'C') => b'S',
        (b'A', b'G', _) => b'R',
        (b'G', b'T', _) => b'V',
        (b'G', b'C', _) => b'A',
        (b'G', b'A', b'T' | b'C') => b'D',
        (b'G', b'A', _) => b'E',
        (b'G', b'G', _) => b'G',
        _ => b'X',
    }
}

pub fn validate_recombination(r: &Recombination, index: &VdjIndex) -> Result<()> {
    let v = index
        .segment(r.v)
        .ok_or_else(|| anyhow::anyhow!("missing V segment"))?;
    let j = index
        .segment(r.j)
        .ok_or_else(|| anyhow::anyhow!("missing J segment"))?;
    if v.chain != r.chain
        || v.kind != SegmentKind::V
        || j.chain != r.chain
        || j.kind != SegmentKind::J
    {
        bail!("recombination segment/chain mismatch")
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn aln(query_start: usize, query_end: usize) -> LocalAlignment {
        LocalAlignment {
            score: 1,
            query_start,
            query_end,
            reference_start: 0,
            reference_end: query_end - query_start,
        }
    }

    #[test]
    fn productivity_uses_conserved_anchors_not_sequence_start() {
        // Coding frame starts at nucleotide 1: a one-base 5' UTR prefix must not
        // make an otherwise productive V-J sequence appear out of frame.
        let observed = b"AAAATGTGCCAAATGGGGTAAAGGTGCC";
        let p = assess_productivity(observed, aln(1, 10), aln(13, 28), Some(0));
        assert!(p.productive);
        assert!(p.in_frame);
        assert!(!p.stop_codon);
        assert!(!p.junction_aa.is_empty());
    }

    #[test]
    fn productivity_reports_stop_in_selected_vj_frame() {
        let observed = b"AAAATGTGCCTAATGGGGTAAAGGTGCC";
        let p = assess_productivity(observed, aln(1, 10), aln(13, 28), Some(0));
        assert!(!p.productive);
        assert!(p.in_frame);
        assert!(p.stop_codon);
    }

    #[test]
    fn productivity_rejects_v_and_j_anchors_in_different_frames() {
        // V conserved C is in frame 1, while an inserted nucleotide shifts the
        // J F/W-G-X-G motif into frame 2.
        let observed = b"AAAATGTAAACTGGGGTAAAGGT";
        let p = assess_productivity(observed, aln(1, 10), aln(11, 23), Some(0));
        assert!(!p.productive);
        assert!(!p.in_frame);
        assert!(!p.stop_codon);
    }

    #[test]
    fn productivity_without_v_cds_is_explicitly_unknown() {
        let observed = b"AAAATGTGCCAAATGGGGTAAAGGTGCC";
        let p = assess_productivity(observed, aln(1, 10), aln(13, 28), None);
        assert!(!p.productive);
        assert_eq!(p.status, ProductivityStatus::MissingVCodingStart);
        assert!(p.status.is_unknown());
    }

    #[test]
    fn local_alignment_reports_geometry() {
        let h = local_alignment(b"TTTAACCGGTTGGG", b"AACCGGTT");
        assert_eq!(
            (
                h.score,
                h.query_start,
                h.query_end,
                h.reference_start,
                h.reference_end
            ),
            (16, 3, 11, 0, 8)
        );
    }
}
