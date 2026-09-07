use crate::cellrep::{summarize_chain_work, CellEvidence, ChainSummaryWork, ReceptorSequenceEvidence};
use crate::index::{
    reference_base_matches, reverse_complement, Chain, SegmentId, SegmentKind, VdjIndex,
};
use anyhow::{bail, Result};

mod identifier;
pub use identifier::{ReceptorRole, RecombinationId};

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
    pub supporting_features: u32,
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

    let features = work.cell.features_for_chain(work.index, work.chain);
    let summaries = summarize_chain_work(ChainSummaryWork {
        features: &features,
        index: work.index,
        chain: work.chain,
        min_overlap: work.min_overlap,
    });

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
    let productivity = assess_productivity(&observed, v_aln, j_aln);
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
        supporting_features: assembled.support_features,
        stable_id: RecombinationId::placeholder(chain),
    };
    recomb.stable_id = RecombinationId::from_recombination(&recomb, index).ok()?;
    Some(recomb)
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ProductivityAssessment {
    productive: bool,
    in_frame: bool,
    stop_codon: bool,
}

/// Determine receptor productivity from the observed V->J coding core.
///
/// The reconstructed sequence can start in 5' UTR (or in the middle of V), so
/// sequence position zero is not a coding-frame origin.  Instead, evaluate all
/// three possible frames and keep the one that contains both:
///
/// * a conserved V cysteine close to the 3' end of the V alignment, and
/// * a J F/W-G-X-G motif close to the 5' end of the J alignment.
///
/// A stop codon is then assessed only in that selected coding frame from the
/// conserved V cysteine through the observed J alignment.
fn assess_productivity(
    observed: &[u8],
    v_aln: LocalAlignment,
    j_aln: LocalAlignment,
) -> ProductivityAssessment {
    const V_ANCHOR_WINDOW: usize = 90;
    const J_ANCHOR_WINDOW: usize = 120;

    let rearr_end = j_aln.query_end.min(observed.len());
    if rearr_end < 3 {
        return ProductivityAssessment {
            productive: false,
            in_frame: false,
            stop_codon: false,
        };
    }

    let v_start = v_aln
        .query_start
        .max(v_aln.query_end.saturating_sub(V_ANCHOR_WINDOW));
    let v_end = v_aln.query_end.min(rearr_end);
    let j_start = j_aln.query_start.min(rearr_end);
    let j_end = j_aln
        .query_end
        .min(j_start.saturating_add(J_ANCHOR_WINDOW))
        .min(rearr_end);

    // (distance from expected V/J ends, frame, V-C nucleotide position,
    //  J-motif nucleotide position).  Smaller distance is better.
    let mut best: Option<(usize, usize, usize, usize)> = None;

    for frame in 0..3usize {
        let mut v_c = None;
        let mut q = frame;
        while q + 3 <= v_end {
            if q >= v_start && codon(observed[q], observed[q + 1], observed[q + 2]) == b'C' {
                // The conserved V cysteine is expected toward the 3' end.
                v_c = Some(q);
            }
            q += 3;
        }
        let Some(v_c) = v_c else {
            continue;
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
            continue;
        };

        if v_c >= j_anchor {
            continue;
        }

        let v_distance = v_end.saturating_sub(v_c + 3);
        let j_distance = j_anchor.saturating_sub(j_start);
        let candidate = (v_distance + j_distance, frame, v_c, j_anchor);
        if best.map_or(true, |current| candidate < current) {
            best = Some(candidate);
        }
    }

    let Some((_, frame, v_c, _j_anchor)) = best else {
        return ProductivityAssessment {
            productive: false,
            in_frame: false,
            stop_codon: false,
        };
    };

    debug_assert_eq!(v_c % 3, frame);
    let stop_codon = (v_c..rearr_end.saturating_sub(2))
        .step_by(3)
        .any(|q| codon(observed[q], observed[q + 1], observed[q + 2]) == b'*');

    ProductivityAssessment {
        productive: !stop_codon,
        in_frame: true,
        stop_codon,
    }
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
        let p = assess_productivity(observed, aln(1, 10), aln(13, 28));
        assert_eq!(
            p,
            ProductivityAssessment {
                productive: true,
                in_frame: true,
                stop_codon: false,
            }
        );
    }

    #[test]
    fn productivity_reports_stop_in_selected_vj_frame() {
        let observed = b"AAAATGTGCCTAATGGGGTAAAGGTGCC";
        let p = assess_productivity(observed, aln(1, 10), aln(13, 28));
        assert_eq!(
            p,
            ProductivityAssessment {
                productive: false,
                in_frame: true,
                stop_codon: true,
            }
        );
    }

    #[test]
    fn productivity_rejects_v_and_j_anchors_in_different_frames() {
        // V conserved C is in frame 1, while an inserted nucleotide shifts the
        // J F/W-G-X-G motif into frame 2.
        let observed = b"AAAATGTAAACTGGGGTAAAGGT";
        let p = assess_productivity(observed, aln(1, 10), aln(11, 23));
        assert!(!p.productive);
        assert!(!p.in_frame);
        assert!(!p.stop_codon);
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
