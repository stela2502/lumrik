use super::{BamFeatureEvidence, SequencePart};
use crate::index::{Chain, SegmentId, VdjIndex};

/// Compact germline-aware sequence summary for one connected receptor fragment.
///
/// Raw BAM sequence parts are consumed into these summaries immediately.  We
/// retain A/C/G/T counts and maximum qualities per position, plus compact
/// germline segment support/placement metadata.  No read-level sequence history
/// is kept after a successful merge.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReceptorSequenceEvidence {
    pub base_counts: [Vec<u16>; 4],
    pub base_max_qual: [Vec<u8>; 4],
    /// (segment id, number of contributing BAM features nominating it).
    pub segment_support: Vec<(SegmentId, u32)>,
    /// Placement of germline segment coordinate 0 in summary coordinates.
    /// Missing placements are simply omitted; the segment can still be used as
    /// a merge/indexing identity and for downstream support scoring.
    pub germline_anchors: Vec<GermlineAnchor>,
    pub support_features: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GermlineAnchor {
    pub segment_id: SegmentId,
    pub summary_start: isize,
}

impl ReceptorSequenceEvidence {
    pub fn len(&self) -> usize {
        self.base_counts.iter().map(Vec::len).max().unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn segment_ids(&self) -> impl Iterator<Item = SegmentId> + '_ {
        self.segment_support.iter().map(|(id, _)| *id)
    }

    pub fn segment_support(&self, id: SegmentId) -> u32 {
        self.segment_support
            .iter()
            .find_map(|(x, n)| (*x == id).then_some(*n))
            .unwrap_or(0)
    }

    /// Resolve a consensus only on demand.  Counts win first, maximum base
    /// quality breaks count ties, and an observed tied base matching a germline
    /// anchor breaks a remaining tie.  Germline can never introduce a base that
    /// was not observed.
    pub fn consensus(&self, index: &VdjIndex) -> Vec<u8> {
        let n = self.len();
        let mut out = Vec::with_capacity(n);
        for pos in 0..n {
            let mut max_count = 0u16;
            for b in 0..4 {
                max_count = max_count.max(self.base_counts[b].get(pos).copied().unwrap_or(0));
            }
            if max_count == 0 {
                out.push(b'N');
                continue;
            }

            let mut max_q = 0u8;
            for b in 0..4 {
                let c = self.base_counts[b].get(pos).copied().unwrap_or(0);
                if c == max_count {
                    max_q = max_q.max(self.base_max_qual[b].get(pos).copied().unwrap_or(0));
                }
            }

            let mut tied = [false; 4];
            let mut tied_count = 0usize;
            for b in 0..4 {
                let c = self.base_counts[b].get(pos).copied().unwrap_or(0);
                let q = self.base_max_qual[b].get(pos).copied().unwrap_or(0);
                if c == max_count && q == max_q {
                    tied[b] = true;
                    tied_count += 1;
                }
            }

            let mut chosen = None;
            if tied_count > 1 {
                // Require agreeing germline votes.  A disagreeing/irrelevant
                // germline does not get to manufacture a nucleotide.
                let mut germline_vote = None::<usize>;
                let mut conflicting_vote = false;
                for anchor in &self.germline_anchors {
                    let Some(segment) = index.segment(anchor.segment_id) else {
                        continue;
                    };
                    let gp = pos as isize - anchor.summary_start;
                    if gp < 0 || gp as usize >= segment.sequence.len() {
                        continue;
                    }
                    let Some(b) = base_index(segment.sequence[gp as usize]) else {
                        continue;
                    };
                    if !tied[b] {
                        continue;
                    }
                    match germline_vote {
                        None => germline_vote = Some(b),
                        Some(prev) if prev == b => {}
                        Some(_) => conflicting_vote = true,
                    }
                }
                if !conflicting_vote {
                    chosen = germline_vote;
                }
            }

            let b = chosen.unwrap_or_else(|| {
                // Deterministic fallback: A,C,G,T order among still-tied
                // observed alternatives.
                (0..4).find(|b| tied[*b]).unwrap_or(0)
            });
            out.push(base_from_index(b));
        }
        out
    }

    fn from_part(
        part: &SequencePart,
        segment_ids: &[SegmentId],
        reverse: bool,
        index: &VdjIndex,
        min_overlap: usize,
    ) -> Self {
        let bases = if reverse {
            crate::index::reverse_complement(&part.bases)
        } else {
            part.bases.clone()
        };
        let mut qualities = part.qualities.clone();
        if reverse {
            qualities.reverse();
        }

        let mut out = Self {
            base_counts: std::array::from_fn(|_| vec![0; bases.len()]),
            base_max_qual: std::array::from_fn(|_| vec![0; bases.len()]),
            segment_support: Vec::with_capacity(segment_ids.len()),
            germline_anchors: Vec::with_capacity(segment_ids.len()),
            support_features: 1,
        };
        for (i, &base) in bases.iter().enumerate() {
            let Some(k) = base_index(base) else { continue };
            out.base_counts[k][i] = 1;
            out.base_max_qual[k][i] = qualities.get(i).copied().unwrap_or(0);
        }
        for &id in segment_ids {
            add_segment_support(&mut out.segment_support, id, 1);
            if let Some(segment) = index.segment(id) {
                // This is computed once when the read enters the summary model.
                // A short D may not have enough sequence for a reliable anchor;
                // in that case its identity is retained without a placement.
                let anchor_overlap = min_overlap.min(8).max(4);
                if let Some(read_start) =
                    fast_anchor_offset(&segment.sequence, &bases, anchor_overlap)
                {
                    out.germline_anchors.push(GermlineAnchor {
                        segment_id: id,
                        summary_start: -read_start,
                    });
                }
            }
        }
        out
    }

    fn shares_segment(&self, other: &Self) -> bool {
        self.segment_ids().any(|id| other.segment_support(id) > 0)
    }

    fn has_hard_segment_conflict(&self, other: &Self, index: &VdjIndex) -> bool {
        for kind in [
            crate::index::SegmentKind::V,
            crate::index::SegmentKind::D,
            crate::index::SegmentKind::J,
            crate::index::SegmentKind::C,
        ] {
            let a: Vec<_> = self
                .segment_ids()
                .filter(|id| index.segment(*id).is_some_and(|s| s.kind == kind))
                .collect();
            let b: Vec<_> = other
                .segment_ids()
                .filter(|id| index.segment(*id).is_some_and(|s| s.kind == kind))
                .collect();
            if !a.is_empty() && !b.is_empty() && !a.iter().any(|x| b.contains(x)) {
                return true;
            }
        }
        false
    }

    fn try_consume(&mut self, other: &Self, index: &VdjIndex, min_overlap: usize) -> bool {
        if self.has_hard_segment_conflict(other, index) {
            return false;
        }

        let a = self.consensus(index);
        let b = other.consensus(index);
        if a.is_empty() || b.is_empty() {
            return false;
        }

        // Shared germline anchors narrow the expected offset first.  If no
        // placed shared anchor is available (e.g. short D), fall back to an
        // observed-sequence overlap.  This fallback only runs between compact
        // summaries, never against the original read collection.
        let mut candidates = Vec::<isize>::new();
        for aa in &self.germline_anchors {
            for bb in &other.germline_anchors {
                if aa.segment_id == bb.segment_id {
                    candidates.push(aa.summary_start - bb.summary_start);
                }
            }
        }
        candidates.sort_unstable();
        candidates.dedup();

        let offset = candidates
            .into_iter()
            .find(|off| overlap_is_compatible(&a, &b, *off, min_overlap))
            .or_else(|| best_offset(&a, &b, min_overlap).map(|x| x.0));
        let Some(offset) = offset else {
            return false;
        };

        self.merge_at(other, offset);
        true
    }

    fn merge_at(&mut self, other: &Self, offset: isize) {
        let prepend = (-offset).max(0) as usize;
        if prepend > 0 {
            for v in &mut self.base_counts {
                let mut x = vec![0; prepend];
                x.extend_from_slice(v);
                *v = x;
            }
            for v in &mut self.base_max_qual {
                let mut x = vec![0; prepend];
                x.extend_from_slice(v);
                *v = x;
            }
            for a in &mut self.germline_anchors {
                a.summary_start += prepend as isize;
            }
        }
        let start = if offset < 0 { 0 } else { offset as usize };
        let need = start + other.len();
        for v in &mut self.base_counts {
            if v.len() < need {
                v.resize(need, 0);
            }
        }
        for v in &mut self.base_max_qual {
            if v.len() < need {
                v.resize(need, 0);
            }
        }

        for b in 0..4 {
            for i in 0..other.base_counts[b].len() {
                let dst = start + i;
                self.base_counts[b][dst] =
                    self.base_counts[b][dst].saturating_add(other.base_counts[b][i]);
                self.base_max_qual[b][dst] =
                    self.base_max_qual[b][dst].max(other.base_max_qual[b][i]);
            }
        }
        for &(id, n) in &other.segment_support {
            add_segment_support(&mut self.segment_support, id, n);
        }
        for anchor in &other.germline_anchors {
            let shifted = GermlineAnchor {
                segment_id: anchor.segment_id,
                summary_start: anchor.summary_start + start as isize,
            };
            if !self.germline_anchors.contains(&shifted) {
                self.germline_anchors.push(shifted);
            }
        }
        self.support_features = self.support_features.saturating_add(other.support_features);
    }
}

/// One independent unit of expensive receptor work.  The caller intentionally
/// keeps this as a normal function boundary so a future runner can collect
/// work items in a Vec and feed them to Rayon without changing the assembly
/// implementation.
#[derive(Clone, Copy)]
pub struct ChainSummaryWork<'a> {
    pub features: &'a [&'a BamFeatureEvidence],
    pub index: &'a VdjIndex,
    pub chain: Chain,
    pub min_overlap: usize,
}

/// Consume all sequence-bearing BAM evidence for one cell/locus into a small
/// set of germline-aware summaries.  Reads disappear as soon as they fit an
/// existing summary.  Compatible but unlinked summaries remain separate until
/// actual sequence/germline evidence bridges them.
pub fn summarize_chain_work(work: ChainSummaryWork<'_>) -> Vec<ReceptorSequenceEvidence> {
    let mut summaries = Vec::<ReceptorSequenceEvidence>::new();
    consume_chain_features(
        &mut summaries,
        work.features,
        work.index,
        work.chain,
        work.min_overlap,
    );
    summaries
}

/// Fold another bounded feature batch into an already-persistent chain summary.
/// The caller retains only `summaries`; all raw feature objects may be dropped
/// immediately after this returns.
/// Merge the compact delta produced by one bounded BAM batch into the
/// persistent summaries for one cell/chain.
///
/// The expensive raw-read assembly has already happened against an empty,
/// batch-local summary set.  Persistent state is therefore probed only once
/// per compact incoming summary rather than once per BAM sequence part.
pub(crate) fn merge_summary_batch(
    summaries: &mut Vec<ReceptorSequenceEvidence>,
    batch_summaries: Vec<ReceptorSequenceEvidence>,
    index: &VdjIndex,
    min_overlap: usize,
) {
    for incoming in batch_summaries {
        let mut target = None;

        // Shared germline identity is the normal path and is both cheaper and
        // biologically more constrained than a sequence-only search.
        for i in 0..summaries.len() {
            if !summaries[i].shares_segment(&incoming) {
                continue;
            }
            if summaries[i].try_consume(&incoming, index, min_overlap) {
                target = Some(i);
                break;
            }
        }

        // A genuinely bridging batch summary can connect components that do
        // not yet share a segment identity.  Keep the old observed-overlap
        // fallback, but pay for it once per compact batch summary, not once per
        // original read.
        if target.is_none() {
            for i in 0..summaries.len() {
                if summaries[i].shares_segment(&incoming) {
                    continue;
                }
                if summaries[i].try_consume(&incoming, index, min_overlap) {
                    target = Some(i);
                    break;
                }
            }
        }

        let target = match target {
            Some(i) => i,
            None => {
                summaries.push(incoming);
                summaries.len() - 1
            }
        };

        // One incoming compact component may bridge two persistent components.
        // Reconcile only around the component that changed.
        collapse_from(summaries, target, index, min_overlap);
    }

    summaries.sort_by_key(|s| std::cmp::Reverse(s.support_features));
}

pub(crate) fn consume_chain_features(
    summaries: &mut Vec<ReceptorSequenceEvidence>,
    features: &[&BamFeatureEvidence],
    index: &VdjIndex,
    chain: Chain,
    min_overlap: usize,
) {
    for feature in features {
        let mut segment_ids: Vec<_> = feature
            .mappings
            .iter()
            .filter_map(|m| {
                index
                    .segment(m.segment_id)
                    .filter(|s| s.chain == chain && s.kind != crate::index::SegmentKind::C)
                    .map(|_| m.segment_id)
            })
            .collect();
        segment_ids.sort_unstable();
        segment_ids.dedup();
        if segment_ids.is_empty() {
            continue;
        }

        // The mapper geometry tells us how read orientation relates to the
        // transcript-oriented germline sequence.  If mappings disagree, keep
        // forward orientation and let observed overlap decide later.
        let reverse = feature
            .mappings
            .iter()
            .filter_map(|m| {
                let s = index.segment(m.segment_id)?;
                (s.chain == chain).then_some(
                    m.alignment.is_reverse ^ matches!(s.strand, crate::index::Strand::Minus),
                )
            })
            .reduce(|a, b| if a == b { a } else { false })
            .unwrap_or(false);

        for part in [&feature.sequence.r1, &feature.sequence.r2]
            .into_iter()
            .flatten()
        {
            if part.bases.is_empty() {
                continue;
            }
            let incoming = ReceptorSequenceEvidence::from_part(
                part,
                &segment_ids,
                reverse,
                index,
                min_overlap,
            );

            let mut target = None;

            // First try summaries sharing a germline identity.  This is the
            // common fast path and avoids global sequence searching.
            for i in 0..summaries.len() {
                if !summaries[i].shares_segment(&incoming) {
                    continue;
                }
                if summaries[i].try_consume(&incoming, index, min_overlap) {
                    target = Some(i);
                    break;
                }
            }

            // A bridging fragment may connect different germline identities
            // through observed overlap (V->D, D->J, J->C).  This fallback is
            // only over compact summaries, never over retained BAM fragments.
            if target.is_none() {
                for i in 0..summaries.len() {
                    if summaries[i].shares_segment(&incoming) {
                        continue; // already tested above
                    }
                    if summaries[i].try_consume(&incoming, index, min_overlap) {
                        target = Some(i);
                        break;
                    }
                }
            }

            let target = match target {
                Some(i) => i,
                None => {
                    summaries.push(incoming);
                    summaries.len() - 1
                }
            };

            // Only the summary touched by this new evidence can have acquired a
            // new bridge to another existing component.  Re-test that one
            // component rather than rescanning every pair after every read.
            collapse_from(summaries, target, index, min_overlap);
        }
    }

    summaries.sort_by_key(|s| std::cmp::Reverse(s.support_features));
}

fn collapse_from(
    summaries: &mut Vec<ReceptorSequenceEvidence>,
    mut target: usize,
    index: &VdjIndex,
    min_overlap: usize,
) {
    loop {
        let mut merged = false;
        for other in 0..summaries.len() {
            if other == target {
                continue;
            }
            if summaries[target].has_hard_segment_conflict(&summaries[other], index) {
                continue;
            }

            let candidate = summaries[other].clone();
            if !summaries[target].try_consume(&candidate, index, min_overlap) {
                continue;
            }

            summaries.swap_remove(other);
            if other < summaries.len() && target == summaries.len() {
                // target was the last element and got moved into `other` by
                // swap_remove.
                target = other;
            }
            merged = true;
            break;
        }
        if !merged {
            break;
        }
    }
}

fn add_segment_support(dst: &mut Vec<(SegmentId, u32)>, id: SegmentId, n: u32) {
    if let Some((_, count)) = dst.iter_mut().find(|(x, _)| *x == id) {
        *count = count.saturating_add(n);
    } else {
        dst.push((id, n));
    }
}

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn base_from_index(i: usize) -> u8 {
    [b'A', b'C', b'G', b'T'][i.min(3)]
}

fn overlap_is_compatible(a: &[u8], b: &[u8], off: isize, min_overlap: usize) -> bool {
    let a0 = off.max(0) as usize;
    let b0 = (-off).max(0) as usize;
    if a0 >= a.len() || b0 >= b.len() {
        return false;
    }
    let ov = (a.len() - a0).min(b.len() - b0);
    if ov < min_overlap {
        return false;
    }
    let mut informative = 0usize;
    let mut matches = 0usize;
    for i in 0..ov {
        let x = a[a0 + i];
        let y = b[b0 + i];
        if x == b'N' || y == b'N' {
            continue;
        }
        informative += 1;
        if x == y {
            matches += 1;
        }
    }
    informative >= min_overlap && matches * 100 >= informative * 90
}

pub(super) fn fast_anchor_offset(reference: &[u8], query: &[u8], seed_len: usize) -> Option<isize> {
    let k = seed_len.min(reference.len()).min(query.len());
    if k < 4 {
        return None;
    }

    // Probe a few deterministic seeds across the observed sequence.  This
    // avoids a full dynamic-programming alignment for every BAM fragment while
    // still tolerating junction sequence outside the anchored germline region.
    let max_q = query.len() - k;
    let mut starts = [0usize, max_q / 3, (max_q * 2) / 3, max_q];
    starts.sort_unstable();
    let mut best = None::<(usize, isize)>;

    for &qs in &starts {
        let seed = &query[qs..qs + k];
        for rs in 0..=reference.len() - k {
            if &reference[rs..rs + k] != seed {
                continue;
            }
            let off = rs as isize - qs as isize;
            let r0 = off.max(0) as usize;
            let q0 = (-off).max(0) as usize;
            let ov = (reference.len() - r0).min(query.len() - q0);
            if ov < k {
                continue;
            }
            let matches = (0..ov)
                .filter(|i| reference[r0 + i] == query[q0 + i])
                .count();
            if matches * 100 < ov * 80 {
                continue;
            }
            if best.is_none_or(|x| matches > x.0) {
                best = Some((matches, off));
            }
        }
    }
    best.map(|(_, off)| off)
}

/// Find a high-confidence observed overlap.  This remains as a bridge fallback
/// for components carrying different germline identities, but it is no longer
/// run as an all-reads-vs-all-reads greedy assembler.
pub(super) fn best_offset(a: &[u8], b: &[u8], min_overlap: usize) -> Option<(isize, usize)> {
    if a.len() < min_overlap || b.len() < min_overlap {
        return None;
    }
    let mut best = None;
    let lo = -(b.len() as isize) + min_overlap as isize;
    let hi = a.len() as isize - min_overlap as isize;
    for off in lo..=hi {
        let a0 = off.max(0) as usize;
        let b0 = (-off).max(0) as usize;
        let ov = (a.len() - a0).min(b.len() - b0);
        if ov < min_overlap {
            continue;
        }
        let mut informative = 0usize;
        let mut matches = 0usize;
        for i in 0..ov {
            let x = a[a0 + i];
            let y = b[b0 + i];
            if x == b'N' || y == b'N' {
                continue;
            }
            informative += 1;
            if x == y {
                matches += 1;
            }
        }
        if informative < min_overlap || matches * 100 < informative * 90 {
            continue;
        }
        let cand = (off, matches);
        if best.is_none_or(|x: (isize, usize)| cand.1 > x.1) {
            best = Some(cand);
        }
    }
    best
}
