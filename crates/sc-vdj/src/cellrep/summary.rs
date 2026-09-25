use super::{BamFeatureEvidence, SequencePart};
use crate::index::{Chain, SegmentId, VdjIndex};
use int_to_dna::IntToDna;
use onehot_dna::OneHotSequence;

/// Compact germline-aware sequence summary for one connected receptor fragment.
///
/// Raw BAM sequence parts are consumed into these summaries immediately.  We
/// retain A/C/G/T counts and maximum qualities per position, plus compact
/// germline segment support/placement metadata.  No read-level sequence history
/// is kept after a successful merge.
#[derive(Debug, Clone, Default, PartialEq)]
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
    /// Packed current best sequence. Kept in lock-step with the A/C/G/T
    /// evidence so hot-path matching and the later rediscovery pass do not
    /// have to rebuild a sequence representation from the count vectors.
    /// `None` means the current consensus contains an unobserved/ambiguous N.
    packed_consensus: Option<IntToDna>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GermlineAnchor {
    pub segment_id: SegmentId,
    pub summary_start: isize,
}

const PACKED_BOOTSTRAP_READS: usize = 50;
const PACKED_STAR_WINDOW_BASES: u32 = 16;
const PACKED_STAR_WINDOW_MIN_READS: usize = 3;

#[derive(Debug, Clone)]
struct PendingPackedRead {
    packed: IntToDna,
    qualities: Vec<u8>,
    segment_ids: Vec<SegmentId>,
    germline_anchors: Vec<GermlineAnchor>,
}

impl PendingPackedRead {
    fn from_part(
        part: &SequencePart,
        segment_ids: &[SegmentId],
        reverse: bool,
        index: &VdjIndex,
        min_overlap: usize,
    ) -> Option<Self> {
        let bases = if reverse {
            crate::index::reverse_complement(&part.bases)
        } else {
            part.bases.clone()
        };
        let packed = IntToDna::try_new(&bases).ok()?;
        let mut qualities = part.qualities.clone();
        if reverse {
            qualities.reverse();
        }
        let mut germline_anchors = Vec::with_capacity(segment_ids.len());
        for &id in segment_ids {
            if let Some(segment) = index.segment(id) {
                let anchor_overlap = min_overlap.min(8).max(4);
                if let Some(read_start) =
                    fast_anchor_offset(&segment.sequence, &bases, anchor_overlap)
                {
                    germline_anchors.push(GermlineAnchor {
                        segment_id: id,
                        summary_start: -read_start,
                    });
                }
            }
        }
        Some(Self {
            packed,
            qualities,
            segment_ids: segment_ids.to_vec(),
            germline_anchors,
        })
    }

    #[inline]
    fn mapping_sequence(&self) -> OneHotSequence {
        OneHotSequence::from_2bit_bytes(&self.packed.u8_encoded, self.packed.size)
    }

    fn into_summary(self, index: &VdjIndex) -> ReceptorSequenceEvidence {
        let n = self.packed.size;
        let mut out = ReceptorSequenceEvidence {
            base_counts: std::array::from_fn(|_| vec![0; n]),
            base_max_qual: std::array::from_fn(|_| vec![0; n]),
            segment_support: Vec::with_capacity(self.segment_ids.len()),
            germline_anchors: self.germline_anchors,
            support_features: 1,
            packed_consensus: Some(self.packed.clone()),
        };
        for pos in 0..n {
            let byte = self.packed.u8_encoded[pos >> 2];
            let base = ((byte >> ((pos & 3) * 2)) & 0b11) as usize;
            out.base_counts[base][pos] = 1;
            out.base_max_qual[base][pos] = self.qualities.get(pos).copied().unwrap_or(0);
        }
        for id in self.segment_ids {
            add_segment_support(&mut out.segment_support, id, 1);
        }
        out.refresh_packed_consensus(index);
        out
    }
}

impl ReceptorSequenceEvidence {
    pub fn len(&self) -> usize {
        self.base_counts.iter().map(Vec::len).max().unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn packed_consensus(&self) -> Option<&IntToDna> {
        self.packed_consensus.as_ref()
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

    /// Packed mapping view of all observed nucleotide states.
    ///
    /// Unlike `consensus()`, this deliberately does not collapse competing
    /// observations. A position with A and G evidence becomes the IUPAC mask R
    /// (A|G), while a position with no observed base becomes a zero mask and is
    /// ignored by overlap scoring. The count/quality evidence remains the
    /// authoritative assembly state; this is only its cheap mapping view.
    pub fn mapping_sequence(&self) -> OneHotSequence {
        let n = self.len();
        let mut masks = Vec::with_capacity(n);
        for pos in 0..n {
            let mut mask = 0u8;
            for base in 0..4 {
                if self.base_counts[base].get(pos).copied().unwrap_or(0) != 0 {
                    mask |= 1u8 << base;
                }
            }
            masks.push(mask);
        }
        OneHotSequence::from_masks(&masks)
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
            packed_consensus: None,
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
        out.refresh_packed_consensus(index);
        out
    }

    #[inline]
    fn refresh_packed_consensus(&mut self, index: &VdjIndex) {
        let consensus = self.consensus(index);
        self.packed_consensus = IntToDna::try_new(&consensus).ok();
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

        // Map against every nucleotide state actually supported by the compact
        // evidence instead of repeatedly materializing a lossy byte consensus.
        // OneHot compatibility is a nibble AND and the overlap scorer consumes
        // 32 positions at a time.
        let a = self.mapping_sequence();
        let b = other.mapping_sequence();
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

        // Mapper/germline geometry has already proposed these offsets. Do not
        // search the alignment space for the common case: verify the known
        // placement directly with packed OneHot comparisons. Up to five
        // incompatible informative positions are accepted here; uglier
        // overlaps fall through to the existing indel-aware rescue.
        let offset = candidates
            .into_iter()
            .find(|off| known_offset_is_compatible_packed(&a, &b, *off, min_overlap, 5))
            .or_else(|| best_offset_packed(&a, &b, min_overlap).map(|x| x.0));
        if let Some(offset) = offset {
            trace_overlap_placement("packed", self, other, index, offset);
            self.merge_at(other, offset);
            self.refresh_packed_consensus(index);
            return true;
        }

        // The cheap mapper above is deliberately ungapped.  A single indel in
        // an otherwise convincing receptor overlap therefore shifts the rest
        // of the read and can incorrectly seed a second summary.  Only after
        // that fast path fails, try an indel-aware overlap alignment.  This is
        // Needleman-Wunsch with free terminal gaps (semi-global / overlap
        // alignment): internal gaps are penalised, while unrelated read ends
        // are not.
        //
        // Keep this rescue conservative: it is only allowed for summaries
        // already linked by a germline segment, and the aligned evidence must
        // still be >=90% identical over a substantial overlap.  The merge is
        // performed through the alignment itself, so an insertion/deletion
        // does not shift downstream evidence into the wrong collector bins.
        if self.shares_segment(other) {
            let ac = self.consensus(index);
            let bc = other.consensus(index);
            if let Some(aln) = needleman_wunsch_overlap(&ac, &bc, min_overlap) {
                self.merge_aligned(other, &aln.columns);
                self.refresh_packed_consensus(index);
                return true;
            }
        }

        false
    }

    fn merge_aligned(&mut self, other: &Self, columns: &[(Option<usize>, Option<usize>)]) {
        // Build a new coordinate system from the gapped alignment.  Columns
        // containing only `self` or only `other` become real collector
        // positions; terminal unaligned sequence is included by the overlap
        // aligner as one-sided columns.
        let n = columns.len();
        let mut counts: [Vec<u16>; 4] = std::array::from_fn(|_| vec![0; n]);
        let mut quals: [Vec<u8>; 4] = std::array::from_fn(|_| vec![0; n]);
        let mut self_to_new = vec![None; self.len()];
        let mut other_to_new = vec![None; other.len()];

        for (dst, &(ai, bi)) in columns.iter().enumerate() {
            if let Some(i) = ai {
                self_to_new[i] = Some(dst);
                for base in 0..4 {
                    counts[base][dst] = counts[base][dst]
                        .saturating_add(self.base_counts[base].get(i).copied().unwrap_or(0));
                    quals[base][dst] =
                        quals[base][dst].max(self.base_max_qual[base].get(i).copied().unwrap_or(0));
                }
            }
            if let Some(i) = bi {
                other_to_new[i] = Some(dst);
                for base in 0..4 {
                    counts[base][dst] = counts[base][dst]
                        .saturating_add(other.base_counts[base].get(i).copied().unwrap_or(0));
                    quals[base][dst] = quals[base][dst]
                        .max(other.base_max_qual[base].get(i).copied().unwrap_or(0));
                }
            }
        }

        self.base_counts = counts;
        self.base_max_qual = quals;

        for anchor in &mut self.germline_anchors {
            if anchor.summary_start >= 0 {
                let old = anchor.summary_start as usize;
                if let Some(Some(new)) = self_to_new.get(old) {
                    anchor.summary_start = *new as isize;
                }
            } else {
                // Preserve the coordinate relation for anchors starting before
                // the observed summary.  The alignment cannot insert columns
                // before self position zero without representing them as
                // one-sided terminal columns, so shift by self(0)'s new start.
                let shift = self_to_new.first().and_then(|x| *x).unwrap_or(0) as isize;
                anchor.summary_start += shift;
            }
        }
        for anchor in &other.germline_anchors {
            let mapped = if anchor.summary_start >= 0 {
                other_to_new
                    .get(anchor.summary_start as usize)
                    .and_then(|x| *x)
                    .map(|x| x as isize)
            } else {
                let shift = other_to_new.first().and_then(|x| *x).unwrap_or(0) as isize;
                Some(anchor.summary_start + shift)
            };
            if let Some(summary_start) = mapped {
                let shifted = GermlineAnchor {
                    segment_id: anchor.segment_id,
                    summary_start,
                };
                if !self.germline_anchors.contains(&shifted) {
                    self.germline_anchors.push(shifted);
                }
            }
        }
        for &(id, n) in &other.segment_support {
            add_segment_support(&mut self.segment_support, id, n);
        }
        self.support_features = self.support_features.saturating_add(other.support_features);
    }

    fn trace_split_dump(&self, label: &str, index: &VdjIndex) {
        let seq = self.consensus(index);
        let segments = self
            .segment_support
            .iter()
            .map(|(id, n)| {
                index
                    .segment(*id)
                    .map(|s| format!("{}:{n}", s.name))
                    .unwrap_or_else(|| format!("#{id:?}:{n}"))
            })
            .collect::<Vec<_>>()
            .join(",");
        eprintln!(
            "[sc-vdj split-trace] {label}: support={} len={} segments=[{}] consensus={}",
            self.support_features,
            self.len(),
            segments,
            String::from_utf8_lossy(&seq),
        );
        // The collector is already represented by the ambiguity-aware consensus
        // above.  Keep split tracing human-readable: the raw A/C/G/T position
        // vectors are useful internally, but obscure the actual sequence being
        // accepted or rejected.
        let anchors = self
            .germline_anchors
            .iter()
            .map(|a| {
                index
                    .segment(a.segment_id)
                    .map(|s| format!("{}@{}", s.name, a.summary_start))
                    .unwrap_or_else(|| format!("#{:?}@{}", a.segment_id, a.summary_start))
            })
            .collect::<Vec<_>>()
            .join(",");
        eprintln!("[sc-vdj split-trace]   anchors: [{anchors}]");
    }

    fn merge_packed_at(
        &mut self,
        other: &PendingPackedRead,
        offset: isize,
        index: &VdjIndex,
    ) {
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
            for anchor in &mut self.germline_anchors {
                anchor.summary_start += prepend as isize;
            }
        }
        let start = if offset < 0 { 0 } else { offset as usize };
        let need = start + other.packed.size;
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
        for pos in 0..other.packed.size {
            let byte = other.packed.u8_encoded[pos >> 2];
            let base = ((byte >> ((pos & 3) * 2)) & 0b11) as usize;
            let dst = start + pos;
            self.base_counts[base][dst] = self.base_counts[base][dst].saturating_add(1);
            self.base_max_qual[base][dst] = self.base_max_qual[base][dst]
                .max(other.qualities.get(pos).copied().unwrap_or(0));
        }
        for &id in &other.segment_ids {
            add_segment_support(&mut self.segment_support, id, 1);
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
        self.support_features = self.support_features.saturating_add(1);
        self.refresh_packed_consensus(index);
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
fn trace_new_summary_split(
    context: &str,
    incoming: &ReceptorSequenceEvidence,
    existing: &[ReceptorSequenceEvidence],
    index: &VdjIndex,
) {
    if std::env::var_os("SC_VDJ_TRACE_SUMMARY_SPLITS").is_none() {
        return;
    }
    eprintln!("\n[sc-vdj split-trace] NEW SUMMARY: {context}; incoming was rejected by {} existing model(s)", existing.len());
    incoming.trace_split_dump("INCOMING", index);
    for (i, model) in existing.iter().enumerate() {
        model.trace_split_dump(&format!("REJECTING MODEL #{i}"), index);
    }
    eprintln!("[sc-vdj split-trace] END NEW SUMMARY\n");
}

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
                trace_new_summary_split("merge_summary_batch", &incoming, summaries, index);
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
    let mut pending = Vec::<PendingPackedRead>::with_capacity(PACKED_BOOTSTRAP_READS);
    let mut pending_star_window: Option<(i32, u32)> = None;

    let flush_pending = |pending: &mut Vec<PendingPackedRead>,
                         summaries: &mut Vec<ReceptorSequenceEvidence>| {
        if pending.is_empty() {
            return;
        }

        // The bootstrap reads have not entered ReceptorSequenceEvidence yet.
        // Seed one collector, then use shared germline geometry to propose an
        // exact coordinate for the remaining packed reads. OneHot verifies only
        // that coordinate; accepted reads splat directly into the collector.
        let first = pending.remove(0);
        let mut collector = first.into_summary(index);
        let mut unresolved = Vec::new();

        for read in pending.drain(..) {
            let a = collector.mapping_sequence();
            let b = read.mapping_sequence();
            let mut offsets = Vec::<isize>::new();
            for aa in &collector.germline_anchors {
                for bb in &read.germline_anchors {
                    if aa.segment_id == bb.segment_id {
                        offsets.push(aa.summary_start - bb.summary_start);
                    }
                }
            }
            offsets.sort_unstable();
            offsets.dedup();

            if let Some(offset) = offsets.into_iter().find(|off| {
                known_offset_is_compatible_packed(&a, &b, *off, min_overlap, 5)
            }) {
                collector.merge_packed_at(&read, offset, index);
            } else {
                unresolved.push(read);
            }
        }

        // Only the assembled stack and genuinely unresolved reads enter the old
        // general summary machinery. Successful packed reads never become
        // individual ReceptorSequenceEvidence objects.
        let mut incoming = Vec::with_capacity(1 + unresolved.len());
        incoming.push(collector);
        incoming.extend(unresolved.into_iter().map(|read| read.into_summary(index)));
        merge_summary_batch(summaries, incoming, index, min_overlap);
    };

    // STAR has already done the expensive genomic placement.  Feed nearby
    // alignments to the packed collector together instead of preserving BAM
    // arrival order: this makes the existing shared-anchor/direct packed check
    // see the reads in genomic order and lets one collector absorb a local pile
    // before persistent summaries are touched.  Unmapped/rescued evidence sorts
    // last and retains the old fallback behaviour.
    let mut ordered_features = features.to_vec();
    ordered_features.sort_unstable_by_key(|feature| {
        feature
            .mappings
            .iter()
            .filter(|m| m.alignment.tid >= 0)
            .map(|m| (m.alignment.tid, m.alignment.start))
            .min()
            .unwrap_or((i32::MAX, u32::MAX))
    });

    for feature in ordered_features {
        let feature_star_start = feature
            .mappings
            .iter()
            .filter(|m| m.alignment.tid >= 0)
            .map(|m| (m.alignment.tid, m.alignment.start))
            .min();

        // Keep the packed collector local in STAR coordinate space.  Sixteen
        // bases is enough to amortize a real placement, but never flush a tiny
        // pile: retain at least three observations so sparse HC evidence is not
        // accidentally reduced to a single anchor.  Every read still enters
        // either the packed collector or the exact old fallback below.
        if let (Some((window_tid, window_start)), Some((tid, start))) =
            (pending_star_window, feature_star_start)
        {
            let outside_window = tid != window_tid
                || start.saturating_sub(window_start) >= PACKED_STAR_WINDOW_BASES;
            if outside_window && pending.len() >= PACKED_STAR_WINDOW_MIN_READS {
                flush_pending(&mut pending, summaries);
                pending_star_window = None;
            }
        }

        if pending_star_window.is_none() {
            pending_star_window = feature_star_start;
        }

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
            if let Some(read) =
                PendingPackedRead::from_part(part, &segment_ids, reverse, index, min_overlap)
            {
                pending.push(read);
                if pending.len() == PACKED_BOOTSTRAP_READS {
                    flush_pending(&mut pending, summaries);
                    pending_star_window = None;
                }
            } else {
                // IUPAC/ambiguous input cannot be represented losslessly in
                // the 2-bit staging buffer. Preserve the old path for it.
                flush_pending(&mut pending, summaries);
                pending_star_window = None;
                let incoming = ReceptorSequenceEvidence::from_part(
                    part,
                    &segment_ids,
                    reverse,
                    index,
                    min_overlap,
                );
                merge_summary_batch(summaries, vec![incoming], index, min_overlap);
            }
        }
    }

    flush_pending(&mut pending, summaries);
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

#[derive(Debug)]
struct OverlapAlignment {
    columns: Vec<(Option<usize>, Option<usize>)>,
}

/// Indel-aware overlap alignment used only as a rescue after the packed,
/// ungapped fast path failed.  This is Needleman-Wunsch with free terminal
/// gaps: choose the best endpoint on the last row/column, then traceback to an
/// edge.  We require >=90% identity among aligned A/C/G/T pairs and at least
/// `max(min_overlap, 24)` informative aligned bases.
fn needleman_wunsch_overlap(a: &[u8], b: &[u8], min_overlap: usize) -> Option<OverlapAlignment> {
    if a.is_empty() || b.is_empty() {
        return None;
    }
    let n = a.len();
    let m = b.len();
    let width = m + 1;
    let mut score = vec![0i32; (n + 1) * width];
    let mut trace = vec![0u8; (n + 1) * width]; // 1 diag, 2 up, 3 left

    // Free leading terminal gaps: first row/column stay zero.
    for i in 1..=n {
        for j in 1..=m {
            let x = a[i - 1].to_ascii_uppercase();
            let y = b[j - 1].to_ascii_uppercase();
            let pair = if x == y && matches!(x, b'A' | b'C' | b'G' | b'T') {
                3
            } else {
                -3
            };
            let diag = score[(i - 1) * width + (j - 1)] + pair;
            let up = score[(i - 1) * width + j] - 4;
            let left = score[i * width + (j - 1)] - 4;
            let (best, dir) = if diag >= up && diag >= left {
                (diag, 1)
            } else if up >= left {
                (up, 2)
            } else {
                (left, 3)
            };
            score[i * width + j] = best;
            trace[i * width + j] = dir;
        }
    }

    // Free trailing terminal gaps: end on whichever last-row/last-column cell
    // has the best score.
    let mut end = (n, m);
    let mut best = i32::MIN;
    for j in 1..=m {
        let s = score[n * width + j];
        if s > best {
            best = s;
            end = (n, j);
        }
    }
    for i in 1..=n {
        let s = score[i * width + m];
        if s > best {
            best = s;
            end = (i, m);
        }
    }

    let (mut i, mut j) = end;
    let mut core = Vec::<(Option<usize>, Option<usize>)>::new();
    let mut informative = 0usize;
    let mut matches = 0usize;
    let mut indel_bases = 0usize;
    while i > 0 && j > 0 {
        match trace[i * width + j] {
            1 => {
                let ai = i - 1;
                let bj = j - 1;
                let x = a[ai].to_ascii_uppercase();
                let y = b[bj].to_ascii_uppercase();
                if matches!(x, b'A' | b'C' | b'G' | b'T') && matches!(y, b'A' | b'C' | b'G' | b'T')
                {
                    informative += 1;
                    if x == y {
                        matches += 1;
                    }
                }
                core.push((Some(ai), Some(bj)));
                i -= 1;
                j -= 1;
            }
            2 => {
                core.push((Some(i - 1), None));
                indel_bases += 1;
                i -= 1;
            }
            3 => {
                core.push((None, Some(j - 1)));
                indel_bases += 1;
                j -= 1;
            }
            _ => break,
        }
    }
    let start_i = i;
    let start_j = j;
    core.reverse();

    let required = min_overlap.max(24);
    let compared = informative + indel_bases;
    if compared < required || matches * 100 < compared * 90 {
        return None;
    }

    let mut columns = Vec::with_capacity(n + m);
    // At most one side has an unaligned leading terminal prefix at traceback
    // termination.  Preserve it as one-sided evidence.
    for ai in 0..start_i {
        columns.push((Some(ai), None));
    }
    for bj in 0..start_j {
        columns.push((None, Some(bj)));
    }
    columns.extend(core);
    for ai in end.0..n {
        columns.push((Some(ai), None));
    }
    for bj in end.1..m {
        columns.push((None, Some(bj)));
    }

    Some(OverlapAlignment { columns })
}

fn trace_overlap_placement(
    method: &str,
    a: &ReceptorSequenceEvidence,
    b: &ReceptorSequenceEvidence,
    index: &VdjIndex,
    offset: isize,
) {
    if std::env::var_os("SC_VDJ_TRACE_OVERLAPS").is_none() {
        return;
    }

    let a_seq = a.consensus(index);
    let b_seq = b.consensus(index);
    let left = 0isize.min(offset);
    let a_pad = (0isize - left) as usize;
    let b_pad = (offset - left) as usize;
    let width = (a_pad + a_seq.len()).max(b_pad + b_seq.len());

    let mut a_line = vec![b' '; width];
    let mut b_line = vec![b' '; width];
    a_line[a_pad..a_pad + a_seq.len()].copy_from_slice(&a_seq);
    b_line[b_pad..b_pad + b_seq.len()].copy_from_slice(&b_seq);

    let mut marks = vec![b' '; width];
    for pos in 0..width {
        let aa = a_line[pos];
        let bb = b_line[pos];
        if aa == b' ' || bb == b' ' {
            continue;
        }
        marks[pos] = if aa == bb && aa != b'N' { b'|' } else { b'.' };
    }

    eprintln!("\\n[sc-vdj overlap-trace] method={method} offset={offset}");
    eprintln!("[sc-vdj overlap-trace] A: {}", String::from_utf8_lossy(&a_seq));
    eprintln!("[sc-vdj overlap-trace] B: {}", String::from_utf8_lossy(&b_seq));
    eprintln!("[sc-vdj overlap-trace] placed:");
    eprintln!("[sc-vdj overlap-trace] A  {}", String::from_utf8_lossy(&a_line));
    eprintln!("[sc-vdj overlap-trace]    {}", String::from_utf8_lossy(&marks));
    eprintln!("[sc-vdj overlap-trace] B  {}", String::from_utf8_lossy(&b_line));
}

fn known_offset_is_compatible_packed(
    a: &OneHotSequence,
    b: &OneHotSequence,
    off: isize,
    min_overlap: usize,
    max_mismatches: usize,
) -> bool {
    let a0 = off.max(0) as usize;
    let b0 = (-off).max(0) as usize;
    if a0 >= a.len() || b0 >= b.len() {
        return false;
    }
    let ov = (a.len() - a0).min(b.len() - b0);
    if ov < min_overlap {
        return false;
    }
    let Some((informative, compatible)) = a.compatibility_counts(a0, b, b0, ov) else {
        return false;
    };
    informative >= min_overlap && informative.saturating_sub(compatible) <= max_mismatches
}

fn overlap_is_compatible_packed(
    a: &OneHotSequence,
    b: &OneHotSequence,
    off: isize,
    min_overlap: usize,
) -> bool {
    let a0 = off.max(0) as usize;
    let b0 = (-off).max(0) as usize;
    if a0 >= a.len() || b0 >= b.len() {
        return false;
    }
    let ov = (a.len() - a0).min(b.len() - b0);
    if ov < min_overlap {
        return false;
    }
    let Some((informative, compatible)) = a.compatibility_counts(a0, b, b0, ov) else {
        return false;
    };
    informative >= min_overlap && compatible * 100 >= informative * 90
}

fn best_offset_packed(
    a: &OneHotSequence,
    b: &OneHotSequence,
    min_overlap: usize,
) -> Option<(isize, usize)> {
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
        let Some((informative, compatible)) = a.compatibility_counts(a0, b, b0, ov) else {
            continue;
        };
        if informative < min_overlap || compatible * 100 < informative * 90 {
            continue;
        }
        let cand = (off, compatible);
        if best.is_none_or(|x: (isize, usize)| cand.1 > x.1) {
            best = Some(cand);
        }
    }
    best
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
