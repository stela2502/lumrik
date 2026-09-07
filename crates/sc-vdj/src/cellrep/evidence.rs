use super::{BamFeatureEvidence, CellEvidence};
use crate::index::{Chain, SegmentId, SegmentKind, Strand, VdjIndex};
use std::collections::{BTreeSet, VecDeque};
use std::fmt;


/// Human-readable view of the compact receptor summaries produced by
/// `summarize_chain_work()`. This is the exact state consumed by
/// `identify_summary()`; no V/D/J/C re-alignment is performed here.
pub struct SummarizedEvidenceDisplay<'a> {
    evidence: &'a CellEvidence,
    index: &'a VdjIndex,
    chain: Chain,
    min_overlap: usize,
}

/// Human-readable, pre-merge view of the raw BAM evidence for one cell/locus.
///
/// Reads are oriented exactly as `summarize_chain_work()` orients them. Their
/// relative placement is inferred from the same germline-anchor and observed
/// sequence-overlap helpers used by the real assembler. Nothing is merged.
pub struct RawEvidenceDisplay<'a> {
    evidence: &'a CellEvidence,
    index: &'a VdjIndex,
    chain: Chain,
    min_overlap: usize,
}

impl CellEvidence {
    pub fn display_raw_chain<'a>(
        &'a self,
        index: &'a VdjIndex,
        chain: Chain,
        min_overlap: usize,
    ) -> RawEvidenceDisplay<'a> {
        RawEvidenceDisplay {
            evidence: self,
            index,
            chain,
            min_overlap,
        }
    }

    pub fn display_summarized_chain<'a>(
        &'a self,
        index: &'a VdjIndex,
        chain: Chain,
        min_overlap: usize,
    ) -> SummarizedEvidenceDisplay<'a> {
        SummarizedEvidenceDisplay {
            evidence: self,
            index,
            chain,
            min_overlap,
        }
    }
}

#[derive(Clone)]
struct RawRead {
    label: String,
    bases: Vec<u8>,
    segment_ids: Vec<SegmentId>,
}

#[derive(Clone, Copy)]
enum NodeKind {
    Segment(SegmentId),
    Read(usize),
}

#[derive(Clone, Copy)]
struct Edge {
    to: usize,
    delta: isize,
}

impl fmt::Display for SummarizedEvidenceDisplay<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let features = self.evidence.features_for_chain(self.index, self.chain);
        let summaries = super::summary::summarize_chain_work(super::summary::ChainSummaryWork {
            features: &features,
            index: self.index,
            chain: self.chain,
            min_overlap: self.min_overlap,
        });

        writeln!(f, "{} POST-MERGE SUMMARIES", self.chain)?;
        if summaries.is_empty() {
            writeln!(f, "  <no summaries>")?;
            return Ok(());
        }

        for (summary_no, summary) in summaries.iter().enumerate() {
            writeln!(f)?;
            writeln!(
                f,
                "summary {}  len={}  supporting_features={}",
                summary_no + 1,
                summary.len(),
                summary.support_features
            )?;

            let consensus = summary.consensus(self.index);
            writeln!(f, "consensus {}", String::from_utf8_lossy(&consensus))?;

            writeln!(f, "segment_support:")?;
            let mut support = summary.segment_support.clone();
            support.sort_by_key(|(id, _)| {
                self.index
                    .segment(*id)
                    .map(|s| (s.kind, s.name.clone()))
            });
            for (id, count) in support {
                if let Some(segment) = self.index.segment(id) {
                    writeln!(
                        f,
                        "  {:?} {}  support={}",
                        segment.kind, segment.name, count
                    )?;
                } else {
                    writeln!(f, "  {:?}  support={}", id, count)?;
                }
            }

            writeln!(f, "germline_anchors:")?;
            if summary.germline_anchors.is_empty() {
                writeln!(f, "  <none>")?;
            } else {
                let mut anchors = summary.germline_anchors.clone();
                anchors.sort_by_key(|a| {
                    self.index
                        .segment(a.segment_id)
                        .map(|s| (s.kind, s.name.clone(), a.summary_start))
                });
                for anchor in anchors {
                    if let Some(segment) = self.index.segment(anchor.segment_id) {
                        writeln!(
                            f,
                            "  {:?} {}  summary_start={}  germline_len={}",
                            segment.kind,
                            segment.name,
                            anchor.summary_start,
                            segment.sequence.len()
                        )?;
                    } else {
                        writeln!(
                            f,
                            "  {:?}  summary_start={}",
                            anchor.segment_id, anchor.summary_start
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}

impl fmt::Display for RawEvidenceDisplay<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let reads = collect_reads(self.evidence, self.index, self.chain);
        writeln!(f, "{} PRE-MERGE RAW EVIDENCE", self.chain)?;

        if reads.is_empty() {
            writeln!(f, "  <no sequence-bearing evidence>")?;
            return Ok(());
        }

        let segment_ids: Vec<_> = reads
            .iter()
            .flat_map(|r| r.segment_ids.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let mut nodes = Vec::<NodeKind>::new();
        for id in &segment_ids {
            nodes.push(NodeKind::Segment(*id));
        }
        for i in 0..reads.len() {
            nodes.push(NodeKind::Read(i));
        }
        let segment_node_count = segment_ids.len();
        let mut graph = vec![Vec::<Edge>::new(); nodes.len()];

        // Read <-> germline anchors: identical placement calculation to the
        // one used when raw reads enter ReceptorSequenceEvidence.
        for (ri, read) in reads.iter().enumerate() {
            let read_node = segment_node_count + ri;
            let anchor_overlap = self.min_overlap.min(8).max(4);
            for &segment_id in &read.segment_ids {
                let Some(segment) = self.index.segment(segment_id) else {
                    continue;
                };
                let Some(off) = super::summary::fast_anchor_offset(
                    &segment.sequence,
                    &read.bases,
                    anchor_overlap,
                ) else {
                    continue;
                };
                let Some(segment_node) = segment_ids.iter().position(|id| *id == segment_id) else {
                    continue;
                };

                // fast_anchor_offset returns reference_start - query_start.
                // Therefore segment_start = read_start - off.
                add_constraint(&mut graph, read_node, segment_node, -off);
            }
        }

        // Raw read <-> read sequence bridges. This is intentionally debug-only
        // O(n^2) work and never participates in production assembly.
        for a in 0..reads.len() {
            for b in (a + 1)..reads.len() {
                if let Some((off, _matches)) = super::summary::best_offset(
                    &reads[a].bases,
                    &reads[b].bases,
                    self.min_overlap,
                ) {
                    add_constraint(
                        &mut graph,
                        segment_node_count + a,
                        segment_node_count + b,
                        off,
                    );
                }
            }
        }

        let components = solve_components(&graph);
        for (component_no, component) in components.iter().enumerate() {
            writeln!(f)?;
            writeln!(f, "component {}", component_no + 1)?;
            writeln!(f, "-----------")?;
            render_component(
                f,
                component,
                &nodes,
                &reads,
                self.index,
                self.chain,
            )?;
        }
        Ok(())
    }
}

fn collect_reads(evidence: &CellEvidence, index: &VdjIndex, chain: Chain) -> Vec<RawRead> {
    let mut out = Vec::new();
    for feature in evidence.features_for_chain(index, chain) {
        let mut segment_ids: Vec<_> = feature
            .mappings
            .iter()
            .filter_map(|m| {
                index
                    .segment(m.segment_id)
                    .filter(|s| s.chain == chain)
                    .map(|_| m.segment_id)
            })
            .collect();
        segment_ids.sort_unstable();
        segment_ids.dedup();
        if segment_ids.is_empty() {
            continue;
        }

        // Keep this orientation rule byte-for-byte equivalent to
        // summarize_chain_work().
        let reverse = feature
            .mappings
            .iter()
            .filter_map(|m| {
                let s = index.segment(m.segment_id)?;
                (s.chain == chain).then_some(
                    m.alignment.is_reverse ^ matches!(s.strand, Strand::Minus),
                )
            })
            .reduce(|a, b| if a == b { a } else { false })
            .unwrap_or(false);

        add_part(&mut out, feature, "R1", feature.sequence.r1.as_ref(), &segment_ids, reverse);
        add_part(&mut out, feature, "R2", feature.sequence.r2.as_ref(), &segment_ids, reverse);
    }
    out
}

fn add_part(
    out: &mut Vec<RawRead>,
    feature: &BamFeatureEvidence,
    mate: &str,
    part: Option<&super::SequencePart>,
    segment_ids: &[SegmentId],
    reverse: bool,
) {
    let Some(part) = part else { return };
    if part.bases.is_empty() {
        return;
    }
    let bases = if reverse {
        crate::index::reverse_complement(&part.bases)
    } else {
        part.bases.clone()
    };
    out.push(RawRead {
        label: format!("read {}:{}/{}", feature.id.flush, feature.id.entry, mate),
        bases,
        segment_ids: segment_ids.to_vec(),
    });
}

fn add_constraint(graph: &mut [Vec<Edge>], a: usize, b: usize, b_minus_a: isize) {
    graph[a].push(Edge {
        to: b,
        delta: b_minus_a,
    });
    graph[b].push(Edge {
        to: a,
        delta: -b_minus_a,
    });
}

/// Return connected components as `(node, solved_start)` pairs.
fn solve_components(graph: &[Vec<Edge>]) -> Vec<Vec<(usize, isize)>> {
    let mut seen = vec![false; graph.len()];
    let mut out = Vec::new();
    for root in 0..graph.len() {
        if seen[root] {
            continue;
        }
        let mut q = VecDeque::new();
        let mut component = Vec::new();
        seen[root] = true;
        q.push_back((root, 0isize));
        while let Some((node, pos)) = q.pop_front() {
            component.push((node, pos));
            for edge in &graph[node] {
                if !seen[edge.to] {
                    seen[edge.to] = true;
                    q.push_back((edge.to, pos + edge.delta));
                }
            }
        }
        out.push(component);
    }
    out
}

fn render_component(
    f: &mut fmt::Formatter<'_>,
    component: &[(usize, isize)],
    nodes: &[NodeKind],
    reads: &[RawRead],
    index: &VdjIndex,
    chain: Chain,
) -> fmt::Result {
    #[derive(Clone)]
    struct RenderRow {
        order: u8,
        start: isize,
        label: String,
        bases: Vec<u8>,
        germline: bool,
    }

    let mut rows = Vec::<RenderRow>::new();
    for &(node, start) in component {
        match nodes[node] {
            NodeKind::Segment(id) => {
                let Some(s) = index.segment(id) else { continue };
                if s.chain != chain {
                    continue;
                }
                rows.push(RenderRow {
                    order: segment_kind_order(s.kind),
                    start,
                    label: format!("germline {}", s.name),
                    bases: s.sequence.clone(),
                    germline: true,
                });
            }
            NodeKind::Read(ri) => {
                let r = &reads[ri];
                rows.push(RenderRow {
                    order: 10,
                    start,
                    label: r.label.clone(),
                    bases: r.bases.clone(),
                    germline: false,
                });
            }
        }
    }
    if rows.is_empty() {
        return Ok(());
    }

    let min_start = rows.iter().map(|r| r.start).min().unwrap_or(0);
    let max_end = rows
        .iter()
        .map(|r| r.start + r.bases.len() as isize)
        .max()
        .unwrap_or(min_start + 1);
    let label_width = rows.iter().map(|r| r.label.len()).max().unwrap_or(0).min(32);

    writeln!(f, "bases {}..{}  (1 bp/column)", min_start, max_end)?;

    // Sort first, then borrow germline sequences from the stable row storage.
    // This avoids cloning germlines purely for rendering while also avoiding
    // the borrow conflict caused by sorting after taking references into rows.
    rows.sort_by(|a, b| (a.order, a.start, &a.label).cmp(&(b.order, b.start, &b.label)));

    // Keep germline sequence only as borrowed comparison masks. For a read,
    // '-' means that its base agrees with at least one germline covering that
    // canvas position. Literal A/C/G/T therefore highlights either a mismatch
    // or sequence outside all germline tracks (for example junction sequence).
    let germlines: Vec<(isize, &[u8])> = rows
        .iter()
        .filter(|r| r.germline)
        .map(|r| (r.start, r.bases.as_slice()))
        .collect();

    for row in &rows {
        let left = (row.start - min_start).max(0) as usize;

        write!(
            f,
            "{label:label_width$}  {padding}",
            label = row.label,
            padding = " ".repeat(left),
        )?;

        if row.germline {
            for _ in 0..row.bases.len() {
                write!(f, "=")?;
            }
        } else {
            // Allocate only one byte per read base, not one byte per complete
            // component column. A very wide/disconnected layout therefore does
            // not multiply memory use by the number of reads.
            let mut rendered = Vec::with_capacity(row.bases.len());
            for (i, &base) in row.bases.iter().enumerate() {
                let absolute = row.start + i as isize;
                let base = base.to_ascii_uppercase();

                let matches_germline = germlines.iter().any(|(germline_start, germline)| {
                    let offset = absolute - *germline_start;
                    if offset < 0 || offset >= germline.len() as isize {
                        return false;
                    }
                    germline[offset as usize].to_ascii_uppercase() == base
                });

                rendered.push(if matches_germline { b'-' } else { base });
            }
            write!(f, "{}", String::from_utf8_lossy(&rendered))?;
        }

        writeln!(
            f,
            "  [{start}..{}]",
            row.start + row.bases.len() as isize,
            start = row.start,
        )?;
    }
    Ok(())
}

fn segment_kind_order(kind: SegmentKind) -> u8 {
    match kind {
        SegmentKind::V => 0,
        SegmentKind::D => 1,
        SegmentKind::J => 2,
        SegmentKind::C => 3,
    }
}
