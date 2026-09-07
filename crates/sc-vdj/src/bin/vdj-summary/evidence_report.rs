use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use sc_vdj::{Chain, EvidenceId, Recombination, SegmentKind, VdjRunner};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt::Write as _;
use std::path::Path;

#[derive(Debug, Default)]
struct BamRecordCounts {
    total: usize,
    mapped: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct LinkSignature {
    v: Vec<String>,
    d: Vec<String>,
    j: Vec<String>,
    c: Vec<String>,
}

impl LinkSignature {
    fn from_segments(segments: &BTreeSet<(SegmentKind, String)>) -> Self {
        let mut out = Self {
            v: Vec::new(),
            d: Vec::new(),
            j: Vec::new(),
            c: Vec::new(),
        };
        for (kind, name) in segments {
            match kind {
                SegmentKind::V => out.v.push(name.clone()),
                SegmentKind::D => out.d.push(name.clone()),
                SegmentKind::J => out.j.push(name.clone()),
                SegmentKind::C => out.c.push(name.clone()),
            }
        }
        out
    }

    fn has_link(&self) -> bool {
        [self.v.len(), self.d.len(), self.j.len(), self.c.len()]
            .into_iter()
            .filter(|n| *n > 0)
            .count()
            >= 2
    }

    fn render(&self) -> String {
        fn names(xs: &[String]) -> String {
            if xs.is_empty() {
                "-".to_string()
            } else {
                xs.join(",")
            }
        }
        format!(
            "V={} -> D={} -> J={} -> C={}",
            names(&self.v),
            names(&self.d),
            names(&self.j),
            names(&self.c)
        )
    }
}

#[derive(Debug, Clone)]
struct FragmentLinkAudit {
    signature: LinkSignature,
    record_signatures: Vec<LinkSignature>,
}

impl FragmentLinkAudit {
    fn crosses_records(&self) -> bool {
        !self.record_signatures.iter().any(|record| {
            (!record.v.is_empty() || !record.j.is_empty()) && !record.c.is_empty()
        })
    }
}

fn evidence_id_label(id: EvidenceId) -> String {
    format!("{}:{}", id.flush, id.entry)
}

pub fn describe_vdj_run<P: AsRef<Path>>(
    bam_path: P,
    runner: &VdjRunner,
    calls: &[(u64, Vec<Recombination>)],
) -> Result<String> {
    let bam_counts = count_bam_records(bam_path.as_ref())?;
    let mut out = String::new();

    let accepted_records: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| cell.features.len())
        .sum();
    let accepted_fragments: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| {
            cell.features
                .iter()
                .map(|feature| feature.id)
                .collect::<HashSet<_>>()
                .len()
        })
        .sum();
    let single_segment_records: usize = runner
        .evidence
        .cells()
        .flat_map(|(_, cell)| &cell.features)
        .filter(|feature| feature.mappings.len() == 1)
        .count();
    let multi_segment_records: usize = runner
        .evidence
        .cells()
        .flat_map(|(_, cell)| &cell.features)
        .filter(|feature| feature.mappings.len() > 1)
        .count();
    let linked_fragments: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| {
            let mut fragments = BTreeMap::<EvidenceId, BTreeSet<SegmentKind>>::new();
            for feature in &cell.features {
                let kinds = fragments.entry(feature.id).or_default();
                for mapping in &feature.mappings {
                    if let Some(segment) = runner.index.segment(mapping.segment_id) {
                        kinds.insert(segment.kind);
                    }
                }
            }
            fragments.values().filter(|kinds| kinds.len() > 1).count()
        })
        .sum();

    writeln!(out, "VDJ EVIDENCE SUMMARY")?;
    writeln!(out, "====================")?;
    writeln!(out, "BAM records:                 {}", bam_counts.total)?;
    writeln!(out, "mapped BAM records:          {}", bam_counts.mapped)?;
    writeln!(out, "VDJ-overlapping records:     {accepted_records}")?;
    writeln!(out, "physical read fragments:     {accepted_fragments}")?;
    writeln!(out, "single-segment BAM records:  {single_segment_records}")?;
    writeln!(out, "multi-segment BAM records:   {multi_segment_records}")?;
    writeln!(out, "multi-segment fragments:     {linked_fragments}")?;
    writeln!(out, "cells with VDJ evidence:     {}", runner.evidence.cell_count())?;

    for (cell_id, cell) in runner.evidence.cells() {
        let cell_name = runner
            .cell_names
            .get(&cell_id)
            .map(String::as_str)
            .unwrap_or("unknown");
        writeln!(out)?;
        writeln!(out, "CELL {cell_name}")?;
        writeln!(out, "{}", "-".repeat(5 + cell_name.len()))?;

        let cell_calls = calls
            .iter()
            .find(|(id, _)| *id == cell_id)
            .map(|(_, calls)| calls.as_slice())
            .unwrap_or(&[]);

        for chain in cell.chains(&runner.index) {
            describe_chain(&mut out, runner, cell_id, chain, cell_calls)?;
        }
    }

    Ok(out)
}

fn describe_chain(
    out: &mut String,
    runner: &VdjRunner,
    cell_id: u64,
    chain: Chain,
    calls: &[Recombination],
) -> Result<()> {
    let cell = runner.evidence.get(&cell_id).expect("cell disappeared");
    let features = cell.features_for_chain(&runner.index, chain);
    let mut fragment_segments =
        BTreeMap::<EvidenceId, BTreeSet<(SegmentKind, String)>>::new();
    let mut fragment_records = BTreeMap::<EvidenceId, Vec<LinkSignature>>::new();
    let mut record_links = BTreeMap::<LinkSignature, usize>::new();
    let mut kind_records = BTreeMap::<SegmentKind, usize>::new();

    for feature in &features {
        let mut record_segments = BTreeSet::<(SegmentKind, String)>::new();
        for mapping in &feature.mappings {
            let Some(segment) = runner.index.segment(mapping.segment_id) else {
                continue;
            };
            if segment.chain != chain {
                continue;
            }
            record_segments.insert((segment.kind, segment.name.clone()));
            fragment_segments
                .entry(feature.id)
                .or_default()
                .insert((segment.kind, segment.name.clone()));
            *kind_records.entry(segment.kind).or_insert(0) += 1;
        }

        let record_signature = LinkSignature::from_segments(&record_segments);
        if record_signature.has_link() {
            *record_links.entry(record_signature.clone()).or_insert(0) += 1;
        }
        fragment_records
            .entry(feature.id)
            .or_default()
            .push(record_signature);
    }

    let mut fragments = BTreeMap::<EvidenceId, FragmentLinkAudit>::new();
    let mut fragment_link_counts = BTreeMap::<LinkSignature, usize>::new();
    for (id, segments) in &fragment_segments {
        let signature = LinkSignature::from_segments(segments);
        if signature.has_link() {
            *fragment_link_counts.entry(signature.clone()).or_insert(0) += 1;
        }
        fragments.insert(
            *id,
            FragmentLinkAudit {
                signature,
                record_signatures: fragment_records.remove(id).unwrap_or_default(),
            },
        );
    }

    writeln!(out, "  {chain}")?;
    writeln!(
        out,
        "    evidence: {} records / {} physical fragments",
        features.len(),
        fragment_segments.len()
    )?;
    writeln!(
        out,
        "    segment mappings: V={} D={} J={} C={}",
        kind_records.get(&SegmentKind::V).copied().unwrap_or(0),
        kind_records.get(&SegmentKind::D).copied().unwrap_or(0),
        kind_records.get(&SegmentKind::J).copied().unwrap_or(0),
        kind_records.get(&SegmentKind::C).copied().unwrap_or(0),
    )?;

    if record_links.is_empty() {
        writeln!(out, "    direct links within one BAM record: none")?;
    } else {
        writeln!(out, "    direct links within one BAM record:")?;
        for (signature, count) in record_links {
            writeln!(out, "      {count:>4} x {}", signature.render())?;
        }
    }

    if fragment_link_counts.is_empty() {
        writeln!(out, "    fragment-level links (same EvidenceId): none")?;
    } else {
        writeln!(out, "    fragment-level links (same EvidenceId):")?;
        for (signature, count) in fragment_link_counts {
            writeln!(out, "      {count:>4} x {}", signature.render())?;
        }
    }

    let chain_calls: Vec<_> = calls.iter().filter(|call| call.chain == chain).collect();
    if chain_calls.is_empty() {
        writeln!(out, "    recombination calls: none")?;
        return Ok(());
    }

    writeln!(out, "    recombination calls:")?;
    for call in chain_calls {
        let v = runner
            .index
            .segment(call.v)
            .map(|segment| segment.name.as_str())
            .unwrap_or("?");
        let d = call
            .d
            .and_then(|id| runner.index.segment(id))
            .map(|segment| segment.name.as_str())
            .unwrap_or("-");
        let j = runner
            .index
            .segment(call.j)
            .map(|segment| segment.name.as_str())
            .unwrap_or("?");
        let c = call
            .constant
            .as_ref()
            .and_then(|evidence| runner.index.segment(evidence.segment))
            .map(|segment| segment.name.as_str())
            .unwrap_or("-");
        writeln!(
            out,
            "      {}: V={v} D={d} J={j} C={c} support={}",
            call.stable_id,
            call.supporting_features
        )?;

        let direct_c = direct_constant_links(&fragments, v, j);
        if direct_c.is_empty() {
            writeln!(out, "        V/J -> C fragment links used by rescue: none")?;
        } else {
            writeln!(out, "        V/J -> C fragment links used by rescue:")?;
            for link in direct_c {
                writeln!(
                    out,
                    "          fragment {} [{}]: {}",
                    evidence_id_label(link.id),
                    if link.crosses_records {
                        "across BAM records"
                    } else {
                        "within one BAM record"
                    },
                    link.signature.render()
                )?;
            }
        }
    }

    Ok(())
}

fn direct_constant_links(
    fragments: &BTreeMap<EvidenceId, FragmentLinkAudit>,
    called_v: &str,
    called_j: &str,
) -> Vec<ConstantLinkAudit> {
    let mut out = Vec::new();
    for (id, fragment) in fragments {
        let supports_call = fragment.signature.v.iter().any(|name| name == called_v)
            || fragment.signature.j.iter().any(|name| name == called_j);
        if !supports_call || fragment.signature.c.is_empty() {
            continue;
        }
        out.push(ConstantLinkAudit {
            id: *id,
            signature: fragment.signature.clone(),
            crosses_records: fragment.crosses_records(),
        });
    }
    out
}

#[derive(Debug, Clone)]
struct ConstantLinkAudit {
    id: EvidenceId,
    signature: LinkSignature,
    crosses_records: bool,
}

fn count_bam_records(path: &Path) -> Result<BamRecordCounts> {
    let mut reader = bam::Reader::from_path(path)
        .with_context(|| format!("opening {} for VDJ evidence summary", path.display()))?;
    let mut counts = BamRecordCounts::default();
    for record in reader.records() {
        let record = record?;
        counts.total += 1;
        if !record.is_unmapped() {
            counts.mapped += 1;
        }
    }
    Ok(counts)
}

pub fn describe_recombination_rescan(
    runner: &VdjRunner,
    before: &[(u64, Vec<Recombination>)],
    after: &[(u64, Vec<Recombination>)],
    report: &sc_vdj::recombination::RecombinationEvidenceRescanReport,
) -> Result<String> {
    let mut out = String::new();
    writeln!(out, "RECOMBINATION EVIDENCE RESCAN")?;
    writeln!(out, "===========================")?;
    writeln!(out, "BAM records scanned:        {}", report.bam_records_scanned)?;
    writeln!(out, "records from wanted cells:  {}", report.wanted_cell_records)?;
    writeln!(out, "CDR3/J bait-hit records:    {}", report.receptor_hit_records)?;
    writeln!(out, "constant-region hit records:{}", report.constant_hit_records)?;
    writeln!(out, "linked physical fragments:  {}", report.linked_fragments)?;
    writeln!(out, "rescued constant calls:     {}", report.rescued)?;

    for call_report in &report.calls {
        let cell_name = runner
            .cell_names
            .get(&call_report.cell_id)
            .map(String::as_str)
            .unwrap_or("unknown");
        let before_call = find_call(before, call_report.cell_id, &call_report.stable_id);
        let after_call = find_call(after, call_report.cell_id, &call_report.stable_id);

        writeln!(out)?;
        writeln!(out, "CELL {cell_name}")?;
        writeln!(out, "  recombination {}", call_report.stable_id)?;
        writeln!(
            out,
            "  CDR3/J bait: {}",
            styled_bait(runner, before_call.or(after_call), &call_report.bait)
        )?;
        if let Some(call) = before_call {
            writeln!(out, "  initial call: {}", render_call(runner, call))?;
        }
        if call_report.candidates.is_empty() {
            writeln!(out, "  constant candidates: none")?;
        } else {
            writeln!(out, "  constant candidates:")?;
            for candidate in &call_report.candidates {
                let name = runner
                    .index
                    .segment(candidate.segment)
                    .map(|segment| segment.name.as_str())
                    .unwrap_or("unknown");
                writeln!(
                    out,
                    "    {name}: fragments={} same-read={}",
                    candidate.fragments, candidate.same_read
                )?;
            }
        }
        if let Some(segment) = call_report.rescued_segment {
            let name = runner
                .index
                .segment(segment)
                .map(|segment| segment.name.as_str())
                .unwrap_or("unknown");
            writeln!(out, "  rescued C: {name}")?;
        } else {
            writeln!(out, "  rescued C: none")?;
        }
        if let Some(call) = after_call {
            writeln!(out, "  final call:   {}", render_call(runner, call))?;
        }
    }

    Ok(out)
}

fn styled_bait(runner: &VdjRunner, call: Option<&Recombination>, bait: &[u8]) -> String {
    let mut rendered = bait
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    let Some(call) = call else {
        return String::from_utf8_lossy(&rendered).into_owned();
    };

    // The bait is cut from the reconstructed receptor, so use the same called
    // germline V/D/J segments to mark which bait positions are germline-derived.
    // Junction positions not covered by any called germline alignment remain upper-case.
    for segment_id in [Some(call.v), call.d, Some(call.j)].into_iter().flatten() {
        let Some(segment) = runner.index.segment(segment_id) else {
            continue;
        };
        let alignment = sc_vdj::recombination::local_alignment(bait, &segment.sequence);
        if alignment.score <= 0 {
            continue;
        }
        let start = alignment.query_start.min(rendered.len());
        let end = alignment.query_end.min(rendered.len());
        for base in &mut rendered[start..end] {
            *base = base.to_ascii_lowercase();
        }
    }

    String::from_utf8_lossy(&rendered).into_owned()
}

fn find_call<'a>(
    calls: &'a [(u64, Vec<Recombination>)],
    cell_id: u64,
    stable_id: &sc_vdj::RecombinationId,
) -> Option<&'a Recombination> {
    calls
        .iter()
        .find(|(id, _)| *id == cell_id)?
        .1
        .iter()
        .find(|call| &call.stable_id == stable_id)
}

fn render_call(runner: &VdjRunner, call: &Recombination) -> String {
    let v = runner
        .index
        .segment(call.v)
        .map(|segment| segment.name.as_str())
        .unwrap_or("?");
    let d = call
        .d
        .and_then(|id| runner.index.segment(id))
        .map(|segment| segment.name.as_str())
        .unwrap_or("-");
    let j = runner
        .index
        .segment(call.j)
        .map(|segment| segment.name.as_str())
        .unwrap_or("?");
    let c = call
        .constant
        .as_ref()
        .and_then(|constant| runner.index.segment(constant.segment))
        .map(|segment| segment.name.as_str())
        .unwrap_or("-");
    format!("V={v} D={d} J={j} C={c}")
}
