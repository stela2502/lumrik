use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use sc_vdj::{Chain, Recombination, SegmentKind, VdjRunner};
use std::fmt::Write as _;
use std::path::Path;

#[derive(Debug, Default)]
struct BamRecordCounts {
    total: usize,
    mapped: usize,
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
        .map(|(_, cell)| cell.accepted_records())
        .sum();
    let accepted_fragments: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| cell.physical_fragments())
        .sum();
    let single_segment_records: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| cell.single_segment_records())
        .sum();
    let multi_segment_records: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| cell.multi_segment_records())
        .sum();
    let linked_fragments: usize = runner
        .evidence
        .cells()
        .map(|(_, cell)| cell.linked_fragments())
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
    writeln!(
        out,
        "cells with VDJ evidence:     {}",
        runner.evidence.cell_count()
    )?;
    writeln!(
        out,
        "raw read-level evidence:     discarded after each 20,000-record batch"
    )?;

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

    writeln!(out, "  {chain}")?;
    writeln!(
        out,
        "    evidence: {} records / compact fragment summaries",
        cell.chain_records(chain)
    )?;
    writeln!(
        out,
        "    segment mappings: V={} D={} J={} C={}",
        cell.segment_mappings(chain, SegmentKind::V),
        cell.segment_mappings(chain, SegmentKind::D),
        cell.segment_mappings(chain, SegmentKind::J),
        cell.segment_mappings(chain, SegmentKind::C),
    )?;

    let mut links: Vec<_> = cell
        .fragment_link_support()
        .iter()
        .filter(|(signature, _)| signature.chain == chain)
        .collect();
    links.sort_by_key(|(signature, _)| {
        (
            signature.receptor_segments.clone(),
            signature.constant_segments.clone(),
        )
    });
    if links.is_empty() {
        writeln!(out, "    fragment-level receptor -> C links: none")?;
    } else {
        writeln!(out, "    fragment-level receptor -> C links:")?;
        for (signature, count) in links {
            let receptor = signature
                .receptor_segments
                .iter()
                .filter_map(|id| runner.index.segment(*id).map(|s| s.name.as_str()))
                .collect::<Vec<_>>()
                .join(",");
            let constant = signature
                .constant_segments
                .iter()
                .filter_map(|id| runner.index.segment(*id).map(|s| s.name.as_str()))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(out, "      {count:>4} x {receptor} -> {constant}")?;
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
            call.stable_id, call.supporting_features
        )?;
    }

    Ok(())
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
    writeln!(
        out,
        "BAM records scanned:        {}",
        report.bam_records_scanned
    )?;
    writeln!(
        out,
        "records from wanted cells:  {}",
        report.wanted_cell_records
    )?;
    writeln!(
        out,
        "CDR3/J bait-hit records:    {}",
        report.receptor_hit_records
    )?;
    writeln!(
        out,
        "constant-region hit records:{}",
        report.constant_hit_records
    )?;
    writeln!(
        out,
        "linked physical fragments:  {}",
        report.linked_fragments
    )?;
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
