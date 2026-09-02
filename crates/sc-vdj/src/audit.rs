use crate::align::local_alignment;
use crate::posterior::{CellVdjSummary, RearrangementCall, RearrangementSupportingRead};
use crate::reference::VdjReference;
use crate::sequence::reverse_complement;
use crate::types::SegmentKind;
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;

const MIN_CANDIDATE_LOCAL_SCORE: i32 = 18;
const MIN_OVERLAP: usize = 20;
const MAX_VJ_ALIGNMENT_OVERLAP: usize = 12;
const MAX_V_END_DISTANCE: usize = 35;
const MAX_J_START_DISTANCE: usize = 35;

#[derive(Debug, Clone)]
struct OrientedRead {
    sequence: Vec<u8>,
    supports_v: bool,
    supports_j: bool,
    reads_used: usize,
}

#[derive(Debug, Clone)]
struct CandidateSequence {
    sequence: Vec<u8>,
    method: &'static str,
    reads_used: usize,
}

pub struct AuditFastaWriter<'a> {
    reference: &'a VdjReference,
    germline: BufWriter<File>,
    supporting: BufWriter<File>,
    candidates: BufWriter<File>,
    proteins: BufWriter<File>,
}

impl<'a> AuditFastaWriter<'a> {
    pub fn create<P: AsRef<Path>>(dir: P, reference: &'a VdjReference) -> Result<Self> {
        let dir = dir.as_ref();
        fs::create_dir_all(dir)?;
        Ok(Self {
            reference,
            germline: BufWriter::new(
                File::create(dir.join("vdj_germline_segments.fasta"))
                    .context("creating vdj_germline_segments.fasta")?,
            ),
            supporting: BufWriter::new(
                File::create(dir.join("vdj_supporting_reads.fasta"))
                    .context("creating vdj_supporting_reads.fasta")?,
            ),
            candidates: BufWriter::new(
                File::create(dir.join("vdj_rearrangement_candidates.fasta"))
                    .context("creating vdj_rearrangement_candidates.fasta")?,
            ),
            proteins: BufWriter::new(
                File::create(dir.join("vdj_rearrangement_candidate_proteins.fasta"))
                    .context("creating vdj_rearrangement_candidate_proteins.fasta")?,
            ),
        })
    }

    pub fn write_cells(&mut self, cells: &[CellVdjSummary]) -> Result<()> {
        for cell in cells {
            for call in &cell.rearrangements {
                write_germline_call(&mut self.germline, &cell.cell, call, self.reference)?;
                write_supporting_call(&mut self.supporting, &cell.cell, call)?;
                if let Some(candidate) = reconstruct_candidate(call, self.reference) {
                    let v = support_name(call, SegmentKind::V);
                    let d = support_name(call, SegmentKind::D);
                    let j = support_name(call, SegmentKind::J);
                    let c = support_name(call, SegmentKind::C);
                    let header = format!(
                        "cell={}|chain={}|V={}|D={}|J={}|C={}|support_umis={}|method={}|reads={}",
                        header_value(&cell.cell),
                        call.chain,
                        header_value(v),
                        header_value(d),
                        header_value(j),
                        header_value(c),
                        call.total_supporting_umis,
                        candidate.method,
                        candidate.reads_used
                    );
                    writeln!(self.candidates, ">{header}")?;
                    write_fasta_sequence(&mut self.candidates, &candidate.sequence)?;
                    write_candidate_proteins(&mut self.proteins, &header, &candidate.sequence)?;
                }
            }
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.germline.flush()?;
        self.supporting.flush()?;
        self.candidates.flush()?;
        self.proteins.flush()?;
        Ok(())
    }
}

pub fn write_vdj_audit_fastas<P: AsRef<Path>>(
    dir: P,
    cells: &[CellVdjSummary],
    reference: &VdjReference,
) -> Result<()> {
    let mut writer = AuditFastaWriter::create(dir, reference)?;
    writer.write_cells(cells)?;
    writer.finish()
}

fn write_germline_call<W: Write>(
    out: &mut W,
    cell: &str,
    call: &RearrangementCall,
    reference: &VdjReference,
) -> Result<()> {
    for (kind, support) in [
        (SegmentKind::V, call.v.as_ref()),
        (SegmentKind::D, call.d.as_ref()),
        (SegmentKind::J, call.j.as_ref()),
        (SegmentKind::C, call.c.as_ref()),
    ] {
        let Some(support) = support else {
            continue;
        };
        let Some(segment) = reference.segments.get(support.segment_index) else {
            continue;
        };
        writeln!(
            out,
            ">{}|{}|{:?}|{}",
            header_value(cell),
            call.chain,
            kind,
            header_value(&segment.name)
        )?;
        write_fasta_sequence(out, &segment.sequence)?;
    }
    Ok(())
}

fn write_supporting_call<W: Write>(out: &mut W, cell: &str, call: &RearrangementCall) -> Result<()> {
    let v = support_name(call, SegmentKind::V);
    let j = support_name(call, SegmentKind::J);
    for read in &call.supporting_reads {
        writeln!(
            out,
            ">cell={}|chain={}|umi={}|read={}|V={}|J={}|bam_reverse={}|supplementary={}|supports={}|V_align={}|J_align={}",
            header_value(cell),
            call.chain,
            header_value(&read.umi),
            header_value(&read.read_name),
            header_value(v),
            header_value(j),
            read.bam_is_reverse,
            read.is_supplementary,
            support_labels(read),
            alignment_label(read.v_alignment.as_ref()),
            alignment_label(read.j_alignment.as_ref())
        )?;
        write_fasta_sequence(out, &read.sequence)?;
    }
    Ok(())
}

fn reconstruct_candidate(
    call: &RearrangementCall,
    reference: &VdjReference,
) -> Option<CandidateSequence> {
    let v = reference.segments.get(call.v.as_ref()?.segment_index)?;
    let j = reference.segments.get(call.j.as_ref()?.segment_index)?;

    let mut reads = orient_and_deduplicate_reads(&call.supporting_reads, &v.sequence, &j.sequence);
    if reads.is_empty() {
        return None;
    }

    reads.sort_by(|a, b| {
        b.sequence
            .len()
            .cmp(&a.sequence.len())
            .then_with(|| a.sequence.cmp(&b.sequence))
    });

    // Prefer a literal observed read that already spans a geometrically valid V->J
    // junction as the assembly anchor.  Unlike the previous implementation, do not
    // stop there: every other coherent supporting read may extend that anchor at
    // either end through observed sequence overlap.
    if let Some(anchor_index) = reads.iter().position(|read| {
        read.supports_v
            && read.supports_j
            && valid_candidate_geometry(&read.sequence, &v.sequence, &j.sequence)
    }) {
        let anchor = reads.remove(anchor_index);
        let extended = extend_observed_contig(anchor, reads);
        return Some(CandidateSequence {
            sequence: extended.sequence,
            method: if extended.reads_used == 1 {
                "single_read"
            } else {
                "extended_consensus"
            },
            reads_used: extended.reads_used,
        });
    }

    // If no individual read spans V->J, assemble only through observed overlaps.
    // Once a geometrically valid V/J contig appears, continue extending it rather
    // than returning the first (usually shortest) bridge.
    let mut contigs = reads;
    loop {
        let mut best: Option<(usize, usize, usize, usize, Vec<u8>)> = None;
        for i in 0..contigs.len() {
            for j_idx in 0..contigs.len() {
                if i == j_idx {
                    continue;
                }
                let Some((overlap, mismatches, merged)) =
                    overlap_merge(&contigs[i].sequence, &contigs[j_idx].sequence)
                else {
                    continue;
                };
                let replace = best.as_ref().map_or(true, |old| {
                    overlap > old.2
                        || (overlap == old.2
                            && (mismatches < old.3
                                || (mismatches == old.3 && (i, j_idx) < (old.0, old.1))))
                });
                if replace {
                    best = Some((i, j_idx, overlap, mismatches, merged));
                }
            }
        }

        let Some((i, j_idx, _, _, sequence)) = best else {
            break;
        };
        let supports_v = contigs[i].supports_v || contigs[j_idx].supports_v;
        let supports_j = contigs[i].supports_j || contigs[j_idx].supports_j;
        let reads_used = contigs[i].reads_used + contigs[j_idx].reads_used;
        let merged = OrientedRead {
            sequence,
            supports_v,
            supports_j,
            reads_used,
        };

        let hi = i.max(j_idx);
        let lo = i.min(j_idx);
        contigs.remove(hi);
        contigs.remove(lo);
        contigs.push(merged);

        let candidate_index = contigs.iter().position(|contig| {
            contig.supports_v
                && contig.supports_j
                && valid_candidate_geometry(&contig.sequence, &v.sequence, &j.sequence)
        });
        if let Some(candidate_index) = candidate_index {
            let anchor = contigs.remove(candidate_index);
            let extended = extend_observed_contig(anchor, contigs);
            return Some(CandidateSequence {
                sequence: extended.sequence,
                method: "overlap_consensus",
                reads_used: extended.reads_used,
            });
        }
    }

    None
}

fn extend_observed_contig(mut contig: OrientedRead, mut reads: Vec<OrientedRead>) -> OrientedRead {
    loop {
        let mut best: Option<(usize, usize, usize, usize, usize, Vec<u8>)> = None;

        for (idx, read) in reads.iter().enumerate() {
            // read -> contig extends the 5' end; contig -> read extends the 3' end.
            for (direction, left, right) in [
                (0usize, read.sequence.as_slice(), contig.sequence.as_slice()),
                (1usize, contig.sequence.as_slice(), read.sequence.as_slice()),
            ] {
                let Some((overlap, mismatches, merged)) = overlap_merge(left, right) else {
                    continue;
                };
                if merged.len() <= contig.sequence.len() {
                    continue;
                }
                let extension = merged.len() - contig.sequence.len();
                let replace = best.as_ref().map_or(true, |old| {
                    extension > old.1
                        || (extension == old.1
                            && (overlap > old.2
                                || (overlap == old.2
                                    && (mismatches < old.3
                                        || (mismatches == old.3
                                            && (idx, direction) < (old.0, old.4))))))
                });
                if replace {
                    best = Some((idx, extension, overlap, mismatches, direction, merged));
                }
            }
        }

        let Some((idx, _, _, _, _, sequence)) = best else {
            break;
        };
        let read = reads.remove(idx);
        contig.sequence = sequence;
        contig.supports_v |= read.supports_v;
        contig.supports_j |= read.supports_j;
        contig.reads_used += read.reads_used;
    }

    contig
}

fn orient_and_deduplicate_reads(
    reads: &[RearrangementSupportingRead],
    _v: &[u8],
    _j: &[u8],
) -> Vec<OrientedRead> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for read in reads {
        let orientation = match (read.v_alignment.as_ref(), read.j_alignment.as_ref()) {
            (Some(v), Some(j)) if v.reverse_complement == j.reverse_complement => {
                Some(v.reverse_complement)
            }
            (Some(v), None) => Some(v.reverse_complement),
            (None, Some(j)) => Some(j.reverse_complement),
            // A read whose V and J assignments demand opposite orientations cannot
            // establish a transcript-contiguous V->J candidate.
            (Some(_), Some(_)) => None,
            (None, None) => None,
        };
        let Some(reverse) = orientation else { continue };
        let sequence = if reverse {
            reverse_complement(&read.sequence)
        } else {
            read.sequence.clone()
        };
        let key = (sequence.clone(), read.supports_v, read.supports_j);
        if !seen.insert(key) {
            continue;
        }
        out.push(OrientedRead {
            sequence,
            supports_v: read.supports_v,
            supports_j: read.supports_j,
            reads_used: 1,
        });
    }
    out
}

fn overlap_merge(left: &[u8], right: &[u8]) -> Option<(usize, usize, Vec<u8>)> {
    let max_overlap = left.len().min(right.len());
    if max_overlap < MIN_OVERLAP {
        return None;
    }
    for overlap in (MIN_OVERLAP..=max_overlap).rev() {
        let a = &left[left.len() - overlap..];
        let b = &right[..overlap];
        let mismatches = a.iter().zip(b).filter(|(x, y)| x != y).count();
        let allowed = (overlap / 20).max(1);
        if mismatches > allowed {
            continue;
        }
        let mut merged = Vec::with_capacity(left.len() + right.len() - overlap);
        merged.extend_from_slice(&left[..left.len() - overlap]);
        for (&x, &y) in a.iter().zip(b) {
            merged.push(if x == y { x } else { b'N' });
        }
        merged.extend_from_slice(&right[overlap..]);
        return Some((overlap, mismatches, merged));
    }
    None
}

fn valid_candidate_geometry(sequence: &[u8], v: &[u8], j: &[u8]) -> bool {
    let va = local_alignment(sequence, v);
    let ja = local_alignment(sequence, j);
    if va.score < MIN_CANDIDATE_LOCAL_SCORE || ja.score < MIN_CANDIDATE_LOCAL_SCORE {
        return false;
    }
    if v.len().saturating_sub(va.reference_end) > MAX_V_END_DISTANCE
        || ja.reference_start > MAX_J_START_DISTANCE
        || va.query_start > ja.query_start
    {
        return false;
    }
    va.query_end.saturating_sub(ja.query_start) <= MAX_VJ_ALIGNMENT_OVERLAP
}

fn alignment_label(alignment: Option<&crate::posterior::ReadSegmentAlignment>) -> String {
    let Some(a) = alignment else { return "-".to_string() };
    format!(
        "score:{};q:{}-{};ref:{}-{};rc:{}",
        a.score,
        a.query_start,
        a.query_end,
        a.reference_start,
        a.reference_end,
        a.reverse_complement
    )
}

fn write_candidate_proteins<W: Write>(out: &mut W, header: &str, sequence: &[u8]) -> Result<()> {
    for frame in 0..3 {
        let protein = translate_frame(sequence, frame);
        let stops = protein.iter().filter(|&&aa| aa == b'*').count();
        writeln!(out, ">{header}|frame={frame}|aa_len={}|stops={stops}", protein.len())?;
        write_fasta_sequence(out, &protein)?;
    }
    Ok(())
}

fn translate_frame(sequence: &[u8], frame: usize) -> Vec<u8> {
    sequence
        .get(frame..)
        .unwrap_or_default()
        .chunks_exact(3)
        .map(translate_codon)
        .collect()
}

fn translate_codon(codon: &[u8]) -> u8 {
    if codon.len() != 3 {
        return b'X';
    }
    let a = codon[0].to_ascii_uppercase();
    let b = codon[1].to_ascii_uppercase();
    let c = codon[2].to_ascii_uppercase();
    match (a, b, c) {
        (b'T', b'T', b'T' | b'C') => b'F',
        (b'T', b'T', b'A' | b'G') => b'L',
        (b'T', b'C', _) => b'S',
        (b'T', b'A', b'T' | b'C') => b'Y',
        (b'T', b'A', b'A' | b'G') => b'*',
        (b'T', b'G', b'T' | b'C') => b'C',
        (b'T', b'G', b'A') => b'*',
        (b'T', b'G', b'G') => b'W',
        (b'C', b'T', _) => b'L',
        (b'C', b'C', _) => b'P',
        (b'C', b'A', b'T' | b'C') => b'H',
        (b'C', b'A', b'A' | b'G') => b'Q',
        (b'C', b'G', _) => b'R',
        (b'A', b'T', b'T' | b'C' | b'A') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'A', b'C', _) => b'T',
        (b'A', b'A', b'T' | b'C') => b'N',
        (b'A', b'A', b'A' | b'G') => b'K',
        (b'A', b'G', b'T' | b'C') => b'S',
        (b'A', b'G', b'A' | b'G') => b'R',
        (b'G', b'T', _) => b'V',
        (b'G', b'C', _) => b'A',
        (b'G', b'A', b'T' | b'C') => b'D',
        (b'G', b'A', b'A' | b'G') => b'E',
        (b'G', b'G', _) => b'G',
        _ => b'X',
    }
}

fn support_name(call: &RearrangementCall, kind: SegmentKind) -> &str {
    let support = match kind {
        SegmentKind::V => call.v.as_ref(),
        SegmentKind::D => call.d.as_ref(),
        SegmentKind::J => call.j.as_ref(),
        SegmentKind::C => call.c.as_ref(),
    };
    support.map_or("?", |x| x.id.as_str())
}

fn support_labels(read: &RearrangementSupportingRead) -> String {
    let mut labels = Vec::new();
    if read.supports_v {
        labels.push("V");
    }
    if read.supports_d {
        labels.push("D");
    }
    if read.supports_j {
        labels.push("J");
    }
    if read.supports_c {
        labels.push("C");
    }
    labels.join(",")
}

fn header_value(value: &str) -> String {
    value
        .chars()
        .map(|c| if c.is_whitespace() || c == '|' || c == '>' { '_' } else { c })
        .collect()
}

fn write_fasta_sequence<W: Write>(out: &mut W, sequence: &[u8]) -> Result<()> {
    for line in sequence.chunks(80) {
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::posterior::{GermlineSegmentSupport, RecombinationStage};
    use crate::types::{Chain, VdjSegment};

    fn support(index: usize, id: &str, kind: SegmentKind) -> GermlineSegmentSupport {
        GermlineSegmentSupport {
            segment_index: index,
            id: id.to_string(),
            kind,
            local_alignment_score: 100,
            supporting_umis: 1,
            supporting_reads: 1,
            locus_fraction: 0.0,
            distance_to_recombination_center: 0,
        }
    }

    fn segment(name: &str, kind: SegmentKind, sequence: &[u8]) -> VdjSegment {
        VdjSegment {
            name: name.to_string(),
            transcript_id: format!("{name}_tx"),
            gene_id: name.to_string(),
            chain: Chain::Igk,
            kind,
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
    fn exact_overlap_merges_observed_sequence_without_inventing_gap() {
        let left = b"AAAACCCCGGGGTTTTAAAACCCCGGGG";
        let right = b"GGGGTTTTAAAACCCCGGGGTTTTAAAA";
        let (_, _, merged) = overlap_merge(left, right).expect("overlap");
        assert_eq!(merged, b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA");
    }

    #[test]
    fn spanning_anchor_is_extended_on_both_observed_ends() {
        let anchor = OrientedRead {
            sequence: b"CCCCGGGGTTTTAAAACCCCGGGG".to_vec(),
            supports_v: true,
            supports_j: true,
            reads_used: 1,
        };
        let reads = vec![
            OrientedRead {
                sequence: b"AAAACCCCGGGGTTTTAAAACCCC".to_vec(),
                supports_v: true,
                supports_j: false,
                reads_used: 1,
            },
            OrientedRead {
                sequence: b"GGGGTTTTAAAACCCCGGGGTTTT".to_vec(),
                supports_v: false,
                supports_j: true,
                reads_used: 1,
            },
        ];

        let extended = extend_observed_contig(anchor, reads);
        assert_eq!(
            extended.sequence,
            b"AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"
        );
        assert_eq!(extended.reads_used, 3);
        assert!(extended.supports_v);
        assert!(extended.supports_j);
    }

    #[test]
    fn no_overlap_means_no_rearrangement_candidate() {
        let reference = VdjReference {
            segments: vec![
                segment("IGKV1", SegmentKind::V, b"AACCGGTTAACCGGTTAACCGGTT"),
                segment("IGKJ1", SegmentKind::J, b"TTGGAACCTTGGAACCTTGGAACC"),
            ],
        };
        let call = RearrangementCall {
            chain: Chain::Igk,
            stage: RecombinationStage::Vj,
            v: Some(support(0, "IGKV1", SegmentKind::V)),
            d: None,
            j: Some(support(1, "IGKJ1", SegmentKind::J)),
            c: None,
            total_supporting_umis: 1,
            supporting_reads: vec![
                RearrangementSupportingRead {
                    umi: "u1".into(),
                    read_name: "vread".into(),
                    sequence: b"AACCGGTTAACCGGTTAACCGGTT".to_vec(),
                    bam_is_reverse: false,
                    is_supplementary: false,
                    supports_v: true,
                    supports_j: false,
                    supports_d: false,
                    supports_c: false,
                    v_alignment: None,
                    j_alignment: None,
                    d_alignment: None,
                    c_alignment: None,
                },
                RearrangementSupportingRead {
                    umi: "u1".into(),
                    read_name: "jread".into(),
                    sequence: b"TTGGAACCTTGGAACCTTGGAACC".to_vec(),
                    bam_is_reverse: false,
                    is_supplementary: false,
                    supports_v: false,
                    supports_j: true,
                    supports_d: false,
                    supports_c: false,
                    v_alignment: None,
                    j_alignment: None,
                    d_alignment: None,
                    c_alignment: None,
                },
            ],
            notation: "IGK:IGKV1-IGKJ1-?".into(),
        };
        assert!(reconstruct_candidate(&call, &reference).is_none());
    }

    #[test]
    fn translates_all_standard_codons_used_in_simple_example() {
        assert_eq!(translate_frame(b"ATGGCTTAA", 0), b"MA*");
        assert_eq!(translate_frame(b"NATGCT", 0), b"XA");
    }

    #[test]
    fn rejects_v_and_j_that_explain_same_query_region() {
        let sequence = b"TTTAACCGGTTAACCGGTTAAA";
        let germline = b"AACCGGTTAACCGGTT";
        assert!(!valid_candidate_geometry(sequence, germline, germline));
    }

}
