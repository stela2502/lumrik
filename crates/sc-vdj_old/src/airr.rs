use crate::posterior::RearrangementCall;
use crate::sequence::{translate_codon, translate_frame};

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AirrAnnotation {
    pub sequence: Vec<u8>,
    pub sequence_aa: Vec<u8>,
    pub productive: Option<bool>,
    pub vj_in_frame: Option<bool>,
    pub stop_codon: Option<bool>,
    pub junction: Vec<u8>,
    pub junction_aa: Vec<u8>,
    pub cdr3: Vec<u8>,
    pub cdr3_aa: Vec<u8>,
    pub cdr3_start: Option<usize>,
    pub cdr3_end: Option<usize>,
    pub v_cys_anchor: Option<bool>,
    pub j_anchor: Option<bool>,
}

pub fn annotate_rearrangement(call: &RearrangementCall) -> AirrAnnotation {
    let Some(j) = call.junction.as_ref() else {
        return AirrAnnotation::default();
    };

    let sequence = observed_rearrangement(j);
    if sequence.is_empty() || j.observed_v.is_empty() || j.observed_j.is_empty() {
        return AirrAnnotation { sequence, ..AirrAnnotation::default() };
    }

    let Some(v_anchor_ref) = terminal_v_cys(&j.naive_v) else {
        return AirrAnnotation { sequence, ..AirrAnnotation::default() };
    };
    let Some(j_anchor_ref) = initial_j_anchor(&j.naive_j) else {
        return AirrAnnotation { sequence, ..AirrAnnotation::default() };
    };

    // `naive_v` ends at the same V-side recombination boundary as `observed_v`,
    // so anchor placement is most stable when measured backwards from that end.
    let v_suffix = j.naive_v.len().saturating_sub(v_anchor_ref);
    if v_suffix > j.observed_v.len() {
        return AirrAnnotation { sequence, ..AirrAnnotation::default() };
    }
    let v_anchor = j.observed_v.len() - v_suffix;

    let j_offset = sequence.len().saturating_sub(j.observed_j.len());
    if j_anchor_ref + 3 > j.observed_j.len() {
        return AirrAnnotation { sequence, ..AirrAnnotation::default() };
    }
    let j_anchor = j_offset + j_anchor_ref;
    if v_anchor + 3 > sequence.len() || j_anchor + 3 > sequence.len() || v_anchor >= j_anchor {
        return AirrAnnotation { sequence, ..AirrAnnotation::default() };
    }

    let junction = sequence[v_anchor..j_anchor + 3].to_vec();
    let junction_aa = translate_frame(&junction, 0);
    let v_anchor_ok = junction_aa.first().copied() == Some(b'C');
    let j_anchor_ok = matches!(junction_aa.last().copied(), Some(b'W') | Some(b'F'));
    let vj_in_frame = junction.len() % 3 == 0;

    // The V cysteine establishes the coding phase for the entire reconstructed
    // sequence. This lets us report a protein even when the sequence begins
    // within the V gene rather than at its first coding base.
    let frame = v_anchor % 3;
    let sequence_aa = translate_frame(&sequence, frame);
    let stop_codon = sequence_aa.contains(&b'*');
    let productive = vj_in_frame && !stop_codon;

    // AIRR coordinates are 1-based closed intervals and CDR3 excludes the two
    // conserved codons included by `junction`.
    let cdr3 = if junction.len() >= 6 {
        junction[3..junction.len() - 3].to_vec()
    } else {
        Vec::new()
    };
    let cdr3_aa = if junction_aa.len() >= 2 {
        junction_aa[1..junction_aa.len() - 1].to_vec()
    } else {
        Vec::new()
    };
    let cdr3_start = (v_anchor + 3 < j_anchor).then_some(v_anchor + 4);
    let cdr3_end = (v_anchor + 3 < j_anchor).then_some(j_anchor);

    AirrAnnotation {
        sequence,
        sequence_aa,
        productive: Some(productive),
        vj_in_frame: Some(vj_in_frame),
        stop_codon: Some(stop_codon),
        junction,
        junction_aa,
        cdr3,
        cdr3_aa,
        cdr3_start,
        cdr3_end,
        v_cys_anchor: Some(v_anchor_ok),
        j_anchor: Some(j_anchor_ok),
    }
}

fn observed_rearrangement(j: &crate::junction::JunctionMeasurement) -> Vec<u8> {
    let mut out = Vec::with_capacity(
        j.observed_v.len()
            + j.p_v3.len()
            + j.n1.len()
            + j.p_d5.len()
            + j.observed_d.len()
            + j.p_d3.len()
            + j.n2.len()
            + j.p_j5.len()
            + j.observed_j.len(),
    );
    out.extend_from_slice(&j.observed_v);
    out.extend_from_slice(&j.p_v3);
    out.extend_from_slice(&j.n1);
    out.extend_from_slice(&j.p_d5);
    out.extend_from_slice(&j.observed_d);
    out.extend_from_slice(&j.p_d3);
    out.extend_from_slice(&j.n2);
    out.extend_from_slice(&j.p_j5);
    out.extend_from_slice(&j.observed_j);
    out
}

fn terminal_v_cys(v: &[u8]) -> Option<usize> {
    let window_start = v.len().saturating_sub(90);
    (window_start..v.len().saturating_sub(2))
        .rev()
        .find(|&pos| translate_codon(&v[pos..pos + 3]) == b'C')
}

fn initial_j_anchor(j: &[u8]) -> Option<usize> {
    let limit = j.len().min(90);
    (0..limit.saturating_sub(8)).find(|&pos| {
        let first = translate_codon(&j[pos..pos + 3]);
        let third = translate_codon(&j[pos + 6..pos + 9]);
        matches!(first, b'W' | b'F') && third == b'G'
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::junction::JunctionMeasurement;
    use crate::posterior::{GermlineSegmentSupport, RearrangementCall, RecombinationStage};
    use crate::types::{Chain, SegmentKind};

    fn support(id: &str, kind: SegmentKind) -> GermlineSegmentSupport {
        GermlineSegmentSupport {
            segment_index: 0,
            id: id.to_string(),
            kind,
            local_alignment_score: 100,
            supporting_umis: 2,
            supporting_reads: 2,
            locus_fraction: 0.0,
            distance_to_recombination_center: 0,
        }
    }

    fn call_with_junction(junction: JunctionMeasurement) -> RearrangementCall {
        RearrangementCall {
            chain: Chain::Igh,
            stage: RecombinationStage::Vdj,
            v: Some(support("IGHV1", SegmentKind::V)),
            d: Some(support("IGHD1", SegmentKind::D)),
            d_inferred_from_vj_junction: true,
            d_hypothesis_margin: Some(10),
            j: Some(support("IGHJ1", SegmentKind::J)),
            c: None,
            total_supporting_umis: 2,
            supporting_reads: Vec::new(),
            junction: Some(junction),
            notation: "IGH:IGHV1-IGHD1-IGHJ1".into(),
        }
    }

    #[test]
    fn airr_junction_includes_conserved_cys_and_j_anchor() {
        // V ends ...C, J begins WG; reconstructed junction translates C-A-W.
        let call = call_with_junction(JunctionMeasurement {
            v_del_3: 0,
            p_v3: Vec::new(),
            n1: b"GCT".to_vec(),
            p_d5: Vec::new(),
            d_del_5: Some(0),
            d_retained_len: Some(0),
            d_del_3: Some(0),
            p_d3: Vec::new(),
            n2: Vec::new(),
            p_j5: Vec::new(),
            j_del_5: 0,
            pn_alternative: false,
            observed_v: b"AAATGT".to_vec(),
            observed_d: Vec::new(),
            observed_j: b"TGGGCTGGTAAA".to_vec(),
            naive_v: b"AAATGT".to_vec(),
            naive_d: Vec::new(),
            naive_j: b"TGGGCTGGTAAA".to_vec(),
            observed_sequence: Vec::new(),
            inferred_naive_sequence: Vec::new(),
        });
        let a = annotate_rearrangement(&call);
        assert_eq!(a.junction, b"TGTGCTTGG");
        assert_eq!(a.junction_aa, b"CAW");
        assert_eq!(a.cdr3, b"GCT");
        assert_eq!(a.cdr3_aa, b"A");
        assert_eq!(a.productive, Some(true));
        assert_eq!(a.v_cys_anchor, Some(true));
        assert_eq!(a.j_anchor, Some(true));
        assert_eq!(a.cdr3_start, Some(7));
        assert_eq!(a.cdr3_end, Some(9));
    }

    #[test]
    fn stop_codon_makes_rearrangement_nonproductive() {
        let call = call_with_junction(JunctionMeasurement {
            v_del_3: 0,
            p_v3: Vec::new(),
            n1: b"TAA".to_vec(),
            p_d5: Vec::new(),
            d_del_5: Some(0),
            d_retained_len: Some(0),
            d_del_3: Some(0),
            p_d3: Vec::new(),
            n2: Vec::new(),
            p_j5: Vec::new(),
            j_del_5: 0,
            pn_alternative: false,
            observed_v: b"AAATGT".to_vec(),
            observed_d: Vec::new(),
            observed_j: b"TGGGCTGGTAAA".to_vec(),
            naive_v: b"AAATGT".to_vec(),
            naive_d: Vec::new(),
            naive_j: b"TGGGCTGGTAAA".to_vec(),
            observed_sequence: Vec::new(),
            inferred_naive_sequence: Vec::new(),
        });
        let a = annotate_rearrangement(&call);
        assert_eq!(a.junction_aa, b"C*W");
        assert_eq!(a.stop_codon, Some(true));
        assert_eq!(a.productive, Some(false));
    }
}
