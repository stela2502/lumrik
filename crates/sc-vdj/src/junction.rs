use crate::align::OrientedLocalAlignment;
use crate::sequence::reverse_complement;

const MAX_PALINDROME_LEN: usize = 12;

/// Inferred coding-end geometry for one observed V(D)J or VJ molecule.
///
/// Coordinates are measured against the selected germline segments. N bases are
/// retained for audit/reconstruction, but callers that need a mutation-stable
/// identity should use only their lengths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JunctionMeasurement {
    pub v_del_3: u16,
    pub p_v3: Vec<u8>,
    pub n1: Vec<u8>,
    pub p_d5: Vec<u8>,
    pub d_del_5: Option<u16>,
    pub d_retained_len: Option<u16>,
    pub d_del_3: Option<u16>,
    pub p_d3: Vec<u8>,
    pub n2: Vec<u8>,
    pub p_j5: Vec<u8>,
    pub j_del_5: u16,
    /// True when at least one alternative P/N decomposition also fits the same
    /// observed junction. The stored decomposition is the deterministic maximum-P one.
    pub pn_alternative: bool,
    /// Observed aligned V/D/J pieces in transcript orientation. These retain
    /// somatic mutations present in the molecule. `observed_d` is empty for VJ chains.
    pub observed_v: Vec<u8>,
    pub observed_d: Vec<u8>,
    pub observed_j: Vec<u8>,
    /// Retained germline V/D/J pieces after the inferred coding-end deletions.
    /// These are the mutation-free segment contributions used to reconstruct the
    /// naive recombination. `naive_d` is empty for VJ chains.
    pub naive_v: Vec<u8>,
    pub naive_d: Vec<u8>,
    pub naive_j: Vec<u8>,
    /// Observed molecule in transcript orientation.
    pub observed_sequence: Vec<u8>,
    /// Inferred mutation-free recombination: clipped germline V/D/J pieces plus
    /// the detected P and N nucleotides.
    pub inferred_naive_sequence: Vec<u8>,
}

impl JunctionMeasurement {
    pub fn p_v3_len(&self) -> u16 {
        self.p_v3.len() as u16
    }
    pub fn p_d5_len(&self) -> u16 {
        self.p_d5.len() as u16
    }
    pub fn p_d3_len(&self) -> u16 {
        self.p_d3.len() as u16
    }
    pub fn p_j5_len(&self) -> u16 {
        self.p_j5.len() as u16
    }
    pub fn n1_len(&self) -> u16 {
        self.n1.len() as u16
    }
    pub fn n2_len(&self) -> u16 {
        self.n2.len() as u16
    }
}

#[derive(Debug, Clone)]
pub struct JunctionInput<'a> {
    pub observed: &'a [u8],
    pub v: &'a [u8],
    pub v_alignment: OrientedLocalAlignment,
    pub d: Option<&'a [u8]>,
    pub d_alignment: Option<OrientedLocalAlignment>,
    pub j: &'a [u8],
    pub j_alignment: OrientedLocalAlignment,
}

/// Infer coding-end deletions, P additions and residual N sequence from one
/// continuous molecule. Returns None when the selected segments are not observed
/// in one orientation and in biological V-(D)-J order.
pub fn measure_junction(input: JunctionInput<'_>) -> Option<JunctionMeasurement> {
    let orientation = input.v_alignment.reverse_complement;
    if input.j_alignment.reverse_complement != orientation
        || input
            .d_alignment
            .is_some_and(|x| x.reverse_complement != orientation)
    {
        return None;
    }

    let observed = if orientation {
        reverse_complement(input.observed)
    } else {
        input.observed.to_vec()
    };
    let v = input.v_alignment.alignment;
    let j = input.j_alignment.alignment;
    let d = input.d_alignment.map(|x| x.alignment);

    if v.reference_end > input.v.len()
        || j.reference_start > input.j.len()
        || j.reference_end > input.j.len()
        || v.query_end > observed.len()
        || j.query_start > observed.len()
    {
        return None;
    }

    let v_del_3 = checked_u16(input.v.len().checked_sub(v.reference_end)?)?;
    let j_del_5 = checked_u16(j.reference_start)?;

    let (p_v3, n1, p_d5, d_del_5, d_retained_len, d_del_3, p_d3, n2, p_j5, ambiguous) =
        if let (Some(d_ref), Some(d)) = (input.d, d) {
            if d.reference_start > d.reference_end
                || d.reference_end > d_ref.len()
                || v.query_end > d.query_start
                || d.query_end > j.query_start
                || d.query_start > observed.len()
                || d.query_end > observed.len()
            {
                return None;
            }
            let left = &observed[v.query_end..d.query_start];
            let right = &observed[d.query_end..j.query_start];
            let retained_v = &input.v[..v.reference_end];
            let retained_d = &d_ref[d.reference_start..d.reference_end];
            let retained_j = &input.j[j.reference_start..];

            let left_split = split_junction(left, retained_v, retained_d);
            let right_split = split_junction(right, retained_d, retained_j);
            (
                left_split.left_p,
                left_split.n,
                left_split.right_p,
                Some(checked_u16(d.reference_start)?),
                Some(checked_u16(d.reference_end - d.reference_start)?),
                Some(checked_u16(d_ref.len() - d.reference_end)?),
                right_split.left_p,
                right_split.n,
                right_split.right_p,
                left_split.ambiguous || right_split.ambiguous,
            )
        } else {
            if v.query_end > j.query_start {
                return None;
            }
            let junction = &observed[v.query_end..j.query_start];
            let retained_v = &input.v[..v.reference_end];
            let retained_j = &input.j[j.reference_start..];
            let split = split_junction(junction, retained_v, retained_j);
            (
                split.left_p,
                split.n,
                Vec::new(),
                None,
                None,
                None,
                Vec::new(),
                Vec::new(),
                split.right_p,
                split.ambiguous,
            )
        };

    let observed_v = observed[v.query_start..v.query_end].to_vec();
    let observed_d = d
        .map(|d| observed[d.query_start..d.query_end].to_vec())
        .unwrap_or_default();
    let observed_j = observed[j.query_start..j.query_end].to_vec();
    let naive_v = input.v[..v.reference_end].to_vec();
    let naive_d = if let (Some(d_ref), Some(d)) = (input.d, d) {
        d_ref[d.reference_start..d.reference_end].to_vec()
    } else {
        Vec::new()
    };
    let naive_j = input.j[j.reference_start..].to_vec();

    let mut naive = Vec::new();
    naive.extend_from_slice(&naive_v);
    naive.extend_from_slice(&p_v3);
    naive.extend_from_slice(&n1);
    naive.extend_from_slice(&p_d5);
    if input.d.is_some() && d.is_some() {
        naive.extend_from_slice(&naive_d);
        naive.extend_from_slice(&p_d3);
        naive.extend_from_slice(&n2);
        naive.extend_from_slice(&p_j5);
    }
    naive.extend_from_slice(&naive_j);

    Some(JunctionMeasurement {
        v_del_3,
        p_v3,
        n1,
        p_d5,
        d_del_5,
        d_retained_len,
        d_del_3,
        p_d3,
        n2,
        p_j5,
        j_del_5,
        pn_alternative: ambiguous,
        observed_v,
        observed_d,
        observed_j,
        naive_v,
        naive_d,
        naive_j,
        observed_sequence: observed,
        inferred_naive_sequence: naive,
    })
}

#[derive(Debug)]
struct JunctionSplit {
    left_p: Vec<u8>,
    n: Vec<u8>,
    right_p: Vec<u8>,
    ambiguous: bool,
}

/// Split one observed inter-segment junction into coding-end P from the left,
/// residual N, and coding-end P from the right. Every exact P-compatible split
/// is considered. The canonical result maximizes total P length, then left P;
/// ambiguity is retained when more than one split fits.
fn split_junction(junction: &[u8], left_retained: &[u8], right_retained: &[u8]) -> JunctionSplit {
    let left_max = MAX_PALINDROME_LEN
        .min(left_retained.len())
        .min(junction.len());
    let right_max = MAX_PALINDROME_LEN
        .min(right_retained.len())
        .min(junction.len());

    let mut candidates = Vec::new();
    for lp in 0..=left_max {
        let expected_left = coding_end_palindrome_3(left_retained, lp);
        if junction.get(..lp) != Some(expected_left.as_slice()) {
            continue;
        }
        for rp in 0..=right_max.min(junction.len() - lp) {
            let expected_right = coding_end_palindrome_5(right_retained, rp);
            if junction.get(junction.len() - rp..) != Some(expected_right.as_slice()) {
                continue;
            }
            candidates.push((lp, rp));
        }
    }

    candidates.sort_by(|a, b| {
        (b.0 + b.1)
            .cmp(&(a.0 + a.1))
            .then_with(|| b.0.cmp(&a.0))
            .then_with(|| b.1.cmp(&a.1))
    });
    let (lp, rp) = candidates.first().copied().unwrap_or((0, 0));
    JunctionSplit {
        left_p: junction[..lp].to_vec(),
        n: junction[lp..junction.len() - rp].to_vec(),
        right_p: junction[junction.len() - rp..].to_vec(),
        ambiguous: candidates.len() > 1,
    }
}

fn coding_end_palindrome_3(retained: &[u8], len: usize) -> Vec<u8> {
    reverse_complement(&retained[retained.len() - len..])
}

fn coding_end_palindrome_5(retained: &[u8], len: usize) -> Vec<u8> {
    reverse_complement(&retained[..len])
}

fn checked_u16(value: usize) -> Option<u16> {
    u16::try_from(value).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{LocalAlignment, OrientedLocalAlignment};

    fn oa(qs: usize, qe: usize, rs: usize, re: usize) -> OrientedLocalAlignment {
        OrientedLocalAlignment {
            alignment: LocalAlignment {
                score: 100,
                query_start: qs,
                query_end: qe,
                reference_start: rs,
                reference_end: re,
            },
            reverse_complement: false,
        }
    }

    #[test]
    fn splits_four_coding_end_palindromes_and_n() {
        let v = b"AAAACG";
        let d = b"TAGCCT";
        let j = b"GATTTT";
        // retain V ...CG => P=CG; retain D starts AG => P=CT;
        // retain D ends CC => P=GG; retain J starts AT => P=AT.
        let observed = b"AAAACGCGAACTAGCCGGTCATATTTT";
        let m = measure_junction(JunctionInput {
            observed,
            v,
            v_alignment: oa(0, 6, 0, 6),
            d: Some(d),
            d_alignment: Some(oa(12, 16, 1, 5)),
            j,
            j_alignment: oa(22, 27, 1, 6),
        })
        .unwrap();
        assert_eq!(m.v_del_3, 0);
        assert_eq!(m.p_v3, b"CG");
        assert_eq!(m.n1, b"AA");
        assert_eq!(m.p_d5, b"CT");
        assert_eq!(m.d_del_5, Some(1));
        assert_eq!(m.d_retained_len, Some(4));
        assert_eq!(m.d_del_3, Some(1));
        assert_eq!(m.p_d3, b"GG");
        assert_eq!(m.n2, b"TC");
        assert_eq!(m.p_j5, b"AT");
        assert!(m.pn_alternative);
    }

    #[test]
    fn vj_measurement_reconstructs_germline_derived_sequence() {
        let v = b"AACCGG";
        let j = b"TTGGAA";
        let observed = b"AACCGGCCATCCGGAA";
        let m = measure_junction(JunctionInput {
            observed,
            v,
            v_alignment: oa(0, 6, 0, 6),
            d: None,
            d_alignment: None,
            j,
            j_alignment: oa(12, 16, 2, 6),
        })
        .unwrap();
        assert_eq!(m.j_del_5, 2);
        assert_eq!(m.p_v3, b"CC");
        assert_eq!(m.n1, b"A");
        assert_eq!(m.p_j5, b"TCC");
        assert!(m.pn_alternative);
        assert_eq!(m.observed_v, b"AACCGG");
        assert!(m.observed_d.is_empty());
        assert_eq!(m.observed_j, b"GGAA");
        assert_eq!(m.naive_v, b"AACCGG");
        assert!(m.naive_d.is_empty());
        assert_eq!(m.naive_j, b"GGAA");
        assert_eq!(m.inferred_naive_sequence, b"AACCGGCCATCCGGAA");
    }
}
