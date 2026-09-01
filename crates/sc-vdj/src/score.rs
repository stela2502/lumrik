use crate::gex::ExpressionMatrix;

#[derive(Debug, Clone)]
pub struct MarkerContribution {
    pub gene: &'static str,
    pub expression: f64,
    pub weight: f64,
    pub contribution: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DevelopmentProgram {
    None,
    BLineageRecombining,
    TLineageRecombining,
}

#[derive(Debug, Clone)]
pub struct DevelopmentEvidence {
    pub program: DevelopmentProgram,
    pub probability_like_score: f64,
    pub rag_activity: f64,
    pub contributions: Vec<MarkerContribution>,
    /// Exact, reversible compact rationale. High 16 bits are confidence (0..65535),
    /// then program and marker-presence bits. Sorting this integer is therefore
    /// primarily sorting by confidence, while decode() recovers the rationale flags.
    pub evidence_code: u64,
}

const B_MARKERS: [(&str, f64); 8] = [
    ("RAG1", 2.0),
    ("RAG2", 2.0),
    ("DNTT", 1.0),
    ("PAX5", 1.0),
    ("EBF1", 1.0),
    ("VPREB1", 1.5),
    ("IGLL1", 1.5),
    ("CD79A", 0.5),
];
const T_MARKERS: [(&str, f64); 8] = [
    ("RAG1", 2.0),
    ("RAG2", 2.0),
    ("DNTT", 1.0),
    ("PTCRA", 1.5),
    ("NOTCH1", 1.0),
    ("LCK", 0.5),
    ("CD3D", 0.5),
    ("CD3E", 0.5),
];

pub fn score_development<E: ExpressionMatrix>(
    gex: &E,
    cell: &str,
    b_locus: f64,
    t_locus: f64,
) -> DevelopmentEvidence {
    let (b, bc, bmask) = marker_score(gex, cell, &B_MARKERS);
    let (t, tc, tmask) = marker_score(gex, cell, &T_MARKERS);
    let rag1 = squash(gex.expression(cell, "RAG1"));
    let rag2 = squash(gex.expression(cell, "RAG2"));
    let rag = (rag1 * rag2).sqrt();
    let b_total = 0.55 * b + 0.30 * b_locus + 0.15 * rag;
    let t_total = 0.55 * t + 0.30 * t_locus + 0.15 * rag;
    let (program, score, contrib, mask) = if b_total < 0.08 && t_total < 0.08 {
        (
            DevelopmentProgram::None,
            b_total.max(t_total),
            Vec::new(),
            0u16,
        )
    } else if b_total >= t_total {
        (DevelopmentProgram::BLineageRecombining, b_total, bc, bmask)
    } else {
        (DevelopmentProgram::TLineageRecombining, t_total, tc, tmask)
    };
    let conf = (score.clamp(0.0, 1.0) * 65535.0).round() as u64;
    let pbits = match program {
        DevelopmentProgram::None => 0u64,
        DevelopmentProgram::BLineageRecombining => 1,
        DevelopmentProgram::TLineageRecombining => 2,
    };
    let evidence_code = (conf << 48) | (pbits << 40) | mask as u64;
    DevelopmentEvidence {
        program,
        probability_like_score: score.clamp(0.0, 1.0),
        rag_activity: rag,
        contributions: contrib,
        evidence_code,
    }
}

pub fn decode_evidence_code(code: u64) -> (f64, DevelopmentProgram, u16) {
    let confidence = ((code >> 48) & 0xffff) as f64 / 65535.0;
    let program = match (code >> 40) & 0xff {
        1 => DevelopmentProgram::BLineageRecombining,
        2 => DevelopmentProgram::TLineageRecombining,
        _ => DevelopmentProgram::None,
    };
    (confidence, program, (code & 0xffff) as u16)
}

fn marker_score<E: ExpressionMatrix>(
    gex: &E,
    cell: &str,
    markers: &[(&'static str, f64)],
) -> (f64, Vec<MarkerContribution>, u16) {
    let total_w: f64 = markers.iter().map(|x| x.1).sum();
    let mut sum = 0.0;
    let mut out = Vec::new();
    let mut mask = 0u16;
    for (i, (gene, w)) in markers.iter().enumerate() {
        let x = gex.expression(cell, gene);
        let s = squash(x);
        let c = s * *w;
        sum += c;
        if x > 0.0 {
            mask |= 1u16 << i;
        }
        out.push(MarkerContribution {
            gene,
            expression: x,
            weight: *w,
            contribution: c / total_w,
        });
    }
    (sum / total_w, out, mask)
}
fn squash(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else {
        1.0 - (-x / 2.0).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn code_roundtrips_headline() {
        let code = (50000u64 << 48) | (1u64 << 40) | 0x55;
        let (c, p, m) = decode_evidence_code(code);
        assert!(c > 0.7);
        assert_eq!(p, DevelopmentProgram::BLineageRecombining);
        assert_eq!(m, 0x55);
    }
}
