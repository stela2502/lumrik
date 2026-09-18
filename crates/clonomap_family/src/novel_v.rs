//! Provisional V-segment registration for receptors that cannot be explained by
//! the supplied germline reconstruction. This is analysis-local evidence: no
//! reference index is modified and no database allele name is claimed.

use crate::{align_fragment, Receptor};
use sc_primer::{Chemistry, Grammar, PrimerDetector};

#[derive(Debug, Clone)]
pub struct NovelVEntry {
    pub name: String,
    pub chain: String,
    pub original_v: String,
    /// Best currently available observed fragment for this provisional V.
    pub fragment: String,
    pub observations: usize,
    pub evidence_reads: usize,
}

#[derive(Debug, Default, Clone)]
pub struct NovelVFinder {
    entries: Vec<NovelVEntry>,
}

impl NovelVFinder {
    fn prefix(chain: &str) -> &'static str {
        match chain { "IGH" => "IGHV", "IGK" => "IGKV", "IGL" => "IGLV", _ => "V" }
    }

    fn resolve(&mut self, receptor: &mut Receptor, hard_identity: f64, technical_detector: Option<&PrimerDetector>) -> bool {
        // A provisional V is an escape hatch for a well-supported receptor whose
        // own supplied naive reconstruction cannot explain its observed sequence.
        // Single-read reconstructions never create registry entries.
        if receptor.evidence_reads < 2 { return false; }

        let original_v = receptor.v.clone();
        let Some(candidate_v) = candidate_v_fragment(receptor, technical_detector) else { return false; };

        // Once a provisional V exists, test it before asking the supplied
        // reference to veto discovery again.  A later receptor that matches an
        // already registered replacement is evidence for that replacement; it
        // must not be lost merely because the coverage-aware reference walker
        // can find a locally compatible patch in the old naive reconstruction.
        if let Some(entry) = self.entries.iter_mut().find(|entry| {
            entry.chain == receptor.chain
                && entry.original_v == original_v
                && align_fragment(&entry.fragment, &candidate_v)
                    .is_some_and(|m| m.identity >= hard_identity)
        }) {
            entry.observations += 1;
            entry.evidence_reads += receptor.evidence_reads;
            // Until we add a base-voting merger, never throw away span: retain
            // the longest supported fragment as the current provisional model.
            if candidate_v.len() > entry.fragment.len() {
                entry.fragment = candidate_v.clone();
            }
            receptor.v = entry.name.clone();
            return true;
        }

        // No registered replacement explains this receptor.  Only now decide
        // whether its supplied reference is bad enough to seed a new
        // provisional V identity.
        let Some(reference_fit) = align_fragment(&receptor.naive, &candidate_v) else { return false; };
        if reference_fit.identity >= hard_identity { return false; }

        let name = format!("{}-novel{}", Self::prefix(&receptor.chain), self.entries.len() + 1);
        self.entries.push(NovelVEntry {
            name: name.clone(),
            chain: receptor.chain.clone(),
            original_v,
            fragment: candidate_v,
            observations: 1,
            evidence_reads: receptor.evidence_reads,
        });
        receptor.v = name;
        true
    }

    fn get(&self, name: &str) -> Option<&NovelVEntry> {
        self.entries.iter().find(|entry| entry.name == name)
    }

    pub fn entries(&self) -> &[NovelVEntry] { &self.entries }
}

#[derive(Debug, Clone)]
pub struct NovelVRegistries {
    pub hc: NovelVFinder,
    pub lc: NovelVFinder,
    technical_detector: Option<PrimerDetector>,
}

impl Default for NovelVRegistries {
    fn default() -> Self {
        Self { hc: NovelVFinder::default(), lc: NovelVFinder::default(), technical_detector: None }
    }
}

impl NovelVRegistries {
    pub fn with_chemistries(chemistries: &[Chemistry]) -> Result<Self, String> {
        // ClonoMap sees reconstructed receptor fragments, so the original read
        // boundary is gone.  Reuse sc_primer's real BD grammars but deliberately
        // widen their leading search window: known technical structure may begin
        // anywhere in the first 150 nt of the reconstructed fragment.
        let grammars = chemistries.iter().copied()
            .map(clonomap_technical_grammar)
            .collect::<Result<Vec<_>, _>>()?;
        let detector = PrimerDetector::from_grammars(grammars)
            .map_err(|e| e.to_string())?;
        Ok(Self { hc: NovelVFinder::default(), lc: NovelVFinder::default(), technical_detector: Some(detector) })
    }

    pub fn resolve_receptor(&mut self, receptor: &mut Receptor, hard_identity: f64) -> bool {
        let detector = self.technical_detector.as_ref();
        match receptor.chain.as_str() {
            "IGH" => self.hc.resolve(receptor, hard_identity, detector),
            "IGK" | "IGL" => self.lc.resolve(receptor, hard_identity, detector),
            _ => false,
        }
    }

    pub fn entry(&self, name: &str) -> Option<&NovelVEntry> {
        self.hc.get(name).or_else(|| self.lc.get(name))
    }

    /// Return the naive reconstruction to use *now*. Ordinary receptors keep
    /// the sc-vdj naive verbatim. A provisional V replaces only the prefix up
    /// to a rightmost exact anchor shared by the old naive and the registry's
    /// current fragment; the original junction/D/J suffix is retained.
    pub fn effective_naive(&self, receptor: &Receptor) -> String {
        let Some(entry) = self.entry(&receptor.v) else { return receptor.naive.clone(); };
        splice_v_prefix(&receptor.naive, &entry.fragment, &receptor.cdr3_nt).unwrap_or_else(|| receptor.naive.clone())
    }

    pub fn total(&self) -> usize { self.hc.entries.len() + self.lc.entries.len() }
}



fn clonomap_technical_grammar(chemistry: Chemistry) -> Result<Grammar, String> {
    match chemistry {
        Chemistry::BdV2_384 => Grammar::parse(
            "clonomap-bd-v2-384-wide",
            "TYPE:GEX+SEARCH:0..150+BD_CELL:v2.384+POLYT:min=0",
        ).map_err(|e| e.to_string()),
        Chemistry::BdV2_384Vdj => Grammar::parse(
            "clonomap-bd-v2-384-vdj-wide",
            "TYPE:VDJ+SEARCH:0..150+FIXED:ACAGGAAACTCATGGTGCGT:mm=2+BD_CELL:v2.384-vdj",
        ).map_err(|e| e.to_string()),
        other => other.grammar().map_err(|e| e.to_string()),
    }
}

const HOMOPOLYMER_CLIP: usize = 10;
const MIN_ORF_NT: usize = 30;

/// Extract only the observed V-side fragment that is allowed to seed/reuse a
/// provisional V identity.  The CDR3 start is our right-hand boundary.  Long
/// homopolymers are treated as upstream contamination and clipped away; if the
/// remaining V-adjacent sequence does not contain a stop-free coding suffix in
/// any frame, it is not credible enough to register as a V segment.
fn candidate_v_fragment(receptor: &Receptor, technical_detector: Option<&PrimerDetector>) -> Option<String> {
    if receptor.cdr3_nt.is_empty() { return None; }
    let boundary = receptor.observed.find(&receptor.cdr3_nt)?;
    let raw = &receptor.observed[..boundary];
    if raw.len() < MIN_ORF_NT { return None; }

    // sc_primer is the authority for known assay structure. If a complete
    // selected-chemistry primer is embedded on the V side, discard everything
    // through that technical structure and retain only the sequence toward CDR3.
    // Failure to recognize a primer is not itself a reason to reject biological
    // evidence: the remaining sanity checks still apply.
    let technical_start = technical_prefix_end(raw.as_bytes(), technical_detector).unwrap_or(0);
    let raw = &raw[technical_start..];
    if raw.len() < MIN_ORF_NT { return None; }

    // A long homopolymer is not a novel germline feature.  Because the V
    // contribution must end immediately before CDR3, salvage only the suffix
    // after the last such run.  If that leaves too little sequence, reject it.
    let start = last_homopolymer_end(raw.as_bytes(), HOMOPOLYMER_CLIP).unwrap_or(0);
    let clipped = &raw[start..];
    if clipped.len() < MIN_ORF_NT { return None; }

    // We do not assume the frame supplied by the old reference is trustworthy.
    // Ask only whether *some* frame has a stop-free coding suffix adjacent to
    // the recombination boundary.  Keep that suffix; this simultaneously clips
    // non-coding 5' junk without inventing genomic coordinates.
    let start = longest_stop_free_suffix_start(clipped.as_bytes())?;
    let candidate = &clipped[start..];
    (candidate.len() >= MIN_ORF_NT).then(|| candidate.to_string())
}

fn technical_prefix_end(seq: &[u8], detector: Option<&PrimerDetector>) -> Option<usize> {
    let detector = detector?;
    let qual = vec![b'I'; seq.len()];
    detector.detect_all(seq, &qual).ok()?.into_iter()
        // We are clipping an upstream technical structure, not arbitrary
        // sequence near the recombination boundary. Keep the rightmost end of
        // any complete primer that starts in the upstream half of the candidate.
        .filter(|hit| hit.primer_start < seq.len() / 2)
        .map(|hit| hit.primer_end)
        .filter(|&end| end <= seq.len())
        .max()
}

fn last_homopolymer_end(seq: &[u8], min_run: usize) -> Option<usize> {
    if seq.is_empty() { return None; }
    let mut last = None;
    let mut run_start = 0usize;
    for i in 1..=seq.len() {
        if i == seq.len() || seq[i] != seq[run_start] {
            if i - run_start >= min_run {
                last = Some(i);
            }
            run_start = i;
        }
    }
    last
}

fn longest_stop_free_suffix_start(seq: &[u8]) -> Option<usize> {
    let mut best: Option<(usize, usize)> = None; // (length, start)
    for frame in 0..3 {
        if frame >= seq.len() { continue; }
        let mut last_stop_end = frame;
        let mut i = frame;
        while i + 3 <= seq.len() {
            if matches!(&seq[i..i + 3], b"TAA" | b"TAG" | b"TGA") {
                last_stop_end = i + 3;
            }
            i += 3;
        }
        let len = seq.len().saturating_sub(last_stop_end);
        if best.is_none_or(|(best_len, _)| len > best_len) {
            best = Some((len, last_stop_end));
        }
    }
    best.map(|(_, start)| start)
}

fn splice_v_prefix(old_naive: &str, novel_fragment: &str, cdr3_nt: &str) -> Option<String> {
    if cdr3_nt.is_empty() { return None; }
    let naive_boundary = old_naive.find(cdr3_nt)?;
    let mut out = String::with_capacity(novel_fragment.len() + old_naive.len() - naive_boundary);
    out.push_str(novel_fragment);
    out.push_str(&old_naive[naive_boundary..]);
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn receptor(v: &str, naive: &str, observed: &str, reads: usize) -> Receptor {
        Receptor { id:"r".into(), chain:"IGL".into(), v:v.into(), d:String::new(), j:"J".into(), c:String::new(), cdr3_nt:"TTTTGGGG".into(), naive:naive.into(), observed:observed.into(), alternative_reconstructions:0, evidence_reads:reads }
    }

    #[test]
    fn bad_supported_reference_changes_identity_and_is_reused() {
        let mut regs = NovelVRegistries::default();
        let mut a = receptor("Iglv3", "AAAACCCCGGGGAAAACCGCCGCCGCCGCCGCCGCCGTTTTGGGGCCCC", "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGTTTTGGGGCCCC", 3);
        let mut b = receptor("Iglv3", "AAAACCCCGGGGAAAACCGCCGCCGCCGCCGCCGCCGTTTTGGGGCCCC", "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCATTTTGGGGCCCC", 4);
        assert!(regs.resolve_receptor(&mut a, 0.75));
        assert_eq!(a.v, "IGLV-novel1");
        assert!(regs.resolve_receptor(&mut b, 0.75));
        assert_eq!(b.v, "IGLV-novel1");
        assert_eq!(regs.lc.entries()[0].observations, 2);
        assert_eq!(regs.lc.entries()[0].evidence_reads, 7);
    }

    #[test]
    fn single_read_cannot_register_novel_v() {
        let mut regs = NovelVRegistries::default();
        let mut r = receptor("Iglv3", "AAAACCCCGGGGAAAACCGCCGCCGCCGCCGCCGCCGTTTTGGGGCCCC", "GCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGTTTTGGGGCCCC", 1);
        assert!(!regs.resolve_receptor(&mut r, 0.75));
        assert_eq!(r.v, "Iglv3");
    }

    #[test]
    fn candidate_v_clips_upstream_homopolymer_and_keeps_orf_suffix() {
        let r = receptor("Iglv3", "AAAACCCCGGGGAAAATTTTGGGGCCCC", "AAAAAAAAAAAACCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGTTTTGGGG", 3);
        let v = candidate_v_fragment(&r, None).unwrap();
        assert!(!v.starts_with("AAAAAAAAAA"));
        assert!(v.ends_with("CCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCG"));
    }

    #[test]
    fn candidate_v_rejects_homopolymer_that_leaves_no_v_evidence() {
        let r = receptor("Iglv3", "AAAACCCCGGGGAAAATTTTGGGGCCCC", "CCGCCGAAAAAAAAAAAAAAAAAAAAATTTTGGGG", 3);
        assert!(candidate_v_fragment(&r, None).is_none());
    }

    #[test]
    fn effective_naive_preserves_downstream_reference_suffix() {
        let old = "AAAACCCCGGGGAAAATTTTGGGGCCCC";
        let novel = "CCGCCGCCGCCGCCGCCGCCG";
        let spliced = splice_v_prefix(old, novel, "TTTTGGGG").unwrap();
        assert!(spliced.ends_with("TTTTGGGGCCCC"));
        assert!(spliced.starts_with(novel));
    }
}
