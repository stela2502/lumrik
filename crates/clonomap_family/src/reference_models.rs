//! ClonoMap discovery of reference-incompatible receptor fragments.
//! Persistent identity and evidence belong to `reference_curator`; this module
//! only decides which reconstructed V-side fragments are credible candidates.

use crate::{align_fragment, Receptor};
use reference_curator::{Observation, ReferenceCurator};
use sc_primer::{Chemistry, Grammar, PrimerDetector};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct ReferenceModelEntry {
    pub id: String,
    pub chain: String,
    pub original_v: String,
    pub observations: usize,
    pub evidence_reads: usize,
}

#[derive(Debug)]
pub struct ReferenceModels {
    curator: ReferenceCurator,
    entries: Vec<ReferenceModelEntry>,
    technical_detector: Option<PrimerDetector>,
}

impl Default for ReferenceModels {
    fn default() -> Self {
        Self { curator: ReferenceCurator::new(), entries: Vec::new(), technical_detector: None }
    }
}

impl ReferenceModels {
    pub fn open(path: Option<&Path>, chemistries: &[Chemistry]) -> Result<Self, String> {
        let grammars = chemistries.iter().copied().map(clonomap_technical_grammar).collect::<Result<Vec<_>, _>>()?;
        let detector = PrimerDetector::from_grammars(grammars).map_err(|e| e.to_string())?;
        let curator = match path {
            Some(path) => ReferenceCurator::open(path).map_err(|e| e.to_string())?,
            None => ReferenceCurator::new(),
        };
        Ok(Self { curator, entries: Vec::new(), technical_detector: Some(detector) })
    }

    pub fn resolve_receptor(&mut self, receptor: &mut Receptor, hard_identity: f64) -> bool {
        self.resolve_receptor_for_sample(receptor, hard_identity, None)
    }

    /// Register a repertoire-level mutation-distance outlier as a hypothetical
    /// reference model. Unlike `resolve_receptor_for_sample`, the decision that
    /// this receptor deserves a model has already been made by ClonoMap's global
    /// background pass, so there is no second identity threshold here. Sequence
    /// extraction/sanity checks remain exactly the same and persistent identity
    /// is still owned by reference_curator.
    pub fn hypothesize_receptor_for_sample(&mut self, receptor: &mut Receptor, sample: Option<&str>) -> bool {
        if !matches!(receptor.chain.as_str(), "IGH" | "IGK" | "IGL") { return false; }
        let original_v = receptor.v.clone();
        let Some(candidate_v) = candidate_v_fragment(receptor, self.technical_detector.as_ref()) else { return false; };

        if let Some(entry) = self.entries.iter_mut().find(|entry| {
            entry.chain == receptor.chain && entry.original_v == original_v &&
            self.curator.candidate(&entry.id).is_some_and(|candidate| {
                let model = candidate.resolved.as_ref().map(|r| r.sequence.as_slice()).unwrap_or(&candidate.sequence);
                align_fragment(&String::from_utf8_lossy(model), &candidate_v).is_some_and(|m| m.identity >= 0.75)
            })
        }) {
            entry.observations += 1;
            entry.evidence_reads += receptor.evidence_reads;
            receptor.v = entry.id.clone();
            return true;
        }

        let observation = Observation {
            source: "clonomap_family:global_mutation_outlier".into(),
            sample: sample.map(str::to_string),
            run: None,
            kind: format!("hypothetical_v:{}:{}", receptor.chain, original_v),
            support: receptor.evidence_reads as u64,
        };
        let Ok(id) = self.curator.observe(candidate_v.as_bytes(), observation) else { return false; };
        self.entries.push(ReferenceModelEntry {
            id: id.clone(), chain: receptor.chain.clone(), original_v, observations: 1, evidence_reads: receptor.evidence_reads,
        });
        receptor.v = id;
        true
    }

    pub fn resolve_receptor_for_sample(&mut self, receptor: &mut Receptor, hard_identity: f64, sample: Option<&str>) -> bool {
        if receptor.evidence_reads < 2 { return false; }
        if !matches!(receptor.chain.as_str(), "IGH" | "IGK" | "IGL") { return false; }
        let original_v = receptor.v.clone();
        let Some(candidate_v) = candidate_v_fragment(receptor, self.technical_detector.as_ref()) else { return false; };

        // Biological merging remains ClonoMap's responsibility. The curator is
        // deliberately exact: once ClonoMap has selected a model sequence it
        // returns the same stable candidate ID in every run/sample.
        if let Some(entry) = self.entries.iter_mut().find(|entry| {
            entry.chain == receptor.chain && entry.original_v == original_v &&
            self.curator.candidate(&entry.id).is_some_and(|candidate| {
                let model = candidate.resolved.as_ref().map(|r| r.sequence.as_slice()).unwrap_or(&candidate.sequence);
                align_fragment(&String::from_utf8_lossy(model), &candidate_v).is_some_and(|m| m.identity >= hard_identity)
            })
        }) {
            entry.observations += 1;
            entry.evidence_reads += receptor.evidence_reads;
            receptor.v = entry.id.clone();
            return true;
        }

        let Some(reference_fit) = align_fragment(&receptor.naive, &candidate_v) else { return false; };
        if reference_fit.identity >= hard_identity { return false; }

        let observation = Observation {
            source: "clonomap_family".into(),
            sample: sample.map(str::to_string),
            run: None,
            kind: format!("reconstructed_v:{}:{}", receptor.chain, original_v),
            support: receptor.evidence_reads as u64,
        };
        let Ok(id) = self.curator.observe(candidate_v.as_bytes(), observation) else { return false; };
        self.entries.push(ReferenceModelEntry {
            id: id.clone(), chain: receptor.chain.clone(), original_v, observations: 1, evidence_reads: receptor.evidence_reads,
        });
        receptor.v = id;
        true
    }

    pub fn effective_naive(&self, receptor: &Receptor) -> String {
        let Some(candidate) = self.curator.candidate(&receptor.v) else { return receptor.naive.clone(); };
        let model = candidate.resolved.as_ref().map(|r| r.sequence.as_slice()).unwrap_or(&candidate.sequence);
        splice_v_prefix(&receptor.naive, &String::from_utf8_lossy(model), &receptor.cdr3_nt).unwrap_or_else(|| receptor.naive.clone())
    }

    pub fn total(&self) -> usize { self.entries.len() }
    pub fn entries(&self) -> &[ReferenceModelEntry] { &self.entries }
    pub fn curator(&self) -> &ReferenceCurator { &self.curator }
    pub fn save(&self, path: &Path) -> Result<(), String> { self.curator.save(path).map_err(|e| e.to_string()) }
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

    #[test]
    fn unresolved_candidate_is_owned_by_reference_curator() {
        let mut models = ReferenceModels::open(None, &[Chemistry::BdV2_384]).unwrap();
        let mut receptor = Receptor {
            id: "r".into(), chain: "IGL".into(), v: "IGLV-old".into(), d: String::new(), j: "J".into(), c: String::new(),
            cdr3_nt: "TTTTGGGG".into(),
            naive: "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTGGGG".into(),
            observed: "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTTTGGGG".into(),
            alternative_reconstructions: 0, evidence_reads: 3,
        };
        assert!(models.resolve_receptor(&mut receptor, 0.75));
        let candidate = models.curator().candidate(&receptor.v).unwrap();
        assert!(candidate.resolved.is_none());
        assert_eq!(candidate.observations[0].source, "clonomap_family");
    }
}
