use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const STORE_MAGIC: &[u8; 8] = b"LUMREFC1";
const STORE_VERSION: u32 = 1;

pub type CandidateId = String;
pub type ReferenceId = String;
pub type EvidenceId = u64;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReferenceCurator {
    version: u32,
    candidates: HashMap<CandidateId, ReferenceCandidate>,
    sequence_index: HashMap<Vec<u8>, CandidateId>,
    next_candidate_number: u64,
    next_evidence_id: EvidenceId,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReferenceCandidate {
    pub id: CandidateId,
    /// The sequence as reconstructed/observed. It is never replaced by resolution.
    pub sequence: Vec<u8>,
    /// Curated reference representation, if this candidate has been resolved.
    pub resolved: Option<ResolvedReference>,
    pub observations: Vec<Observation>,
    /// Alignment evidence grouped by the exact reference/assembly it was produced against.
    pub annotations: HashMap<ReferenceId, Vec<Annotation>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Observation {
    pub source: String,
    pub sample: Option<String>,
    pub run: Option<String>,
    pub kind: String,
    pub support: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ResolvedReference {
    pub sequence: Vec<u8>,
    pub source: String,
    pub evidence_ids: Vec<EvidenceId>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Annotation {
    pub evidence_id: EvidenceId,
    /// Program that produced the alignment: STAR, BWA, minimap2, BLAST, ...
    pub producer: Producer,
    pub target: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub cigar: Option<String>,
    pub mapq: Option<u8>,
    pub edit_distance: Option<u32>,
    pub query_start: u32,
    pub query_end: u32,
    pub query_len: u32,
    pub secondary: bool,
    pub supplementary: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Producer {
    pub name: String,
    pub version: Option<String>,
    pub parameters: Option<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateFilter {
    All,
    ResolvedOnly,
    UnresolvedOnly,
}

impl Default for ReferenceCurator {
    fn default() -> Self {
        Self::new()
    }
}

impl ReferenceCurator {
    pub fn new() -> Self {
        Self {
            version: STORE_VERSION,
            candidates: HashMap::new(),
            sequence_index: HashMap::new(),
            next_candidate_number: 1,
            next_evidence_id: 1,
        }
    }

    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        if !path.exists() {
            return Ok(Self::new());
        }
        let mut reader = BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
        let magic: [u8; 8] = bincode::deserialize_from(&mut reader).context("read reference-curator magic")?;
        if &magic != STORE_MAGIC {
            bail!("{} is not a reference_curator store", path.display());
        }
        let store: Self = bincode::deserialize_from(&mut reader).context("read reference-curator store")?;
        if store.version != STORE_VERSION {
            bail!("unsupported reference_curator store version {} (expected {})", store.version, STORE_VERSION);
        }
        store.validate()?;
        Ok(store)
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        self.validate()?;
        let path = path.as_ref();
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).with_context(|| format!("create {}", parent.display()))?;
        }
        let tmp = temporary_path(path);
        {
            let mut writer = BufWriter::new(File::create(&tmp).with_context(|| format!("create {}", tmp.display()))?);
            bincode::serialize_into(&mut writer, STORE_MAGIC).context("write reference-curator magic")?;
            bincode::serialize_into(&mut writer, self).context("write reference-curator store")?;
            writer.flush()?;
        }
        fs::rename(&tmp, path).with_context(|| format!("replace {}", path.display()))?;
        Ok(())
    }

    /// Register an observed sequence. Exact sequence identity is deliberately strict:
    /// fuzzy biological merging belongs in the producer, not in the curator.
    pub fn observe(&mut self, sequence: impl AsRef<[u8]>, observation: Observation) -> Result<CandidateId> {
        let sequence = normalize_sequence(sequence.as_ref())?;
        if let Some(id) = self.sequence_index.get(&sequence).cloned() {
            self.candidates.get_mut(&id).expect("sequence index must point to candidate").observations.push(observation);
            return Ok(id);
        }
        let id = format!("RC-{:08}", self.next_candidate_number);
        self.next_candidate_number += 1;
        self.sequence_index.insert(sequence.clone(), id.clone());
        self.candidates.insert(id.clone(), ReferenceCandidate {
            id: id.clone(), sequence, resolved: None, observations: vec![observation], annotations: HashMap::new(),
        });
        Ok(id)
    }

    pub fn add_annotation(&mut self, candidate_id: &str, reference: impl Into<String>, mut annotation: Annotation) -> Result<EvidenceId> {
        let candidate = self.candidates.get_mut(candidate_id).with_context(|| format!("unknown candidate {candidate_id}"))?;
        let reference = reference.into();
        if reference.trim().is_empty() { bail!("reference id must not be empty"); }
        let evidence_id = self.next_evidence_id;
        self.next_evidence_id += 1;
        annotation.evidence_id = evidence_id;
        candidate.annotations.entry(reference).or_default().push(annotation);
        Ok(evidence_id)
    }

    pub fn resolve(&mut self, candidate_id: &str, resolved: ResolvedReference) -> Result<()> {
        let candidate = self.candidates.get_mut(candidate_id).with_context(|| format!("unknown candidate {candidate_id}"))?;
        for evidence_id in &resolved.evidence_ids {
            if !candidate.annotations.values().flatten().any(|a| a.evidence_id == *evidence_id) {
                bail!("candidate {candidate_id} has no evidence id {evidence_id}");
            }
        }
        candidate.resolved = Some(resolved);
        Ok(())
    }

    pub fn candidate(&self, id: &str) -> Option<&ReferenceCandidate> { self.candidates.get(id) }
    pub fn candidates(&self) -> impl Iterator<Item = &ReferenceCandidate> { self.candidates.values() }
    pub fn len(&self) -> usize { self.candidates.len() }
    pub fn is_empty(&self) -> bool { self.candidates.is_empty() }
    pub fn unresolved(&self) -> impl Iterator<Item = &ReferenceCandidate> { self.candidates.values().filter(|c| c.resolved.is_none()) }
    pub fn resolved(&self) -> impl Iterator<Item = &ReferenceCandidate> { self.candidates.values().filter(|c| c.resolved.is_some()) }

    pub fn export_fasta(&self, path: impl AsRef<Path>, filter: CandidateFilter) -> Result<()> {
        let mut candidates: Vec<_> = self.candidates.values().collect();
        candidates.sort_by(|a, b| a.id.cmp(&b.id));
        let mut out = BufWriter::new(File::create(path.as_ref())?);
        for candidate in candidates {
            let sequence = match (&candidate.resolved, filter) {
                (Some(resolved), CandidateFilter::All | CandidateFilter::ResolvedOnly) => &resolved.sequence,
                (None, CandidateFilter::All | CandidateFilter::UnresolvedOnly) => &candidate.sequence,
                (Some(_), CandidateFilter::UnresolvedOnly) | (None, CandidateFilter::ResolvedOnly) => continue,
            };
            writeln!(out, ">{}", candidate.id)?;
            writeln!(out, "{}", String::from_utf8_lossy(sequence))?;
        }
        Ok(())
    }

    fn validate(&self) -> Result<()> {
        if self.version != STORE_VERSION { bail!("invalid store version {}", self.version); }
        if self.sequence_index.len() != self.candidates.len() { bail!("sequence index/candidate count mismatch"); }
        for (sequence, id) in &self.sequence_index {
            let candidate = self.candidates.get(id).with_context(|| format!("sequence index points to missing candidate {id}"))?;
            if &candidate.sequence != sequence { bail!("sequence index mismatch for candidate {id}"); }
        }
        Ok(())
    }
}

impl Annotation {
    pub fn new(producer: Producer, target: impl Into<String>, start: u64, end: u64, strand: Strand, query_len: u32) -> Self {
        Self { evidence_id: 0, producer, target: target.into(), start, end, strand, cigar: None, mapq: None, edit_distance: None,
            query_start: 0, query_end: query_len, query_len, secondary: false, supplementary: false }
    }
}

fn normalize_sequence(sequence: &[u8]) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(sequence.len());
    for &base in sequence {
        if base.is_ascii_whitespace() { continue; }
        let base = base.to_ascii_uppercase();
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') { bail!("unsupported sequence base {:?}", base as char); }
        out.push(base);
    }
    if out.is_empty() { bail!("candidate sequence must not be empty"); }
    Ok(out)
}

fn temporary_path(path: &Path) -> PathBuf {
    let mut tmp = path.as_os_str().to_os_string();
    tmp.push(".tmp");
    PathBuf::from(tmp)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn observation(sample: &str) -> Observation {
        Observation { source: "clonomap_family".into(), sample: Some(sample.into()), run: None, kind: "reconstructed_v".into(), support: 7 }
    }

    #[test]
    fn exact_sequence_reuses_stable_candidate() {
        let mut curator = ReferenceCurator::new();
        let first = curator.observe(b"acgt", observation("a")).unwrap();
        let second = curator.observe(b"ACGT\n", observation("b")).unwrap();
        assert_eq!(first, second);
        assert_eq!(curator.len(), 1);
        assert_eq!(curator.candidate(&first).unwrap().observations.len(), 2);
    }

    #[test]
    fn evidence_is_grouped_by_reference_and_remembers_producer() {
        let mut curator = ReferenceCurator::new();
        let id = curator.observe(b"ACGT", observation("a")).unwrap();
        let ann = Annotation::new(Producer { name: "STAR".into(), version: Some("2.7".into()), parameters: Some("--x y".into()) }, "chr1", 10, 14, Strand::Forward, 4);
        let evidence = curator.add_annotation(&id, "GRCm39_M39", ann).unwrap();
        let stored = &curator.candidate(&id).unwrap().annotations["GRCm39_M39"][0];
        assert_eq!(stored.evidence_id, evidence);
        assert_eq!(stored.producer.name, "STAR");
    }

    #[test]
    fn fasta_export_can_select_unresolved_candidates() {
        let mut curator = ReferenceCurator::new();
        let unresolved = curator.observe(b"AAAA", observation("a")).unwrap();
        let resolved = curator.observe(b"CCCC", observation("b")).unwrap();
        curator.resolve(&resolved, ResolvedReference { sequence: b"CCCT".to_vec(), source: "manual".into(), evidence_ids: vec![] }).unwrap();

        let path = std::env::temp_dir().join(format!("reference-curator-fasta-{}.fa", std::process::id()));
        curator.export_fasta(&path, CandidateFilter::UnresolvedOnly).unwrap();
        let fasta = fs::read_to_string(&path).unwrap();
        fs::remove_file(path).ok();

        assert!(fasta.contains(&format!(">{unresolved}")));
        assert!(fasta.contains("AAAA"));
        assert!(!fasta.contains(&format!(">{resolved}")));
        assert!(!fasta.contains("CCCT"));
    }

    #[test]
    fn unresolved_is_a_valid_persistent_state() {
        let mut curator = ReferenceCurator::new();
        let id = curator.observe(b"AACCGGTT", observation("a")).unwrap();
        let path = std::env::temp_dir().join(format!("reference-curator-{}.bin", std::process::id()));
        curator.save(&path).unwrap();
        let loaded = ReferenceCurator::open(&path).unwrap();
        fs::remove_file(path).ok();
        assert!(loaded.candidate(&id).unwrap().resolved.is_none());
        assert_eq!(loaded.unresolved().count(), 1);
    }
}
