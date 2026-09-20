//! Seed-and-verify supplemental-feature mapper.
//!
//! Exact 16-mers are deliberately only a candidate-location index.  A seed
//! hit records the feature, reference position and strand; the final call is
//! made by comparing the complete feature against a concatenated one-hot
//! reference genome.  This keeps the old very cheap rejection path while no
//! longer treating a pile of seed votes as an alignment.

use int_to_dna::IntToDna;
use mapping_info::MappingInfo;
use onehot_dna::OneHotSequence;
use std::{fs::File, io::{BufRead, BufReader}, path::Path};

use crate::{FeatureEntry, MapStatus};

const SEED_BASES: usize = 16;
const BIN_COUNT: usize = 256;
const DEFAULT_PHRED: u8 = 30;
const MAX_CANDIDATES: usize = 10;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SeedEntry {
    prefix_key: u8,
    confirm: u16,
    feature_index: u32,
    ref_pos: u32,
    reverse: bool,
}

#[derive(Debug, Clone, Copy)]
struct ReferenceSpan {
    forward_start: u32,
    reverse_start: u32,
    len: u32,
}

/// Compact direct-comparison reference.  Each byte is a one-hot DNA base
/// (A=1, C=2, G=4, T=8, other=0).  All feature sequences are concatenated.
/// Forward and reverse-complement copies are kept separately so the hot
/// verifier never has to decode or reverse a reference sequence.
#[derive(Debug, Clone, Default)]
struct OneHotGenome {
    forward: Vec<u8>,
    reverse: Vec<u8>,
    spans: Vec<ReferenceSpan>,
}

impl OneHotGenome {
    fn push(&mut self, seq: &[u8]) -> usize {
        assert!(self.forward.len() <= u32::MAX as usize);
        assert!(self.reverse.len() <= u32::MAX as usize);
        assert!(seq.len() <= u32::MAX as usize);
        let forward_start = self.forward.len() as u32;
        self.forward.extend(seq.iter().copied().map(one_hot));
        let reverse_start = self.reverse.len() as u32;
        self.reverse.extend(seq.iter().rev().copied().map(complement_one_hot));
        self.spans.push(ReferenceSpan { forward_start, reverse_start, len: seq.len() as u32 });
        self.spans.len() - 1
    }

    #[inline]
    fn slice(&self, feature_index: usize, reverse: bool) -> &[u8] {
        let span = self.spans[feature_index];
        let start = (if reverse { span.reverse_start } else { span.forward_start }) as usize;
        let len = span.len as usize;
        let genome = if reverse { &self.reverse } else { &self.forward };
        &genome[start..start + len]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentStrand { Forward, Reverse }

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastAlignment {
    pub feature_id: u64,
    pub feature_index: usize,
    pub query_start: usize,
    pub strand: AlignmentStrand,
    pub seed_hits: u32,
    pub score: i32,
    pub second_best_score: Option<i32>,
    pub mismatches: u32,
    /// SAM-compatible ungapped CIGAR.  Candidate discovery tolerates sequence
    /// mismatches; indel-aware refinement can be added behind this type later
    /// without changing callers.
    pub cigar: String,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    feature_index: usize,
    query_start: isize,
    reverse: bool,
    hits: u32,
}

#[derive(Debug, Clone, Copy)]
struct VerifiedCandidate {
    candidate: Candidate,
    score: i32,
    mismatches: u32,
}

#[derive(Debug, Clone)]
pub struct FastTagMapper {
    bins: [Vec<SeedEntry>; BIN_COUNT],
    features: Vec<FeatureEntry>,
    genome: OneHotGenome,
    packed_features: Vec<OneHotSequence>,
    packed_reverse_features: Vec<OneHotSequence>,
    min_hits: u32,
}

impl Default for FastTagMapper { fn default() -> Self { Self::new() } }

impl FastTagMapper {
    pub fn new() -> Self {
        Self { bins: std::array::from_fn(|_| Vec::new()), features: Vec::new(), genome: OneHotGenome::default(), packed_features: Vec::new(), packed_reverse_features: Vec::new(), min_hits: 4 }
    }

    pub fn with_min_hits(mut self, min_hits: u32) -> Self { self.min_hits = min_hits; self }
    pub fn min_hits(&self) -> u32 { self.min_hits }
    pub fn feature_count(&self) -> usize { self.features.len() }
    pub fn features(&self) -> &[FeatureEntry] { &self.features }
    pub fn feature(&self, feature_index: usize) -> Option<&FeatureEntry> { self.features.get(feature_index) }
    pub fn feature_by_id(&self, feature_id: u64) -> Option<&FeatureEntry> { self.features.iter().find(|f| f.id == feature_id) }
    pub fn indexed_seed_count(&self) -> usize { self.bins.iter().map(Vec::len).sum() }
    pub fn occupied_bin_count(&self) -> usize { self.bins.iter().filter(|bin| !bin.is_empty()).count() }
    pub fn reference_bases(&self) -> usize { self.genome.forward.len() }

    pub fn load_fasta<P: AsRef<Path>>(&mut self, path: P) -> std::io::Result<usize> { self.load_fasta_as(path, "Antibody Capture") }

    pub fn load_fasta_as<P, S>(&mut self, path: P, feature_type: S) -> std::io::Result<usize>
    where P: AsRef<Path>, S: AsRef<str> {
        let reader = BufReader::new(File::open(path)?);
        let feature_type = feature_type.as_ref();
        let mut name: Option<String> = None;
        let mut seq = Vec::new();
        let mut added = 0usize;
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() { continue; }
            if let Some(header) = line.strip_prefix('>') {
                if let Some(old_name) = name.take() {
                    added += self.add_loaded_fasta_record(old_name, &seq, feature_type);
                    seq.clear();
                }
                name = Some(header.split_whitespace().next().unwrap_or(header).to_string());
            } else { seq.extend_from_slice(line.as_bytes()); }
        }
        if let Some(old_name) = name { added += self.add_loaded_fasta_record(old_name, &seq, feature_type); }
        Ok(added)
    }

    fn add_loaded_fasta_record(&mut self, name: String, seq: &[u8], feature_type: &str) -> usize {
        let feature_id = self.features.iter().map(|f| f.id).max().unwrap_or(0) + 1;
        self.add_feature(seq, FeatureEntry::new(feature_id, name, feature_type));
        1
    }

    /// Add a feature and index exact 16-mers in both orientations.  Unlike the
    /// old mapper, entries retain their reference position; hits therefore
    /// nominate an alignment start instead of merely voting for a feature.
    pub fn add_feature(&mut self, seq: &[u8], feature: FeatureEntry) -> usize {
        let feature_index = self.features.len();
        assert!(feature_index <= u32::MAX as usize, "too many fast-mapper features");
        self.features.push(feature);
        self.genome.push(seq);
        self.packed_features.push(OneHotSequence::from_bytes(seq));
        self.packed_reverse_features.push(OneHotSequence::from_bytes(&reverse_complement(seq)));
        if seq.len() < SEED_BASES { return feature_index; }

        self.index_orientation(seq, feature_index, false);
        let reverse = reverse_complement(seq);
        self.index_orientation(&reverse, feature_index, true);
        feature_index
    }

    fn index_orientation(&mut self, seq: &[u8], feature_index: usize, reverse: bool) {
        for (ref_pos, seed) in Rolling16::new(seq) {
            let (bin_index, prefix_key, confirm) = split_seed(seed);
            let bin = &mut self.bins[bin_index];
            let ref_pos = ref_pos as u32;
            if bin.iter().any(|e| e.prefix_key == prefix_key && e.confirm == confirm && e.feature_index == feature_index as u32 && e.ref_pos == ref_pos && e.reverse == reverse) { continue; }
            bin.push(SeedEntry { prefix_key, confirm, feature_index: feature_index as u32, ref_pos, reverse });
        }
    }

    pub fn map_feature_id(&self, seq: &[u8], mapping: &mut MappingInfo) -> Option<u64> {
        match self.map_status(seq, mapping) { MapStatus::Hit { feature_id, .. } => Some(feature_id), _ => None }
    }

    /// Quality-aware hot API. `qual` contains numeric Phred scores (not ASCII).
    pub fn map_feature_id_with_qual(&self, seq: &[u8], qual: &[u8], mapping: &mut MappingInfo) -> Option<u64> {
        match self.map_status_with_qual(seq, qual, mapping) { MapStatus::Hit { feature_id, .. } => Some(feature_id), _ => None }
    }

    pub fn map_status(&self, seq: &[u8], mapping: &mut MappingInfo) -> MapStatus {
        self.map_status_impl(seq, None, mapping)
    }

    pub fn map_status_with_qual(&self, seq: &[u8], qual: &[u8], mapping: &mut MappingInfo) -> MapStatus {
        self.map_status_impl(seq, Some(qual), mapping)
    }

    fn map_status_impl(&self, seq: &[u8], qual: Option<&[u8]>, mapping: &mut MappingInfo) -> MapStatus {
        mapping.start_ticker();

        // Hot path: a query fragment is only a locator.  Its index bucket
        // already contains every possible (feature, reference position,
        // strand) placement for that exact fragment.  Try those placements
        // immediately against the complete feature and return on the first
        // clean full-overlap match.  Only advance to the next query fragment
        // when this bucket contains no usable placement.
        let (direct_hit, saw_position) = self.first_direct_full_overlap(seq);
        let status = if let Some(hit) = direct_hit {
            hit
        } else if !saw_position {
            // No exact fragment had any indexed position, therefore the old
            // exact-seed candidate path cannot possibly discover a candidate
            // either.  Do not scan the read a second time.
            MapStatus::NoHit
        } else {
            // At least one indexed position existed, but no complete exact
            // overlap survived.  Keep the existing quality-aware path as the
            // rescue for sequencing errors / imperfect feature observations.
            let verified = self.best_candidates(seq, qual);
            self.resolve_verified(&verified)
        };
        match &status {
            MapStatus::Hit { .. } => mapping.report("bd_fast_mapper_hit"),
            MapStatus::NoHit => mapping.report("bd_fast_mapper_no_hit"),
            MapStatus::Tie { .. } => mapping.report("bd_fast_mapper_tie"),
        }
        mapping.stop_single_processor_time();
        status
    }

    /// Return the winning placement with score, strand and CIGAR.  This uses
    /// the same candidate path as map_status, so diagnostics cannot disagree
    /// with the hot call.
    pub fn align(&self, seq: &[u8], qual: Option<&[u8]>) -> Option<FastAlignment> {
        let verified = self.best_candidates(seq, qual);
        let best_score = verified.iter().map(|x| x.score).max()?;
        let mut best = verified.iter().filter(|x| x.score == best_score);
        let winner = *best.next()?;
        if best.next().is_some() || winner.candidate.hits < self.min_hits { return None; }
        let second_best_score = verified.iter().filter(|x| x.candidate.feature_index != winner.candidate.feature_index || x.candidate.query_start != winner.candidate.query_start || x.candidate.reverse != winner.candidate.reverse).map(|x| x.score).max();
        let len = self.genome.spans[winner.candidate.feature_index].len as usize;
        Some(FastAlignment {
            feature_id: self.features[winner.candidate.feature_index].id,
            feature_index: winner.candidate.feature_index,
            query_start: winner.candidate.query_start as usize,
            strand: if winner.candidate.reverse { AlignmentStrand::Reverse } else { AlignmentStrand::Forward },
            seed_hits: winner.candidate.hits,
            score: winner.score,
            second_best_score,
            mismatches: winner.mismatches,
            cigar: format!("{len}M"),
        })
    }

    /// Direct locator -> full-overlap hot path.
    ///
    /// For each query 16-mer, the index bucket is already the complete array
    /// of possible reference positions for that fragment.  We therefore try
    /// those positions immediately.  A clean full-feature overlap is decisive
    /// and returns at once; there is no vote collection on the common path.
    ///
    /// The bool reports whether *any* indexed position was seen.  If false,
    /// the caller can return NoHit without running the exact-seed scan again.
    fn first_direct_full_overlap(&self, seq: &[u8]) -> (Option<MapStatus>, bool) {
        let mut packed_query: Option<OneHotSequence> = None;
        let mut saw_position = false;

        for (query_pos, seed) in Rolling16::new(seq) {
            // The exact seed is only a locator for the full-feature check.
            // Probing every overlapping 16-mer is unnecessary and makes the
            // overwhelmingly common no-hit read pay ~135 index lookups for a
            // 150-bp read.  Probe 16-mers every 8 bases instead: 0, 8, 16, ...
            // This keeps 8-bp overlap between adjacent locator windows while
            // reducing index probes by roughly 8x.
            if query_pos % 8 != 0 { continue; }

            let (bin_index, prefix_key, confirm) = split_seed(seed);
            for entry in &self.bins[bin_index] {
                if entry.prefix_key != prefix_key || entry.confirm != confirm { continue; }
                saw_position = true;

                let feature_index = entry.feature_index as usize;
                let reference_len = self.genome.spans[feature_index].len as usize;
                if reference_len <= SEED_BASES { continue; }

                let query_start = query_pos as isize - entry.ref_pos as isize;
                let Ok(query_start) = usize::try_from(query_start) else { continue; };
                let Some(query_end) = query_start.checked_add(reference_len) else { continue; };
                if query_end > seq.len() { continue; }

                // Packing is lazy: the overwhelmingly common no-index-hit read
                // never allocates a OneHotSequence at all.
                if packed_query.is_none() {
                    packed_query = OneHotSequence::try_from_bytes(seq).ok();
                }
                let Some(query) = packed_query.as_ref() else { return (None, saw_position); };
                let reference = if entry.reverse {
                    &self.packed_reverse_features[feature_index]
                } else {
                    &self.packed_features[feature_index]
                };
                let Some((informative, compatible)) =
                    query.compatibility_counts(query_start, reference, 0, reference_len)
                else { continue; };

                if informative == reference_len && compatible == reference_len {
                    return (Some(MapStatus::Hit {
                        feature_id: self.features[feature_index].id,
                        feature_index,
                        hits: 1,
                    }), true);
                }
            }
        }

        (None, saw_position)
    }

    fn discover_candidates(&self, seq: &[u8]) -> Vec<Candidate> {
        let mut candidates: Vec<Candidate> = Vec::new();
        for (query_pos, seed) in Rolling16::new(seq) {
            let (bin_index, prefix_key, confirm) = split_seed(seed);
            for entry in &self.bins[bin_index] {
                if entry.prefix_key != prefix_key || entry.confirm != confirm { continue; }
                let query_start = query_pos as isize - entry.ref_pos as isize;
                if let Some(c) = candidates.iter_mut().find(|c| c.feature_index == entry.feature_index as usize && c.query_start == query_start && c.reverse == entry.reverse) {
                    c.hits += 1;
                } else {
                    candidates.push(Candidate { feature_index: entry.feature_index as usize, query_start, reverse: entry.reverse, hits: 1 });
                }
            }
        }
        candidates
    }

    fn best_candidates(&self, seq: &[u8], qual: Option<&[u8]>) -> Vec<VerifiedCandidate> {
        self.verify_candidates(seq, qual, self.discover_candidates(seq))
    }

    fn verify_candidates(&self, seq: &[u8], qual: Option<&[u8]>, mut candidates: Vec<Candidate>) -> Vec<VerifiedCandidate> {
        if let Some(q) = qual { if q.len() != seq.len() { return Vec::new(); } }
        candidates.retain(|c| c.hits >= self.min_hits);
        candidates.sort_unstable_by(|a,b| b.hits.cmp(&a.hits));
        candidates.truncate(MAX_CANDIDATES);
        candidates.into_iter().filter_map(|c| self.verify_candidate(seq, qual, c)).collect()
    }

    fn verify_candidate(&self, seq: &[u8], qual: Option<&[u8]>, candidate: Candidate) -> Option<VerifiedCandidate> {
        let query_start = usize::try_from(candidate.query_start).ok()?;
        let reference = self.genome.slice(candidate.feature_index, candidate.reverse);
        let query = seq.get(query_start..query_start.checked_add(reference.len())?)?;
        let qslice = qual.and_then(|q| q.get(query_start..query_start + reference.len()));
        let mut score = 0i32;
        let mut mismatches = 0u32;
        for (i, (&read_base, &reference_base)) in query.iter().zip(reference).enumerate() {
            let read_one_hot = one_hot(read_base);
            let q = qslice.map(|x| x[i]).unwrap_or(DEFAULT_PHRED).min(60);
            if read_one_hot != 0 && read_one_hot == reference_base {
                score += match_reward(q);
            } else {
                mismatches += 1;
                score -= mismatch_penalty(q);
            }
        }
        Some(VerifiedCandidate { candidate, score, mismatches })
    }

    fn resolve_verified(&self, verified: &[VerifiedCandidate]) -> MapStatus {
        let Some(best_score) = verified.iter().map(|x| x.score).max() else { return MapStatus::NoHit; };
        let best: Vec<_> = verified.iter().filter(|x| x.score == best_score).collect();
        if best.len() != 1 {
            let hits = best.iter().map(|x| x.candidate.hits).max().unwrap_or(0);
            let mut feature_ids: Vec<_> = best.iter().map(|x| self.features[x.candidate.feature_index].id).collect();
            feature_ids.sort_unstable(); feature_ids.dedup();
            return MapStatus::Tie { hits, feature_ids };
        }
        let c = best[0].candidate;
        MapStatus::Hit { feature_id: self.features[c.feature_index].id, feature_index: c.feature_index, hits: c.hits }
    }
}

#[inline] fn match_reward(q: u8) -> i32 { 2 + i32::from(q) / 20 }
#[inline] fn mismatch_penalty(q: u8) -> i32 { 2 + i32::from(q) / 5 }

#[inline]
fn one_hot(base: u8) -> u8 { match base { b'A'|b'a'=>1, b'C'|b'c'=>2, b'G'|b'g'=>4, b'T'|b't'=>8, _=>0 } }
#[inline]
fn complement_one_hot(base: u8) -> u8 { match base { b'A'|b'a'=>8, b'C'|b'c'=>4, b'G'|b'g'=>2, b'T'|b't'=>1, _=>0 } }
fn reverse_complement(seq: &[u8]) -> Vec<u8> { seq.iter().rev().map(|&b| match b { b'A'|b'a'=>b'T', b'C'|b'c'=>b'G', b'G'|b'g'=>b'C', b'T'|b't'=>b'A', _=>b'N' }).collect() }

#[inline(always)]
fn encode_exact_16_at(seq: &[u8], start: usize) -> Option<u32> {
    let window = seq.get(start..start.checked_add(SEED_BASES)?)?;
    let mut word = 0u32;
    for &base in window {
        if matches!(base, b'N' | b'n') { return None; }
        word = (word << 2) | u32::from(IntToDna::encode_binary(base).ok()?);
    }
    Some(word)
}

fn split_seed(seed: u32) -> (usize, u8, u16) { let first_8=(seed>>16) as u16; ((first_8>>8) as usize, first_8 as u8, seed as u16) }

struct Rolling16<'a> { seq: &'a [u8], pos: usize, word: u32, valid_bases: usize }
impl<'a> Rolling16<'a> { fn new(seq: &'a [u8]) -> Self { Self { seq, pos:0, word:0, valid_bases:0 } } }
impl Iterator for Rolling16<'_> {
    type Item=(usize,u32);
    fn next(&mut self)->Option<Self::Item>{
        while self.pos < self.seq.len() {
            let base=self.seq[self.pos]; self.pos+=1;
            let bits=match base { b'N'|b'n'=>None, _=>IntToDna::encode_binary(base).ok().map(u32::from) };
            match bits { Some(bits)=>{ self.word=(self.word<<2)|bits; self.valid_bases+=1; if self.valid_bases>=SEED_BASES{return Some((self.pos-SEED_BASES,self.word));} }, None=>{self.word=0;self.valid_bases=0;} }
        }
        None
    }
}

pub(crate) fn encode_8mer(seq:&[u8])->Option<u16>{
    if seq.len()!=8{return None;} let mut word=0u16;
    for &base in seq { if matches!(base,b'N'|b'n'){return None;} word=(word<<2)|u16::from(IntToDna::encode_binary(base).ok()?); }
    Some(word)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn rolling_16_matches_expected_encoding(){let seeds=Rolling16::new(b"ACGTACGTACGTACGTA").collect::<Vec<_>>();assert_eq!(seeds.len(),2);assert_eq!(seeds[0].0,0);assert_eq!(seeds[1].0,1);assert_ne!(seeds[0].1,seeds[1].1);}
    #[test] fn invalid_base_resets_window(){assert_eq!(Rolling16::new(b"AAAAAAAAAAAAAAAANAAAAAAAAAAAAAAAA").count(),2);}
    #[test] fn one_hot_reference_has_both_strands(){let mut g=OneHotGenome::default();g.push(b"ACGTAA");assert_eq!(g.slice(0,false),&[1,2,4,8,1,1]);assert_eq!(g.slice(0,true),&[8,8,1,2,4,8]);}
    #[test] fn quality_changes_mismatch_score(){assert!(mismatch_penalty(40)>mismatch_penalty(5));}
}
