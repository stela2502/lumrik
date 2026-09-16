//! Exact 16-bp supplemental-feature mapper.
//!
//! The first 8 bp of every 16-bp seed select one of 256 bins plus an 8-bit
//! key inside that bin. The following 8 bp are stored as the exact
//! confirmation word. There is no mismatch correction and no alignment-offset
//! inference: an exact 16-bp seed contributes one vote to its feature.

use int_to_str::IntToStr;
use mapping_info::MappingInfo;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::{FeatureEntry, MapStatus};

const SEED_BASES: usize = 16;
const BIN_COUNT: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SeedEntry {
    /// Low byte of the encoded first 8 bp. The high byte selected the bin.
    prefix_key: u8,
    /// Encoded second 8 bp.
    confirm: u16,
    feature_index: u32,
}

#[derive(Debug, Clone)]
pub struct FastTagMapper {
    bins: [Vec<SeedEntry>; BIN_COUNT],
    features: Vec<FeatureEntry>,
    min_hits: u32,
}

impl Default for FastTagMapper {
    fn default() -> Self {
        Self::new()
    }
}

impl FastTagMapper {
    pub fn new() -> Self {
        Self {
            bins: std::array::from_fn(|_| Vec::new()),
            features: Vec::new(),
            min_hits: 4,
        }
    }

    pub fn with_min_hits(mut self, min_hits: u32) -> Self {
        self.min_hits = min_hits;
        self
    }

    pub fn min_hits(&self) -> u32 {
        self.min_hits
    }

    pub fn feature_count(&self) -> usize {
        self.features.len()
    }

    pub fn features(&self) -> &[FeatureEntry] {
        &self.features
    }

    pub fn feature(&self, feature_index: usize) -> Option<&FeatureEntry> {
        self.features.get(feature_index)
    }

    pub fn feature_by_id(&self, feature_id: u64) -> Option<&FeatureEntry> {
        self.features.iter().find(|f| f.id == feature_id)
    }

    pub fn indexed_seed_count(&self) -> usize {
        self.bins.iter().map(Vec::len).sum()
    }

    pub fn occupied_bin_count(&self) -> usize {
        self.bins.iter().filter(|bin| !bin.is_empty()).count()
    }

    pub fn load_fasta<P: AsRef<Path>>(&mut self, path: P) -> std::io::Result<usize> {
        self.load_fasta_as(path, "Antibody Capture")
    }

    pub fn load_fasta_as<P, S>(&mut self, path: P, feature_type: S) -> std::io::Result<usize>
    where
        P: AsRef<Path>,
        S: AsRef<str>,
    {
        let reader = BufReader::new(File::open(path)?);
        let feature_type = feature_type.as_ref();

        let mut name: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();
        let mut added = 0usize;

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() {
                continue;
            }

            if let Some(header) = line.strip_prefix('>') {
                if let Some(old_name) = name.take() {
                    added += self.add_loaded_fasta_record(old_name, &seq, feature_type);
                    seq.clear();
                }

                name = Some(
                    header
                        .split_whitespace()
                        .next()
                        .unwrap_or(header)
                        .to_string(),
                );
            } else {
                seq.extend_from_slice(line.as_bytes());
            }
        }

        if let Some(old_name) = name {
            added += self.add_loaded_fasta_record(old_name, &seq, feature_type);
        }

        Ok(added)
    }

    fn add_loaded_fasta_record(&mut self, name: String, seq: &[u8], feature_type: &str) -> usize {
        let feature_id = self.features().iter().map(|f| f.id).max().unwrap_or(0) + 1;
        self.add_feature(seq, FeatureEntry::new(feature_id, name, feature_type));
        1
    }

    /// Add a feature/sample/FASTA record and index all exact 16-bp seeds.
    ///
    /// Repeated copies of the same 16-mer inside one feature are indexed only
    /// once. A 16-mer shared by different features is retained for each feature
    /// so later seed evidence can resolve the call or report a tie.
    pub fn add_feature(&mut self, seq: &[u8], feature: FeatureEntry) -> usize {
        let feature_index = self.features.len();
        self.features.push(feature);

        if seq.len() < SEED_BASES {
            return feature_index;
        }

        for seed in Rolling16::new(seq) {
            let (bin_index, prefix_key, confirm) = split_seed(seed);
            let bin = &mut self.bins[bin_index];

            if bin.iter().any(|entry| {
                entry.prefix_key == prefix_key
                    && entry.confirm == confirm
                    && entry.feature_index == feature_index as u32
            }) {
                continue;
            }

            bin.push(SeedEntry {
                prefix_key,
                confirm,
                feature_index: feature_index as u32,
            });
        }

        feature_index
    }

    /// Hot API: return the Scdata feature id for a unique feature with at least
    /// `min_hits` exact 16-bp seed matches.
    pub fn map_feature_id(&self, seq: &[u8], mapping: &mut MappingInfo) -> Option<u64> {
        match self.map_status(seq, mapping) {
            MapStatus::Hit { feature_id, .. } => Some(feature_id),
            MapStatus::NoHit | MapStatus::Tie { .. } => None,
        }
    }

    pub fn map_status(&self, seq: &[u8], mapping: &mut MappingInfo) -> MapStatus {
        mapping.start_ticker();

        // Negative reads allocate nothing. This Vec allocates only after the
        // first exact 16-bp seed hit, which is the rare path for supplemental
        // feature mapping.
        let mut votes: Vec<(u32, u32)> = Vec::new();

        for seed in Rolling16::new(seq) {
            let (bin_index, prefix_key, confirm) = split_seed(seed);

            for entry in &self.bins[bin_index] {
                if entry.prefix_key != prefix_key || entry.confirm != confirm {
                    continue;
                }

                if let Some((_, hits)) = votes
                    .iter_mut()
                    .find(|(feature_index, _)| *feature_index == entry.feature_index)
                {
                    *hits += 1;
                } else {
                    votes.push((entry.feature_index, 1));
                }
            }
        }

        let status = self.resolve_votes(votes);

        match &status {
            MapStatus::Hit { .. } => mapping.report("bd_fast_mapper_hit"),
            MapStatus::NoHit => mapping.report("bd_fast_mapper_no_hit"),
            MapStatus::Tie { .. } => mapping.report("bd_fast_mapper_tie"),
        }

        mapping.stop_single_processor_time();
        status
    }

    fn resolve_votes(&self, votes: Vec<(u32, u32)>) -> MapStatus {
        let Some(best_hits) = votes.iter().map(|(_, hits)| *hits).max() else {
            return MapStatus::NoHit;
        };

        if best_hits < self.min_hits {
            return MapStatus::NoHit;
        }

        let best = votes
            .into_iter()
            .filter(|(_, hits)| *hits == best_hits)
            .map(|(feature_index, _)| feature_index as usize)
            .collect::<Vec<_>>();

        if best.len() != 1 {
            let mut feature_ids = best
                .into_iter()
                .map(|feature_index| self.features[feature_index].id)
                .collect::<Vec<_>>();
            feature_ids.sort_unstable();

            return MapStatus::Tie {
                hits: best_hits,
                feature_ids,
            };
        }

        let feature_index = best[0];
        MapStatus::Hit {
            feature_id: self.features[feature_index].id,
            feature_index,
            hits: best_hits,
        }
    }
}

// Kept outside FastTagMapper so the hot representation and sequence encoding
// remain separate modules/concepts rather than accumulating in lib.rs.
fn split_seed(seed: u32) -> (usize, u8, u16) {
    let first_8 = (seed >> 16) as u16;
    ((first_8 >> 8) as usize, first_8 as u8, seed as u16)
}

struct Rolling16<'a> {
    seq: &'a [u8],
    pos: usize,
    word: u32,
    valid_bases: usize,
}

impl<'a> Rolling16<'a> {
    fn new(seq: &'a [u8]) -> Self {
        Self {
            seq,
            pos: 0,
            word: 0,
            valid_bases: 0,
        }
    }
}

impl Iterator for Rolling16<'_> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let base = self.seq[self.pos];
            self.pos += 1;

            // Keep IntToStr as the canonical DNA -> 2-bit implementation.
            // N is deliberately not accepted as an exact seed base: IntToStr
            // maps N to A for compact storage, but this mapper must not turn an
            // ambiguous base into an exact 16-bp match.
            let bits = match base {
                b'N' | b'n' => None,
                _ => IntToStr::encode_binary(base).ok().map(u32::from),
            };

            match bits {
                Some(bits) => {
                    self.word = (self.word << 2) | bits;
                    self.valid_bases += 1;
                    if self.valid_bases >= SEED_BASES {
                        return Some(self.word);
                    }
                }
                None => {
                    self.word = 0;
                    self.valid_bases = 0;
                }
            }
        }
        None
    }
}

/// Shared by FastLocusMapper. This is intentionally only an encoder now; the
/// old FastTagMapper 8-mer table and positional voting API are gone.
pub(crate) fn encode_8mer(seq: &[u8]) -> Option<u16> {
    if seq.len() != 8 {
        return None;
    }

    let mut word = 0u16;
    for &base in seq {
        if matches!(base, b'N' | b'n') {
            return None;
        }
        word = (word << 2) | u16::from(IntToStr::encode_binary(base).ok()?);
    }
    Some(word)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rolling_16_matches_expected_encoding() {
        let seeds = Rolling16::new(b"ACGTACGTACGTACGTA").collect::<Vec<_>>();
        assert_eq!(seeds.len(), 2);
        assert_ne!(seeds[0], seeds[1]);
    }

    #[test]
    fn invalid_base_resets_window() {
        assert_eq!(
            Rolling16::new(b"AAAAAAAAAAAAAAAANAAAAAAAAAAAAAAAA").count(),
            2
        );
    }
}
