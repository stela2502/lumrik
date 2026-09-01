use int_to_str::IntToStr;
use mapping_info::MappingInfo;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::{FeatureEntry, TagEntry};

const K: usize = 8;
const TABLE_SIZE: usize = 1 << 16;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Slot {
    Empty,
    Hit(TagEntry),
    Duplicate,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MapStatus {
    Hit {
        feature_id: u64,
        feature_index: usize,
        start: isize,
        hits: u32,
    },
    NoHit,
    Tie {
        hits: u32,
        feature_ids: Vec<u64>,
    },
}

#[derive(Debug, Clone)]
pub struct FastTagMapper {
    table: Vec<Slot>,
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
            table: vec![Slot::Empty; TABLE_SIZE],
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

    pub fn slot(&self, encoded_8mer: u16) -> Slot {
        self.table[encoded_8mer as usize]
    }

    pub fn load_fasta<P: AsRef<Path>>(&mut self, path: P) -> std::io::Result<usize> {
        self.load_fasta_as(path, "Antibody Capture")
    }

    /// Load FASTA features and assign every record the supplied feature type.
    ///
    /// Nelrune uses this to keep short-feature classes separate: for example
    /// `hto.fa` is loaded with feature type `hto`, while the FASTA record names
    /// remain the actual feature names.
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

    /// Add a feature/sample/FASTA record.
    ///
    /// Returns the internal feature index. The external Scdata feature id is
    /// `feature.id`.
    pub fn add_feature(&mut self, seq: &[u8], feature: FeatureEntry) -> usize {
        let feature_index = self.features.len();
        self.features.push(feature);

        for tag_pos in 0..=seq.len().saturating_sub(K) {
            let Some(kmer) = encode_8mer_with_int_to_str(&seq[tag_pos..tag_pos + K]) else {
                continue;
            };

            let new_entry = TagEntry {
                feature_index,
                tag_pos,
            };

            let slot = &mut self.table[kmer as usize];
            match *slot {
                Slot::Empty => *slot = Slot::Hit(new_entry),
                Slot::Hit(old_entry) => {
                    if old_entry != new_entry {
                        *slot = Slot::Duplicate;
                    }
                }
                Slot::Duplicate => {}
            }
        }

        feature_index
    }

    /// Hot API.
    ///
    /// Returns the Scdata-ready feature id only if a unique feature/start pair
    /// surpassed the internal `min_hits` threshold.
    pub fn map_feature_id(&self, seq: &[u8], mapping: &mut MappingInfo) -> Option<u64> {
        match self.map_status(seq, mapping) {
            MapStatus::Hit { feature_id, .. } => Some(feature_id),
            MapStatus::NoHit | MapStatus::Tie { .. } => None,
        }
    }

    /// Debug/status API.
    ///
    /// The public decision rule is the same as `map_feature_id`.
    pub fn map_status(&self, seq: &[u8], mapping: &mut MappingInfo) -> MapStatus {
        mapping.start_ticker();

        let mut votes: HashMap<(usize, isize), u32> = HashMap::new();

        for query_pos in 0..=seq.len().saturating_sub(K) {
            let Some(kmer) = encode_8mer_with_int_to_str(&seq[query_pos..query_pos + K]) else {
                mapping.report("bd_fast_mapper_invalid_8mer");
                continue;
            };

            if let Slot::Hit(entry) = self.table[kmer as usize] {
                let start = query_pos as isize - entry.tag_pos as isize;
                *votes.entry((entry.feature_index, start)).or_insert(0) += 1;
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

    pub fn map_encoded_positions_feature_id<I>(
        &self,
        encoded: I,
        mapping: &mut MappingInfo,
    ) -> Option<u64>
    where
        I: IntoIterator<Item = (usize, u16)>,
    {
        match self.map_encoded_positions_status(encoded, mapping) {
            MapStatus::Hit { feature_id, .. } => Some(feature_id),
            MapStatus::NoHit | MapStatus::Tie { .. } => None,
        }
    }

    pub fn map_encoded_positions_status<I>(
        &self,
        encoded: I,
        mapping: &mut MappingInfo,
    ) -> MapStatus
    where
        I: IntoIterator<Item = (usize, u16)>,
    {
        mapping.start_ticker();

        let mut votes: HashMap<(usize, isize), u32> = HashMap::new();

        for (query_pos, kmer) in encoded {
            if let Slot::Hit(entry) = self.table[kmer as usize] {
                let start = query_pos as isize - entry.tag_pos as isize;
                *votes.entry((entry.feature_index, start)).or_insert(0) += 1;
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

    fn resolve_votes(&self, votes: HashMap<(usize, isize), u32>) -> MapStatus {
        let Some(best_hits) = votes.values().copied().max() else {
            return MapStatus::NoHit;
        };

        if best_hits < self.min_hits {
            return MapStatus::NoHit;
        }

        let mut best: Vec<((usize, isize), u32)> = votes
            .into_iter()
            .filter(|(_, hits)| *hits == best_hits)
            .collect();

        if best.len() != 1 {
            best.sort_by_key(|((feature_index, start), _)| (*feature_index, *start));

            let feature_ids = best
                .into_iter()
                .map(|((feature_index, _), _)| self.features[feature_index].id)
                .collect();

            return MapStatus::Tie {
                hits: best_hits,
                feature_ids,
            };
        }

        let ((feature_index, start), hits) = best.pop().unwrap();
        let feature_id = self.features[feature_index].id;

        MapStatus::Hit {
            feature_id,
            feature_index,
            start,
            hits,
        }
    }
}

pub fn encode_8mer_with_int_to_str(seq: &[u8]) -> Option<u16> {
    if seq.len() != K {
        return None;
    }

    if !seq.iter().all(|b| {
        matches!(
            b,
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'N' | b'n'
        )
    }) {
        return None;
    }

    Some(IntToStr::new(seq).into_u16())
}

pub fn encode_seq_positions_with_int_to_str(seq: &[u8]) -> Vec<(usize, u16)> {
    let mut ret = Vec::new();

    for pos in 0..=seq.len().saturating_sub(K) {
        if let Some(kmer) = encode_8mer_with_int_to_str(&seq[pos..pos + K]) {
            ret.push((pos, kmer));
        }
    }

    ret
}
