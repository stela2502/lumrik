//! Fast 8-mer mapper with an explicit locus namespace.
//!
//! The ordinary `FastTagMapper` deliberately treats an 8-mer shared by two
//! features as ambiguous.  That is ideal for sample tags, but not for many
//! cell-specific receptor baits: two cells may legitimately share J-derived
//! sequence.  `FastLocusMapper` keeps one global 8-mer index and filters hits
//! by `locus_id` while mapping.  For sc-vdj the locus is the cell id.

use crate::fast_mapper::encode_8mer_with_int_to_str;
use crate::{FeatureEntry, MapStatus};

const TABLE_SIZE: usize = 1 << 16;
const NO_ENTRY: u32 = u32::MAX;

#[derive(Debug, Clone, Copy)]
struct LocusTagEntry {
    feature_index: u32,
    tag_pos: u32,
    next: u32,
}

#[derive(Debug, Clone)]
struct LocusFeature {
    locus_ids: Vec<u64>,
    feature: FeatureEntry,
}

#[derive(Debug, Clone)]
pub struct FastLocusMapper {
    heads: Vec<u32>,
    entries: Vec<LocusTagEntry>,
    features: Vec<LocusFeature>,
    min_hits: u32,
}

impl Default for FastLocusMapper {
    fn default() -> Self {
        Self::new()
    }
}

impl FastLocusMapper {
    pub fn new() -> Self {
        Self {
            heads: vec![NO_ENTRY; TABLE_SIZE],
            entries: Vec::new(),
            features: Vec::new(),
            min_hits: 4,
        }
    }

    pub fn with_min_hits(mut self, min_hits: u32) -> Self {
        self.min_hits = min_hits;
        self
    }

    pub fn feature_count(&self) -> usize {
        self.features.len()
    }

    pub fn feature(&self, feature_index: usize) -> Option<&FeatureEntry> {
        self.features.get(feature_index).map(|x| &x.feature)
    }

    pub fn locus_ids(&self, feature_index: usize) -> Option<&[u64]> {
        self.features
            .get(feature_index)
            .map(|x| x.locus_ids.as_slice())
    }

    /// Add another logical locus alias to an existing sequence feature.
    ///
    /// This is useful when many cells share the same clonotype: the bait
    /// sequence is indexed once while every cell carrying that receptor can
    /// still query it independently.
    pub fn add_locus_alias(&mut self, feature_index: usize, locus_id: u64) {
        let Some(feature) = self.features.get_mut(feature_index) else {
            return;
        };
        match feature.locus_ids.binary_search(&locus_id) {
            Ok(_) => {}
            Err(pos) => feature.locus_ids.insert(pos, locus_id),
        }
    }

    pub fn add_feature(&mut self, locus_id: u64, seq: &[u8], feature: FeatureEntry) -> usize {
        let feature_index = self.features.len();
        assert!(
            feature_index <= u32::MAX as usize,
            "too many locus-mapper features"
        );
        self.features.push(LocusFeature {
            locus_ids: vec![locus_id],
            feature,
        });

        for tag_pos in 0..=seq.len().saturating_sub(8) {
            let Some(kmer) = encode_8mer_with_int_to_str(&seq[tag_pos..tag_pos + 8]) else {
                continue;
            };
            assert!(tag_pos <= u32::MAX as usize, "feature sequence is too long");
            assert!(
                self.entries.len() < u32::MAX as usize,
                "too many locus-mapper 8-mers"
            );

            let slot = kmer as usize;
            let entry_index = self.entries.len() as u32;
            self.entries.push(LocusTagEntry {
                feature_index: feature_index as u32,
                tag_pos: tag_pos as u32,
                next: self.heads[slot],
            });
            self.heads[slot] = entry_index;
        }

        feature_index
    }

    pub fn map_status(&self, locus_id: u64, seq: &[u8]) -> MapStatus {
        // Most locus queries have one or two candidate features. Keep votes in
        // a compact vector so no hash table is allocated for every BAM read.
        let mut votes = Vec::<((usize, isize), u32)>::new();

        for query_pos in 0..=seq.len().saturating_sub(8) {
            let Some(kmer) = encode_8mer_with_int_to_str(&seq[query_pos..query_pos + 8]) else {
                continue;
            };

            let mut entry_index = self.heads[kmer as usize];
            while entry_index != NO_ENTRY {
                let entry = self.entries[entry_index as usize];
                let feature_index = entry.feature_index as usize;
                if self.features[feature_index]
                    .locus_ids
                    .binary_search(&locus_id)
                    .is_ok()
                {
                    let key = (feature_index, query_pos as isize - entry.tag_pos as isize);
                    if let Some((_, count)) =
                        votes.iter_mut().find(|(candidate, _)| *candidate == key)
                    {
                        *count += 1;
                    } else {
                        votes.push((key, 1));
                    }
                }
                entry_index = entry.next;
            }
        }

        self.resolve_votes(votes)
    }

    fn resolve_votes(&self, votes: Vec<((usize, isize), u32)>) -> MapStatus {
        let Some(best_hits) = votes.iter().map(|(_, hits)| *hits).max() else {
            return MapStatus::NoHit;
        };
        if best_hits < self.min_hits {
            return MapStatus::NoHit;
        }

        let mut best: Vec<_> = votes
            .into_iter()
            .filter(|(_, hits)| *hits == best_hits)
            .collect();
        if best.len() != 1 {
            best.sort_by_key(|((feature_index, start), _)| (*feature_index, *start));
            return MapStatus::Tie {
                hits: best_hits,
                feature_ids: best
                    .into_iter()
                    .map(|((feature_index, _), _)| self.features[feature_index].feature.id)
                    .collect(),
            };
        }

        let ((feature_index, start), hits) = best.pop().unwrap();
        MapStatus::Hit {
            feature_id: self.features[feature_index].feature.id,
            feature_index,
            start,
            hits,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_sequence_is_resolved_inside_locus() {
        let mut mapper = FastLocusMapper::new().with_min_hits(3);
        mapper.add_feature(
            11,
            b"AACCGGTTAACCGGTT",
            FeatureEntry::new(1, "cell11", "test"),
        );
        mapper.add_feature(
            22,
            b"AACCGGTTAACCGGTT",
            FeatureEntry::new(2, "cell22", "test"),
        );

        assert!(matches!(
            mapper.map_status(11, b"AACCGGTTAACCGGTT"),
            MapStatus::Hit { feature_id: 1, .. }
        ));
        assert!(matches!(
            mapper.map_status(22, b"AACCGGTTAACCGGTT"),
            MapStatus::Hit { feature_id: 2, .. }
        ));
        assert_eq!(mapper.map_status(33, b"AACCGGTTAACCGGTT"), MapStatus::NoHit);
    }
    #[test]
    fn one_sequence_can_have_multiple_locus_aliases() {
        let mut mapper = FastLocusMapper::new().with_min_hits(3);
        let feature = mapper.add_feature(
            11,
            b"AACCGGTTAACCGGTT",
            FeatureEntry::new(1, "clone", "test"),
        );
        mapper.add_locus_alias(feature, 22);

        assert!(matches!(
            mapper.map_status(11, b"AACCGGTTAACCGGTT"),
            MapStatus::Hit { feature_id: 1, .. }
        ));
        assert!(matches!(
            mapper.map_status(22, b"AACCGGTTAACCGGTT"),
            MapStatus::Hit { feature_id: 1, .. }
        ));
        assert_eq!(mapper.map_status(33, b"AACCGGTTAACCGGTT"), MapStatus::NoHit);
    }
}
