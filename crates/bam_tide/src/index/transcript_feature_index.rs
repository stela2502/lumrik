//transcript_feature_index.rs

use gtf_splice_index::SpliceIndex;
use scdata::FeatureIndex;
use std::collections::HashMap;

pub struct TranscriptFeatureIndex<'a> {
    idx: &'a SpliceIndex,
    name_to_id: HashMap<String, u64>,
}

impl<'a> TranscriptFeatureIndex<'a> {
    pub fn new(idx: &'a SpliceIndex) -> Self {
        let mut name_to_id = HashMap::new();

        for tx in &idx.transcripts {
            for name in &tx.names {
                name_to_id.entry(name.clone()).or_insert(tx.id as u64);
            }
        }

        Self { idx, name_to_id }
    }
}

impl FeatureIndex for TranscriptFeatureIndex<'_> {
    fn feature_name(&self, feature_id: u64) -> &str {
        self.idx
            .transcript_name(feature_id as usize)
            .unwrap_or("NA")
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.idx.transcripts.len() as u64).collect()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let name = self.feature_name(feature_id);
        format!("{name}\t{name}\tGene Expression")
    }
}
