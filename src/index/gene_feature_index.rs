// gene_feature_index.rs

use gtf_splice_index::SpliceIndex;
use scdata::FeatureIndex;
use std::collections::HashMap;

pub struct GeneFeatureIndex<'a> {
    idx: &'a SpliceIndex,
    name_to_id: HashMap<String, u64>,
}

impl<'a> GeneFeatureIndex<'a> {
    pub fn new(idx: &'a SpliceIndex) -> Self {
        let mut name_to_id = HashMap::new();

        for gene in &idx.genes {
            for name in &gene.names {
                name_to_id.entry(name.clone()).or_insert(gene.id as u64);
            }
        }

        Self { idx, name_to_id }
    }
}

impl FeatureIndex for GeneFeatureIndex<'_> {
    fn feature_name(&self, feature_id: u64) -> &str {
        self.idx.gene_name(feature_id as usize).unwrap_or("NA")
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.idx.genes.len() as u64).collect()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let name = self.feature_name(feature_id);
        format!("{name}\t{name}\tGene Expression")
    }
}
