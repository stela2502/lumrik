use std::collections::HashMap;

use fast_tag_mapper::FeatureEntry;
use scdata::FeatureIndex;

#[derive(Debug, Clone)]
pub struct OwnedFeatureIndex {
    features: Vec<FeatureEntry>,
    name_to_id: HashMap<String, u64>,
    id_to_index: HashMap<u64, usize>,
    ordered_ids: Vec<u64>,
}

impl OwnedFeatureIndex {
    pub fn new(features: Vec<FeatureEntry>) -> Self {
        let mut name_to_id = HashMap::new();
        let mut id_to_index = HashMap::new();
        let mut ordered_ids = Vec::new();
        for (i, feature) in features.iter().enumerate() {
            name_to_id.insert(feature.name.clone(), feature.id);
            id_to_index.insert(feature.id, i);
            ordered_ids.push(feature.id);
        }
        Self {
            features,
            name_to_id,
            id_to_index,
            ordered_ids,
        }
    }

    fn feature_by_id(&self, id: u64) -> Option<&FeatureEntry> {
        self.id_to_index
            .get(&id)
            .and_then(|i| self.features.get(*i))
    }

    pub fn split_by_feature_type(&self) -> HashMap<String, OwnedFeatureIndex> {
        let mut groups: HashMap<String, Vec<FeatureEntry>> = HashMap::new();
        for feature in &self.features {
            groups
                .entry(feature.feature_type.clone())
                .or_default()
                .push(feature.clone());
        }
        groups.into_iter().map(|(k, v)| (k, Self::new(v))).collect()
    }
}

impl FeatureIndex for OwnedFeatureIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        self.feature_by_id(feature_id)
            .map(|x| x.name.as_str())
            .unwrap_or("NA")
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        self.ordered_ids.clone()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let f = self.feature_by_id(feature_id).unwrap_or_else(|| {
            panic!("Feature id {feature_id} was not found in the feature index!")
        });
        format!("{}\t{}\t{}", f.name, f.name, f.feature_type)
    }
}
