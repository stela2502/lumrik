use std::collections::HashMap;

use scdata::FeatureIndex;

use fast_tag_mapper::{FastTagMapper, FeatureEntry};

pub struct FastTagFeatureIndex<'a> {
    mapper: &'a FastTagMapper,
    name_to_id: HashMap<String, u64>,
    id_to_feature_index: HashMap<u64, usize>,
    ordered_ids: Vec<u64>,
}

impl<'a> FastTagFeatureIndex<'a> {
    pub fn new(mapper: &'a FastTagMapper) -> Self {
        let mut name_to_id = HashMap::new();
        let mut id_to_feature_index = HashMap::new();
        let mut ordered_ids = Vec::new();

        for (feature_index, feature) in mapper.features().iter().enumerate() {
            name_to_id.insert(feature.name.clone(), feature.id);
            id_to_feature_index.insert(feature.id, feature_index);
            ordered_ids.push(feature.id);
        }

        Self {
            mapper,
            name_to_id,
            id_to_feature_index,
            ordered_ids,
        }
    }

    fn feature_by_id(&self, feature_id: u64) -> Option<&FeatureEntry> {
        self.id_to_feature_index
            .get(&feature_id)
            .and_then(|idx| self.mapper.feature(*idx))
    }
}

impl FeatureIndex for FastTagFeatureIndex<'_> {
    fn feature_name(&self, feature_id: u64) -> &str {
        self.feature_by_id(feature_id)
            .map(|feature| feature.name.as_str())
            .unwrap_or("NA")
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        self.ordered_ids.clone()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let feature = self
            .feature_by_id(feature_id)
            .unwrap_or_else(|| panic!("Feature id {feature_id} was not found in the feature index!"));

        format!("{}	{}	{}", feature.name, feature.name, feature.feature_type)
    }
}
