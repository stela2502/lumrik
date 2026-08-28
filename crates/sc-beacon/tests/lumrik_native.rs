use std::collections::HashMap;

use mapping_info::MappingInfo;
use sc_beacon::GuideDataset;
use scdata::{FeatureIndex, GeneUmiHash, MatrixValueType, Scdata};

struct SparseIdIndex {
    names: HashMap<u64, String>,
    ids: HashMap<String, u64>,
    ordered: Vec<u64>,
}

impl SparseIdIndex {
    fn new() -> Self {
        let ordered = vec![7, 42];
        let names = [
            (7, "SampleTag07_mm".to_string()),
            (42, "Donor_A_HTO".to_string()),
        ]
        .into_iter()
        .collect::<HashMap<_, _>>();
        let ids = names
            .iter()
            .map(|(id, name)| (name.clone(), *id))
            .collect();
        Self { names, ids, ordered }
    }
}

impl FeatureIndex for SparseIdIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        self.names.get(&feature_id).map(String::as_str).unwrap_or("NA")
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.ids.get(name).copied()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let name = self.feature_name(feature_id);
        format!("{name}\t{name}\tTest")
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        self.ordered.clone()
    }
}

#[test]
fn guide_dataset_uses_feature_index_order_not_dense_feature_ids() {
    let index = SparseIdIndex::new();
    let mut data = Scdata::new(1, MatrixValueType::Integer);
    let mut report = MappingInfo::new(None, 0.0, 0);
    let cell = 123_u64;

    data.try_insert_value(&cell, GeneUmiHash(42, 0), 3.0, &mut report);

    let dataset = GuideDataset::from_scdata(&data, &index, 16).unwrap();

    assert_eq!(dataset.feature_ids, vec![7, 42]);
    assert_eq!(dataset.by_guide[0].len(), 0);
    assert_eq!(dataset.by_guide[1].len(), 1);
    assert_eq!(dataset.by_guide[1][0].guide_id, 1);
    assert_eq!(dataset.feature_id(1), 42);
}
