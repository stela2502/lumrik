#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MapStatus {
    Hit {
        feature_id: u64,
        feature_index: usize,
        hits: u32,
    },
    NoHit,
    Tie {
        hits: u32,
        feature_ids: Vec<u64>,
    },
}
