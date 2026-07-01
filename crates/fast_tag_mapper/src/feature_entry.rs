#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureEntry {
    /// Scdata / matrix feature id.
    ///
    /// Built-in BD sample tags use 1..=12.
    pub id: u64,

    /// Clean FASTA/BD/sample name.
    pub name: String,

    /// 10x feature type.
    pub feature_type: String,
}

impl FeatureEntry {
    pub fn new(id: u64, name: impl Into<String>, feature_type: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            feature_type: feature_type.into(),
        }
    }

    pub fn bd_mouse(id: u64) -> Self {
        Self::new(id, format!("SampleTag{id:02}_mm"), "Antibody Capture")
    }

    pub fn bd_human(id: u64) -> Self {
        Self::new(id, format!("SampleTag{id:02}_hs"), "Antibody Capture")
    }
}
