// feature_index.rs

/// Defines the canonical feature space used by `Scdata`.
///
/// # Example
///
/// ```no_run
/// use std::collections::HashMap;
/// use scdata::FeatureIndex;
///
/// struct SimpleIndex {
///     names: Vec<String>,
///     ids: HashMap<String, u64>,
/// }
///
/// impl SimpleIndex {
///     fn new(names: Vec<&str>) -> Self {
///         let names: Vec<String> = names.into_iter().map(|s| s.to_string()).collect();
///         let ids = names
///             .iter()
///             .enumerate()
///             .map(|(i, n)| (n.clone(), i as u64))
///             .collect();
///         Self { names, ids }
///     }
/// }
///
/// impl FeatureIndex for SimpleIndex {
///     fn feature_name(&self, feature_id: u64) -> &str {
///         &self.names[feature_id as usize]
///     }
///
///     fn feature_id(&self, name: &str) -> Option<u64> {
///         self.ids.get(name).copied()
///     }
///
///     fn ordered_feature_ids(&self) -> Vec<u64> {
///         (0..self.names.len() as u64).collect()
///     }
///
///     fn to_10x_feature_line(&self, feature_id: u64) -> String {
///         let name = self.feature_name(feature_id);
///         format!("{name}\t{name}\tGene Expression")
///     }
/// }
///
/// let index = SimpleIndex::new(vec!["GeneA", "GeneB"]);
/// assert_eq!(index.feature_id("GeneA"), Some(0));
/// assert_eq!(index.feature_name(1), "GeneB");
/// ```
pub trait FeatureIndex: Send + Sync {
    /// name to id translation
    fn feature_name(&self, feature_id: u64) -> &str;
    /// id to name translation
    fn feature_id(&self, name: &str) -> Option<u64>;
    /// one line in the features.tsv file
    fn to_10x_feature_line(&self, feature_id: u64) -> String;
    /// Feature ids in canonical export order.
    fn ordered_feature_ids(&self) -> Vec<u64>;
}
