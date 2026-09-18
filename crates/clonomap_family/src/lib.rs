//! ClonoMap: scalable inference of clonal structure in large AIRR-seq clones.
//!
//! ClonoMap infers subclonal structure within B-cell receptor (BCR) clones
//! using mutation-aware geometric clustering rather than full phylogenetic
//! reconstruction. It is designed to scale to clones containing tens of
//! thousands of sequences and integrates seamlessly with Change-O outputs.

mod clone_data;
mod encoder;
mod family;
mod pca;
mod tree;

mod reference_models;

pub use clone_data::CloneData;
pub use encoder::OneHotEncoder;
pub use reference_models::{ReferenceModelEntry, ReferenceModels};
pub use family::{align_fragment, AlignedCell, CellReceptor, Family, FamilyConfig, FamilyMutationReport, IndelEvent, LightClone, LightMember, MutationMeasurement, Receptor};
pub use pca::PcaModel;
pub use tree::{MstTree, rooted_categorical_hex, rooted_continuous_hex};

use ndarray::Array2;
use std::error::Error;

/// Combined PCA + MST pipeline structure.
pub struct ClonoMap {
    pub encoder: OneHotEncoder,
    pub pca: PcaModel,
    pub tree: MstTree,
}

impl ClonoMap {
    /// Build PCA + MST from raw sequences.
    pub fn new(seqs: Vec<String>, k: usize, aa_based: bool) -> Result<Self, Box<dyn Error>> {
        /*let mut map: HashMap<String, usize> = HashMap::new();

        for s in seqs {
            *map.entry(s).or_insert(0) += 1;
        }

        // Unique sequences
        let unique: Vec<String> = map.keys().cloned().collect();
        // Counts for each unique seq
        let counts: Vec<usize> = unique.iter().map(|s| map[s]).collect();
        */

        // Encode sequences numerically
        let mut encoder = OneHotEncoder::new();
        println!("ClonoMap::new - I got {} sequences", seqs.len());

        let encoded = encoder.encode_relative(&seqs, aa_based)?;

        if aa_based {
            println!(
                "      and I encoded {} unique amino acid respresentations ({} hot one columns)",
                encoded.nrows(),
                encoded.ncols()
            );
        } else {
            println!(
                "      and I encoded {} changed DNA respresentations ({} hot one columns)",
                encoded.nrows(),
                encoded.ncols()
            );
        }

        // Fit PCA
        let mut pca = PcaModel::new(k);
        pca.fit_transform(&encoded)?;

        // Build tree in PCA space
        let tree = MstTree::build(&pca);

        Ok(Self { encoder, pca, tree })
    }

    /// Build PCA + MST from a caller-supplied feature matrix.
    ///
    /// This is the escape hatch for paired-receptor models: callers can combine
    /// HC mutation features, LC mutation features and categorical structural
    /// receptor features before PCA instead of trying to annotate an HC-only
    /// manifold after the fact.
    ///
    /// `row_labels` are retained in `encoder.sequences` so existing plotting and
    /// TSV code remains row-aligned. They must be unique and have one entry per
    /// feature-matrix row.
    pub fn from_feature_matrix(
        row_labels: Vec<String>,
        features: Array2<f32>,
        k: usize,
    ) -> Result<Self, Box<dyn Error>> {
        if row_labels.is_empty() {
            return Err("No feature rows provided".into());
        }
        if features.nrows() != row_labels.len() {
            return Err(format!(
                "Feature row mismatch: {} labels for {} matrix rows",
                row_labels.len(),
                features.nrows()
            )
            .into());
        }
        if features.ncols() < 2 {
            return Err("Feature matrix needs at least two columns".into());
        }

        let mut encoder = OneHotEncoder::new();
        encoder.set_external_states(row_labels, &features)?;

        let mut pca = PcaModel::new(k);
        pca.fit_transform(&features)?;
        let tree = MstTree::build(&pca);

        Ok(Self { encoder, pca, tree })
    }

    /// Build one full-feature-space MST per biological subgroup and attach each
    /// subgroup explicitly to a shared root. PCA is still calculated for exported
    /// coordinates/visualization, but it does not determine tree edges.
    ///
    /// `groups` is row-aligned with `features`; the root row is ignored. Within
    /// each group all feature columns participate in the MST. The group's entry
    /// row is the member with the smallest Euclidean distance to the root in the
    /// first `root_distance_cols` columns (HC coordinates in Valkyrn).
    pub fn from_grouped_feature_matrix(
        row_labels: Vec<String>,
        features: Array2<f32>,
        groups: Vec<String>,
        root_row: usize,
        root_distance_cols: usize,
        k: usize,
    ) -> Result<Self, Box<dyn Error>> {
        if row_labels.is_empty() || features.nrows() != row_labels.len() || groups.len() != row_labels.len() {
            return Err("Grouped feature rows, labels and groups must be non-empty and row-aligned".into());
        }
        if root_row >= features.nrows() || root_distance_cols == 0 || root_distance_cols > features.ncols() {
            return Err("Invalid grouped-tree root or HC distance width".into());
        }

        let mut encoder = OneHotEncoder::new();
        encoder.set_external_states(row_labels, &features)?;
        let mut pca = PcaModel::new(k);
        pca.fit_transform(&features)?;

        let mut by_group = std::collections::BTreeMap::<String, Vec<usize>>::new();
        for (row, group) in groups.into_iter().enumerate() {
            if row != root_row {
                by_group.entry(group).or_default().push(row);
            }
        }

        let mut edges = Vec::with_capacity(features.nrows().saturating_sub(1));
        for rows in by_group.into_values() {
            edges.extend(MstTree::build_feature_rows(&features, &rows).edges);
            let entry = rows.iter().copied().min_by(|&a, &b| {
                let da = features.row(a).iter().take(root_distance_cols).map(|x| x * x).sum::<f32>();
                let db = features.row(b).iter().take(root_distance_cols).map(|x| x * x).sum::<f32>();
                da.total_cmp(&db).then_with(|| a.cmp(&b))
            }).expect("non-empty grouped rows");
            let root_distance = features.row(entry).iter().take(root_distance_cols).map(|x| x * x).sum::<f32>().sqrt();
            edges.push((root_row, entry, root_distance));
        }
        let tree = MstTree { edges };
        Ok(Self { encoder, pca, tree })
    }

    /// PCA coordinates accessor
    pub fn coords(&self) -> &Array2<f32> {
        self.pca.coords()
    }

    /// Tree edge list
    pub fn tree(&self) -> &Vec<(usize, usize, f32)> {
        &self.tree.edges
    }
}
