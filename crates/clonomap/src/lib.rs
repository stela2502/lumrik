//! ClonoMap: scalable inference of clonal structure in large AIRR-seq clones.
//!
//! ClonoMap infers subclonal structure within B-cell receptor (BCR) clones
//! using mutation-aware geometric clustering rather than full phylogenetic
//! reconstruction. It is designed to scale to clones containing tens of
//! thousands of sequences and integrates seamlessly with Change-O outputs.

mod clone_data;
mod encoder;
mod pca;
mod tree;

pub use clone_data::CloneData;
pub use encoder::OneHotEncoder;
pub use pca::PcaModel;
pub use tree::MstTree;

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

    /// PCA coordinates accessor
    pub fn coords(&self) -> &Array2<f32> {
        self.pca.coords()
    }

    /// Tree edge list
    pub fn tree(&self) -> &Vec<(usize, usize, f32)> {
        &self.tree.edges
    }
}
