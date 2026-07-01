// thrust_data.rs

use gtf_splice_index::{GeneId, RefBlock, SpliceIndex, Strand, TranscriptId};
use snp_index::{Genome, SnpIndex, VcfReadOptions};

use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

use std::path::Path;

use scdata::cell_data::GeneUmiHash;
use scdata::{MatrixValueType, Scdata, FeatureIndex};

pub struct TruthData {
    pub gene: Scdata,
    pub intron: Scdata,
    pub snp_ref: Scdata,
    pub snp_alt: Scdata,
}

impl TruthData {
    pub fn new() -> Self {
        Self {
            gene: Scdata::new(1, MatrixValueType::Real),
            intron: Scdata::new(1, MatrixValueType::Real),
            snp_ref: Scdata::new(1, MatrixValueType::Real),
            snp_alt: Scdata::new(1, MatrixValueType::Real),
        }
    }
}

impl TruthData {
    pub fn new() -> Self {
        Self {
            gene: Scdata::new(1, MatrixValueType::Real),
            intron: Scdata::new(1, MatrixValueType::Real),
            snp_ref: Scdata::new(1, MatrixValueType::Real),
            snp_alt: Scdata::new(1, MatrixValueType::Real),
        }
    }

    /// Finalize and write all truth matrices to disk.
    ///
    /// This mirrors the behavior of bam-quant output so that
    /// direct comparisons are possible.
    pub fn write<P: AsRef<Path>>(
        &mut self,
        base: P,
        min_cell_counts: usize,
        gene_index: &FeatureIndex,
        snp_index: &FeatureIndex,
    ) -> Result<(), String> {
        let base = base.as_ref();

        let exonic_path = base.join("exonic");
        let intronic_path = base.join("intronic");
        let ref_path = base.join("ref");
        let alt_path = base.join("alt");

        fs::create_dir_all(&exonic_path)
            .map_err(|e| format!("failed to create {:?}: {e}", exonic_path))?;
        fs::create_dir_all(&intronic_path)
            .map_err(|e| format!("failed to create {:?}: {e}", intronic_path))?;
        fs::create_dir_all(&ref_path)
            .map_err(|e| format!("failed to create {:?}: {e}", ref_path))?;
        fs::create_dir_all(&alt_path)
            .map_err(|e| format!("failed to create {:?}: {e}", alt_path))?;        

        // --- finalize ---
        self.gene.finalize_for_export(min_cell_counts, gene_index);

        let cells: std::collections::HashSet<u64> =
            self.gene.export_cell_ids().iter().copied().collect();

        self.gene.finalize_for_cells(cells, gene_index);
        self.intron.finalize_for_cells(cells, gene_index);

        self.snp_alt.finalize_for_cells(cells, snp_index);

        self.snp_ref
            .retain_features(&self.snp_alt.observed_feature_ids());

        self.snp_ref.finalize_for_cells(cells, snp_index);

        // --- write ---
        self.gene
            .write_sparse(&exonic_path, gene_index)
            .map_err(|e| format!("writing exonic truth failed: {e}"))?;

        self.intron
            .write_sparse(&intronic_path, gene_index)
            .map_err(|e| format!("writing intronic truth failed: {e}"))?;

        self.snp_ref
            .write_sparse(&ref_path, snp_index)
            .map_err(|e| format!("writing SNP ref truth failed: {e}"))?;

        self.snp_alt
            .write_sparse(&alt_path, snp_index)
            .map_err(|e| format!("writing SNP alt truth failed: {e}"))?;

        Ok(())
    }

    /// Compare two truth datasets for equality.
    ///
    /// This should be used in end-to-end tests against bam-quant output.
    ///
    /// NOTE:
    /// - Requires Scdata to implement meaningful equality.
    /// - If floating point issues arise, this should be relaxed to
    ///   approximate comparison.
    pub fn equals(&self, other: &Self) -> bool {
        self.gene == other.gene
            && self.intron == other.intron
            && self.snp_ref == other.snp_ref
            && self.snp_alt == other.snp_alt
    }
}
