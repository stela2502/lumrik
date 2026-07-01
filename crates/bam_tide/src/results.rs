// thrust_data.rs
use std::fs;

//use gtf_splice_index::{GeneId, RefBlock, SpliceIndex, Strand, TranscriptId};
//use snp_index::{Genome, SnpIndex, VcfReadOptions};

//use rand::rngs::SmallRng;
//use rand::{Rng, SeedableRng};

use std::path::Path;

use scdata::{FeatureIndex, MatrixValueType, Scdata};

use mapping_info::MappingInfo;

pub struct QuantData {
    pub gene: Scdata,
    pub intron: Scdata,
    pub snp_ref: Scdata,
    pub snp_alt: Scdata,
    pub report: MappingInfo,
}

impl Default for QuantData {
    fn default() -> Self {
        Self::new()
    }
}

impl QuantData {
    pub fn new() -> Self {
        Self {
            gene: Scdata::new(1, MatrixValueType::Real),
            intron: Scdata::new(1, MatrixValueType::Real),
            snp_ref: Scdata::new(1, MatrixValueType::Real),
            snp_alt: Scdata::new(1, MatrixValueType::Real),
            report: MappingInfo::new(None, 20.0, usize::MAX),
        }
    }

    /// Read truth matrices from a truth output directory.
    ///
    /// Expected layout:
    ///
    /// - `<base>/exonic`
    /// - `<base>/intronic`
    /// - `<base>/ref`
    /// - `<base>/alt`
    pub fn from_path<P, G, S>(base: P, gene_index: &G, snp_index: &S) -> Result<Self, String>
    where
        P: AsRef<std::path::Path>,
        G: FeatureIndex,
        S: FeatureIndex,
    {
        let base = base.as_ref();

        Ok(Self {
            gene: Scdata::read_matrix_market(base.join("exonic"), gene_index).map_err(|e| {
                format!(
                    "failed to read exonic truth from {:?}: {e}",
                    base.join("exonic")
                )
            })?,

            intron: Scdata::read_matrix_market(base.join("intronic"), gene_index).map_err(|e| {
                format!(
                    "failed to read intronic truth from {:?}: {e}",
                    base.join("intronic")
                )
            })?,

            snp_ref: Scdata::read_matrix_market(base.join("ref"), snp_index).map_err(|e| {
                format!(
                    "failed to read SNP ref truth from {:?}: {e}",
                    base.join("ref")
                )
            })?,

            snp_alt: Scdata::read_matrix_market(base.join("alt"), snp_index).map_err(|e| {
                format!(
                    "failed to read SNP alt truth from {:?}: {e}",
                    base.join("alt")
                )
            })?,

            report: MappingInfo::new(None, 0.0, usize::MAX),
        })
    }

    pub fn merge(&mut self, other: &Self) {
        self.gene.merge(&other.gene);
        self.intron.merge(&other.intron);
        self.snp_ref.merge(&other.snp_ref);
        self.snp_alt.merge(&other.snp_alt);
        self.report.merge(&other.report);
    }

    /// Finalize and write all truth matrices to disk.
    ///
    /// This mirrors the behavior of bam-quant output so that
    /// direct comparisons are possible.
    pub fn write<P: AsRef<Path>, T: FeatureIndex, F: FeatureIndex>(
        &mut self,
        base: P,
        min_cell_counts: usize,
        gene_index: &T,
        snp_index: Option<&F>,
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

        self.gene.finalize_for_cells(&cells, gene_index);
        self.intron.finalize_for_cells(&cells, gene_index);

        // --- write ---
        self.gene
            .write_sparse(&exonic_path, gene_index)
            .map_err(|e| format!("writing exonic truth failed: {e}"))?;

        self.intron
            .write_sparse(&intronic_path, gene_index)
            .map_err(|e| format!("writing intronic truth failed: {e}"))?;

        if let Some(snp_index) = snp_index {
            self.snp_alt.finalize_for_cells(&cells, snp_index);

            self.snp_ref
                .retain_features(&self.snp_alt.observed_feature_ids());

            self.snp_ref.finalize_for_cells(&cells, snp_index);

            // --- write ---

            self.snp_ref
                .write_sparse(&ref_path, snp_index)
                .map_err(|e| format!("writing SNP ref truth failed: {e}"))?;

            self.snp_alt
                .write_sparse(&alt_path, snp_index)
                .map_err(|e| format!("writing SNP alt truth failed: {e}"))?;
        }

        Ok(())
    }
}

impl QuantData {
    /// Compare two quantification bundles.
    ///
    /// Returns `Ok(())` if all matrices match, otherwise returns the first
    /// useful discrepancy message.
    pub fn compare(&self, other: &Self) -> Result<(), String> {
        self.gene.compare(&other.gene, "exonic")?;
        self.intron.compare(&other.intron, "intronic")?;
        self.snp_ref.compare(&other.snp_ref, "snp_ref")?;
        self.snp_alt.compare(&other.snp_alt, "snp_alt")?;
        Ok(())
    }

    /// Return true if all quantification matrices match.
    pub fn equals(&self, other: &Self) -> bool {
        self.compare(other).is_ok()
    }
}
