// thrust_data.rs
use std::fs;

//use gtf_splice_index::{GeneId, RefBlock, SpliceIndex, Strand, TranscriptId};
//use snp_index::{Genome, SnpIndex, VcfReadOptions};

//use rand::rngs::SmallRng;
//use rand::{Rng, SeedableRng};

use std::collections::HashSet;
use std::path::Path;

use sc_beacon::{KneeCountFit, fit_knee_counts, write_count_qc};

use scdata::{FeatureIndex, MatrixValueType, Scdata};

use mapping_info::MappingInfo;

pub struct QuantData {
    pub gene: Scdata,
    pub intron: Scdata,
    pub snp_ref: Scdata,
    pub snp_alt: Scdata,
    pub report: MappingInfo,
}

#[derive(Debug, Clone)]
pub struct CellAccounting {
    pub exonic_cells: usize,
    pub intronic_cells: usize,
    pub exonic_or_intronic_cells: usize,
    pub exonic_umis_before_filter: usize,
    pub exonic_umis_after_filter: usize,
    pub exonic_thresholds: Vec<(usize, usize)>,
}

impl std::fmt::Display for CellAccounting {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Cell accounting")?;
        writeln!(f, "---------------")?;
        writeln!(f, "exonic cells before cutoff: {}", self.exonic_cells)?;
        writeln!(f, "intronic cells before cutoff: {}", self.intronic_cells)?;
        writeln!(
            f,
            "cells with exonic or intronic evidence: {}",
            self.exonic_or_intronic_cells
        )?;
        writeln!(
            f,
            "unique exonic UMIs before filter: {}",
            self.exonic_umis_before_filter
        )?;
        writeln!(
            f,
            "unique exonic UMIs after filter: {}",
            self.exonic_umis_after_filter
        )?;
        for (threshold, count) in &self.exonic_thresholds {
            writeln!(f, "exonic cells with >= {threshold} UMIs: {count}")?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct BeaconCellCalling {
    pub fit: KneeCountFit,
    pub cells: Vec<(u64, u32, bool)>,
    pub retained: HashSet<u64>,
}

impl BeaconCellCalling {
    pub fn write_qc<P: AsRef<Path>>(&self, out_dir: P) -> Result<(), String> {
        let counts: Vec<u32> = self.cells.iter().map(|(_, count, _)| *count).collect();
        write_count_qc(&counts, &self.fit, out_dir).map_err(|e| e.to_string())
    }

    pub fn write_tsv<P: AsRef<Path>>(&self, path: P) -> Result<(), String> {
        use std::io::Write;

        let path = path.as_ref();
        let mut rows = self.cells.clone();
        rows.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let mut out =
            std::fs::File::create(path).map_err(|e| format!("creating {}: {e}", path.display()))?;
        writeln!(out, "rank\tcell_id\tumi_count\tcalled")
            .map_err(|e| format!("writing {}: {e}", path.display()))?;
        for (rank, (cell_id, count, called)) in rows.into_iter().enumerate() {
            writeln!(out, "{}\t{}\t{}\t{}", rank + 1, cell_id, count, called)
                .map_err(|e| format!("writing {}: {e}", path.display()))?;
        }
        Ok(())
    }
}

impl Default for QuantData {
    fn default() -> Self {
        Self::new()
    }
}

impl QuantData {
    /// Non-destructive cell accounting before the canonical GEX cutoff.
    pub fn cell_accounting(&self) -> CellAccounting {
        const THRESHOLDS: [usize; 7] = [1, 10, 50, 100, 200, 400, 1_000];
        let exonic = self.gene.cell_ids();
        let intronic = self.intron.cell_ids();
        let mut union = exonic.clone();
        union.extend(intronic.iter().copied());
        CellAccounting {
            exonic_cells: exonic.len(),
            intronic_cells: intronic.len(),
            exonic_or_intronic_cells: union.len(),
            exonic_umis_before_filter: self.gene.total_umis(),
            exonic_umis_after_filter: self.gene.total_umis(),
            exonic_thresholds: THRESHOLDS
                .into_iter()
                .map(|threshold| (threshold, self.gene.cells_with_min_umis(threshold)))
                .collect(),
        }
    }

    /// Return exonic cell ids with at least `min_umis` unique UMIs.
    ///
    /// This is the explicit fixed-cutoff counterpart to sc-beacon cell calling.
    pub fn cells_with_min_exonic_umis(&self, min_umis: usize) -> HashSet<u64> {
        self.gene
            .cell_umi_counts()
            .into_iter()
            .filter_map(|(cell_id, umi_count)| {
                ((umi_count as usize) >= min_umis).then_some(cell_id)
            })
            .collect()
    }

    /// Call cells from the complete exonic UMI-count distribution using
    /// sc-beacon's barcode-rank knee detector.
    pub fn beacon_cell_calling(&self) -> Result<BeaconCellCalling, String> {
        let cell_counts = self.gene.cell_umi_counts();
        let counts: Vec<u32> = cell_counts.iter().map(|(_, count)| *count).collect();
        let fit = fit_knee_counts(&counts).map_err(|e| e.to_string())?;
        let cells: Vec<(u64, u32, bool)> = cell_counts
            .into_iter()
            .map(|(cell_id, count)| (cell_id, count, count >= fit.umi_cutoff))
            .collect();
        let retained = cells
            .iter()
            .filter_map(|(cell_id, _, called)| called.then_some(*cell_id))
            .collect();
        Ok(BeaconCellCalling {
            fit,
            cells,
            retained,
        })
    }

    /// Write an evidence-preserving `unfiltered/` snapshot, then apply the
    /// canonical exonic cell cutoff and write the normal filtered matrices.
    ///
    /// `cell_barcode_len` is optional so the standalone bam-quant behavior
    /// remains unchanged (32-base rendering), while Nelrune can keep writing
    /// the chemistry-specific barcode length.
    pub fn write_with_unfiltered<P, T, F>(
        &mut self,
        base: P,
        min_umi_count: usize,
        gene_index: &T,
        snp_index: Option<&F>,
        cell_barcode_len: Option<usize>,
    ) -> Result<(std::collections::HashSet<u64>, CellAccounting), String>
    where
        P: AsRef<Path>,
        T: FeatureIndex,
        F: FeatureIndex,
    {
        let base = base.as_ref();
        let mut accounting = self.cell_accounting();

        // Prepare each signal independently so the evidence snapshot does not
        // force intronic/SNP cells onto the exonic whitelist.
        self.gene.finalize_for_export(1, gene_index);
        self.intron.finalize_for_export(1, gene_index);
        if let Some(snp_index) = snp_index {
            self.snp_ref.finalize_for_export(1, snp_index);
            self.snp_alt.finalize_for_export(1, snp_index);
        }
        self.write_finalized_impl(
            &base.join("unfiltered"),
            gene_index,
            snp_index,
            cell_barcode_len,
        )?;

        let retained = self.finalize_for_export(min_umi_count, gene_index, snp_index);
        accounting.exonic_umis_after_filter = self.gene.total_umis();
        self.write_finalized_impl(base, gene_index, snp_index, cell_barcode_len)?;

        Ok((retained, accounting))
    }

    /// Write an unfiltered snapshot and a filtered export using an explicit
    /// caller-provided cell set. This lets sc-beacon own cell selection while
    /// QuantData remains responsible for consistent matrix filtering/export.
    pub fn write_with_unfiltered_for_cells<P, T, F>(
        &mut self,
        base: P,
        keep: &HashSet<u64>,
        gene_index: &T,
        snp_index: Option<&F>,
        cell_barcode_len: Option<usize>,
    ) -> Result<CellAccounting, String>
    where
        P: AsRef<Path>,
        T: FeatureIndex,
        F: FeatureIndex,
    {
        let base = base.as_ref();
        let mut accounting = self.cell_accounting();

        self.gene.finalize_for_export(1, gene_index);
        self.intron.finalize_for_export(1, gene_index);
        if let Some(snp_index) = snp_index {
            self.snp_ref.finalize_for_export(1, snp_index);
            self.snp_alt.finalize_for_export(1, snp_index);
        }
        self.write_finalized_impl(
            &base.join("unfiltered"),
            gene_index,
            snp_index,
            cell_barcode_len,
        )?;

        self.gene.finalize_for_cells(keep, gene_index);
        accounting.exonic_umis_after_filter = self.gene.total_umis();
        self.intron.finalize_for_cells(keep, gene_index);
        if let Some(snp_index) = snp_index {
            self.snp_alt.finalize_for_cells(keep, snp_index);
            self.snp_ref
                .retain_features(&self.snp_alt.observed_feature_ids());
            self.snp_ref.finalize_for_cells(keep, snp_index);
        }
        self.write_finalized_impl(base, gene_index, snp_index, cell_barcode_len)?;

        Ok(accounting)
    }

    pub fn finalize_for_export<T: FeatureIndex, F: FeatureIndex>(
        &mut self,
        min_umi_count: usize,
        gene_index: &T,
        snp_index: Option<&F>,
    ) -> std::collections::HashSet<u64> {
        self.gene.finalize_for_export(min_umi_count, gene_index);

        let cells: std::collections::HashSet<u64> =
            self.gene.export_cell_ids().iter().copied().collect();

        self.gene.finalize_for_cells(&cells, gene_index);

        self.intron.finalize_for_cells(&cells, gene_index);

        if let Some(snp_index) = snp_index {
            self.snp_alt.finalize_for_cells(&cells, snp_index);

            self.snp_ref
                .retain_features(&self.snp_alt.observed_feature_ids());

            self.snp_ref.finalize_for_cells(&cells, snp_index);
        }

        cells
    }

    pub fn write_finalized<P, T, F>(
        &mut self,
        base: P,
        gene_index: &T,
        snp_index: Option<&F>,
        cell_barcode_len: usize,
    ) -> Result<(), String>
    where
        P: AsRef<Path>,
        T: FeatureIndex,
        F: FeatureIndex,
    {
        self.write_finalized_impl(base.as_ref(), gene_index, snp_index, Some(cell_barcode_len))
    }

    fn write_finalized_impl<T: FeatureIndex, F: FeatureIndex>(
        &mut self,
        base: &Path,
        gene_index: &T,
        snp_index: Option<&F>,
        cell_barcode_len: Option<usize>,
    ) -> Result<(), String> {
        let exonic_path = base.join("exonic");
        let intronic_path = base.join("intronic");
        let ref_path = base.join("ref");
        let alt_path = base.join("alt");

        fs::create_dir_all(&exonic_path)
            .map_err(|e| format!("failed to create {:?}: {e}", exonic_path))?;
        fs::create_dir_all(&intronic_path)
            .map_err(|e| format!("failed to create {:?}: {e}", intronic_path))?;

        match cell_barcode_len {
            Some(len) => {
                self.gene
                    .write_sparse_with_cell_len(&exonic_path, gene_index, len)
                    .map_err(|e| format!("writing exonic truth failed: {e}"))?;
                self.intron
                    .write_sparse_with_cell_len(&intronic_path, gene_index, len)
                    .map_err(|e| format!("writing intronic truth failed: {e}"))?;
            }
            None => {
                self.gene
                    .write_sparse(&exonic_path, gene_index)
                    .map_err(|e| format!("writing exonic truth failed: {e}"))?;
                self.intron
                    .write_sparse(&intronic_path, gene_index)
                    .map_err(|e| format!("writing intronic truth failed: {e}"))?;
            }
        }

        if let Some(snp_index) = snp_index {
            fs::create_dir_all(&ref_path)
                .map_err(|e| format!("failed to create {:?}: {e}", ref_path))?;
            fs::create_dir_all(&alt_path)
                .map_err(|e| format!("failed to create {:?}: {e}", alt_path))?;
            match cell_barcode_len {
                Some(len) => {
                    self.snp_ref
                        .write_sparse_with_cell_len(&ref_path, snp_index, len)
                        .map_err(|e| format!("writing SNP ref truth failed: {e}"))?;
                    self.snp_alt
                        .write_sparse_with_cell_len(&alt_path, snp_index, len)
                        .map_err(|e| format!("writing SNP alt truth failed: {e}"))?;
                }
                None => {
                    self.snp_ref
                        .write_sparse(&ref_path, snp_index)
                        .map_err(|e| format!("writing SNP ref truth failed: {e}"))?;
                    self.snp_alt
                        .write_sparse(&alt_path, snp_index)
                        .map_err(|e| format!("writing SNP alt truth failed: {e}"))?;
                }
            }
        }
        Ok(())
    }

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
        let gene_merge = self.gene.merge(&other.gene);
        let intron_merge = self.intron.merge(&other.intron);
        let snp_ref_merge = self.snp_ref.merge(&other.snp_ref);
        let snp_alt_merge = self.snp_alt.merge(&other.snp_alt);

        self.report.merge(&other.report);
        self.report.merge(&gene_merge);
        self.report.merge(&intron_merge);
        self.report.merge(&snp_ref_merge);
        self.report.merge(&snp_alt_merge);
    }

    /// Finalize and write all truth matrices to disk.
    ///
    /// Kept for callers that explicitly want only the filtered output.
    pub fn write<P: AsRef<Path>, T: FeatureIndex, F: FeatureIndex>(
        &mut self,
        base: P,
        min_umi_count: usize,
        gene_index: &T,
        snp_index: Option<&F>,
    ) -> Result<(), String> {
        self.finalize_for_export(min_umi_count, gene_index, snp_index);
        self.write_finalized_impl(base.as_ref(), gene_index, snp_index, None)
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
