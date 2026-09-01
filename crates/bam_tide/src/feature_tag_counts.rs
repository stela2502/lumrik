use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{Context, Result, bail};
use fast_tag_mapper::{BuiltinTagSet, FastTagMapper};
use scdata::Scdata;

use crate::index::FastTagFeatureIndex;
use crate::ngs_normalizer::NgsNormalizerSupport;

use sc_beacon::{BackgroundConfig, CallConfig, FitConfig, run_from_scdata};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdditionalFeatureSource {
    BdSampleHuman,
    BdSampleMouse,
    Fasta(PathBuf),
}

impl FromStr for AdditionalFeatureSource {
    type Err = anyhow::Error;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "bd_sample_human" => Ok(Self::BdSampleHuman),
            "bd_sample_mouse" => Ok(Self::BdSampleMouse),
            value => {
                let path = PathBuf::from(value);
                let valid_fasta_name = value.ends_with(".fa")
                    || value.ends_with(".fasta")
                    || value.ends_with(".fa.gz")
                    || value.ends_with(".fasta.gz");

                if !valid_fasta_name {
                    bail!(
                        "additional feature source '{}' is neither a known built-in feature set nor a FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)",
                        value
                    );
                }

                if !path.is_file() {
                    bail!(
                        "additional feature FASTA '{}' does not exist or is not a file",
                        value
                    );
                }

                Ok(Self::Fasta(path))
            }
        }
    }
}

impl AdditionalFeatureSource {
    pub fn feature_type(&self) -> Result<String> {
        match self {
            Self::BdSampleHuman => Ok("bd_sample_human".to_string()),
            Self::BdSampleMouse => Ok("bd_sample_mouse".to_string()),
            Self::Fasta(path) => fasta_feature_type(path),
        }
    }
}

/// Additional-feature counts collected during normalization.
///
/// The normalizers own this object. It keeps the short-feature mapper together
/// with the Scdata carrying the corresponding feature/UMI counts, so the data
/// can be finalized later against the canonical cell set selected by GEX
/// quantification.
pub struct FeatureTagCounts {
    data: Scdata,
    mapper: FastTagMapper,
}

impl FeatureTagCounts {
    pub fn empty() -> Self {
        Self {
            data: NgsNormalizerSupport::new_feature_tag_table(),
            mapper: FastTagMapper::new(),
        }
    }

    pub fn from_sources(sources: &[AdditionalFeatureSource], min_hits: u32) -> Result<Self> {
        let mut mapper = FastTagMapper::new();

        for source in sources {
            match source {
                AdditionalFeatureSource::BdSampleHuman => {
                    mapper.add_builtin(BuiltinTagSet::Human);
                }
                AdditionalFeatureSource::BdSampleMouse => {
                    mapper.add_builtin(BuiltinTagSet::Mouse);
                }
                AdditionalFeatureSource::Fasta(path) => {
                    let feature_type = source.feature_type()?;
                    mapper.load_fasta_as(path, &feature_type).with_context(|| {
                        format!(
                            "failed to load additional feature FASTA '{}'",
                            path.display()
                        )
                    })?;
                }
            }
        }

        Ok(Self {
            data: NgsNormalizerSupport::new_feature_tag_table(),
            mapper: mapper.with_min_hits(min_hits),
        })
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn mapper(&self) -> &FastTagMapper {
        &self.mapper
    }

    pub(crate) fn data_mut(&mut self) -> &mut Scdata {
        &mut self.data
    }

    pub(crate) fn merge_table(&mut self, other: &Scdata) {
        self.data.merge(other);
    }

    /// Restrict feature counts to the canonical GEX cells, write one
    /// 10x-style directory per feature type, and run sc-beacon
    /// independently for every feature type.
    pub fn finalize_and_write(
        &mut self,
        cells: &HashSet<u64>,
        cell_barcode_len: usize,
        out: &Path,
    ) -> Result<()> {
        if self.data.is_empty() {
            return Ok(());
        }

        /*
         * Move the raw table out temporarily.
         *
         * selected  = canonical GEX cells
         * background = all remaining droplets/cells
         */
        let raw = std::mem::replace(
            &mut self.data,
            NgsNormalizerSupport::new_feature_tag_table(),
        );

        let (mut filtered, background) = raw.split_by_cells(cells);

        let background_config = BackgroundConfig::default();
        let fit_config = FitConfig::default();
        let call_config = CallConfig::default();

        let feature_index = FastTagFeatureIndex::new(&self.mapper);

        for (feature_type, feature_index) in feature_index.split_by_feature_type() {
            /*
             * ------------------------------------------------------------
             * Raw filtered feature counts
             * ------------------------------------------------------------
             */
            filtered.finalize_for_cells(cells, &feature_index);

            let out_dir = out.join(&feature_type);

            std::fs::create_dir_all(&out_dir)
                .with_context(|| format!("failed to create {}", out_dir.display()))?;

            filtered
                .write_sparse_with_cell_len(&out_dir, &feature_index, cell_barcode_len)
                .map_err(anyhow::Error::msg)
                .with_context(|| format!("writing feature table '{feature_type}'"))?;

            /*
             * ------------------------------------------------------------
             * sc-beacon
             * ------------------------------------------------------------
             *
             * Both Scdata objects may contain several feature classes.
             * The type-specific FeatureIndex restricts Beacon to this block.
             */
            match run_from_scdata(
                &filtered,
                &background,
                cells,
                cell_barcode_len,
                &feature_index,
                &background_config,
                &fit_config,
                &call_config,
            ) {
                Ok(mut beacon) => {
                    let beacon_out = out_dir.join("beacon");

                    if let Err(err) = beacon.write(&beacon_out, &feature_index, cell_barcode_len) {
                        write_beacon_error(&out_dir, &feature_type, "writing results", &err);
                    }
                }

                Err(err) => {
                    write_beacon_error(&out_dir, &feature_type, "running model", &err);
                }
            }
        }

        /*
         * We intentionally discard the background cells now.
         * From this point onward FeatureTagCounts contains only
         * canonical GEX cells.
         */
        self.data = filtered;

        Ok(())
    }

    /// Standalone normalizer output: retain every observed cell.
    pub fn finalize_and_write_all(&mut self, cell_barcode_len: usize, out: &Path) -> Result<()> {
        if self.data.is_empty() {
            return Ok(());
        }
        let cells = self.data.cell_ids();
        self.finalize_and_write(&cells, cell_barcode_len, out)
    }
}

fn write_beacon_error(out_dir: &Path, feature_type: &str, phase: &str, err: &anyhow::Error) {
    let beacon_out = out_dir.join("beacon");
    let error_path = beacon_out.join("sc_beacon.error.log");
    let message =
        format!("sc-beacon failed for feature type '{feature_type}' while {phase}\n\n{err:#}\n");

    eprintln!(
        "[bam-tide] sc-beacon failed for '{}' while {}: {:#}",
        feature_type, phase, err,
    );

    if let Err(write_err) =
        std::fs::create_dir_all(&beacon_out).and_then(|_| std::fs::write(&error_path, message))
    {
        eprintln!(
            "[bam-tide] could not write sc-beacon error log {}: {}",
            error_path.display(),
            write_err,
        );
    }
}

impl Default for FeatureTagCounts {
    fn default() -> Self {
        Self::empty()
    }
}

fn fasta_feature_type(path: &Path) -> Result<String> {
    let filename = path.file_name().and_then(|x| x.to_str()).with_context(|| {
        format!(
            "additional feature FASTA has no usable file name: {}",
            path.display()
        )
    })?;

    let name = filename
        .strip_suffix(".fa.gz")
        .or_else(|| filename.strip_suffix(".fasta.gz"))
        .or_else(|| filename.strip_suffix(".fa"))
        .or_else(|| filename.strip_suffix(".fasta"))
        .with_context(|| {
            format!(
                "additional feature file '{}' is not a supported FASTA file",
                path.display()
            )
        })?;

    if name.is_empty() {
        bail!(
            "could not derive feature type from FASTA '{}'",
            path.display()
        );
    }

    Ok(name.to_string())
}
