use std::path::{Path, PathBuf};
use std::str::FromStr;

use anyhow::{bail, Context, Result};
use fast_tag_mapper::{BuiltinTagSet, FastTagMapper};

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

                let valid_fasta_name =
                    value.ends_with(".fa")
                        || value.ends_with(".fasta")
                        || value.ends_with(".fa.gz")
                        || value.ends_with(".fasta.gz");

                if !valid_fasta_name {
                    bail!(
                        "additional feature source '{}' is neither a known built-in \
                         feature set nor a FASTA file (.fa, .fasta, .fa.gz, .fasta.gz)",
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

pub fn build_mapper(
    sources: &[AdditionalFeatureSource],
    min_hits: u32,
) -> Result<Option<FastTagMapper>> {
    if sources.is_empty() {
        return Ok(None);
    }

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

                mapper
                    .load_fasta_as(path, &feature_type)
                    .with_context(|| {
                        format!(
                            "failed to load additional feature FASTA '{}'",
                            path.display()
                        )
                    })?;
            }
        }
    }

    Ok(Some(mapper.with_min_hits(min_hits)))
}

fn fasta_feature_type(path: &Path) -> Result<String> {
    let filename = path
        .file_name()
        .and_then(|x| x.to_str())
        .with_context(|| {
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