use crate::query::{Entity, Query};
use crate::Ommverse;
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::path::{Path, PathBuf};

pub const OMMVERSE_SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OmmverseSchema {
    pub format: String,
    pub version: u32,
    pub assembly: String,
    pub providers: BTreeMap<String, ProviderSchema>,
    pub searches: BTreeMap<String, SearchSchema>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProviderSchema {
    pub delivers: BTreeMap<String, DeliverableSchema>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub projections: Vec<AxisProjection>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DeliverableSchema {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub identities: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub axes: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AxisProjection {
    pub from: String,
    pub to: String,
    pub via: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SearchSchema {
    pub fields: BTreeSet<String>,
    pub required_any: Vec<BTreeSet<String>>,
}

impl OmmverseSchema {
    pub fn current(assembly: impl Into<String>) -> Self {
        let mut providers = BTreeMap::new();

        let mut gtf = BTreeMap::new();
        gtf.insert(
            "Gene".to_owned(),
            DeliverableSchema {
                identities: vec!["gene".to_owned()],
                axes: vec!["genomic".to_owned()],
            },
        );
        gtf.insert(
            "Transcript".to_owned(),
            DeliverableSchema {
                identities: vec!["transcript".to_owned()],
                axes: vec!["genomic".to_owned(), "transcriptomic".to_owned()],
            },
        );
        gtf.insert(
            "Protein".to_owned(),
            DeliverableSchema {
                identities: vec!["protein".to_owned()],
                axes: vec!["proteomic".to_owned()],
            },
        );
        providers.insert(
            "gtf_splice_index".to_owned(),
            ProviderSchema {
                delivers: gtf,
                projections: vec![
                    AxisProjection { from: "genomic".into(), to: "transcriptomic".into(), via: "Transcript".into() },
                    AxisProjection { from: "transcriptomic".into(), to: "genomic".into(), via: "Transcript".into() },
                    AxisProjection { from: "proteomic".into(), to: "genomic".into(), via: "Transcript".into() },
                ],
            },
        );

        let mut interpro = BTreeMap::new();
        interpro.insert(
            "ProteinDomain".to_owned(),
            DeliverableSchema {
                identities: vec!["protein".to_owned()],
                axes: vec!["proteomic".to_owned()],
            },
        );
        providers.insert(
            "interpro".to_owned(),
            ProviderSchema { delivers: interpro, projections: Vec::new() },
        );

        let mut searches = BTreeMap::new();
        searches.insert(
            "gene".to_owned(),
            SearchSchema {
                fields: ["name", "id", "symbol", "chromosome", "chrom", "chr", "position", "pos"].into_iter().map(str::to_owned).collect(),
                required_any: Vec::new(),
            },
        );
        searches.insert(
            "transcript".to_owned(),
            SearchSchema {
                fields: ["name", "id"].into_iter().map(str::to_owned).collect(),
                required_any: vec![["name", "id"].into_iter().map(str::to_owned).collect()],
            },
        );
        searches.insert(
            "protein".to_owned(),
            SearchSchema {
                fields: ["name", "id", "accession"].into_iter().map(str::to_owned).collect(),
                required_any: vec![["name", "id", "accession"].into_iter().map(str::to_owned).collect()],
            },
        );
        searches.insert(
            "variant".to_owned(),
            SearchSchema {
                fields: ["chromosome", "chrom", "chr", "position", "pos", "clinical_effect"].into_iter().map(str::to_owned).collect(),
                required_any: vec![
                    ["chromosome", "chrom", "chr"].into_iter().map(str::to_owned).collect(),
                    ["position", "pos"].into_iter().map(str::to_owned).collect(),
                ],
            },
        );

        Self {
            format: "ommverse-schema".to_owned(),
            version: OMMVERSE_SCHEMA_VERSION,
            assembly: assembly.into(),
            providers,
            searches,
        }
    }

    /// Build the capability schema for one already-loaded index.
    ///
    /// This is primarily the bridge for pre-schema Ommverse files: they are
    /// deserialized once, their optional payloads are inspected, and a tiny
    /// sidecar is written so subsequent query validation does not need to open
    /// the large index.
    pub fn for_ommverse(omm: &Ommverse) -> Self {
        let mut schema = Self::current(omm.assembly.clone());
        if omm.interpro_entries.is_empty() {
            schema.providers.remove("interpro");
        }

        let genomic = |_name: &str| DeliverableSchema {
            identities: Vec::new(),
            axes: vec!["genomic".to_owned()],
        };
        let provider = |name: &str, deliverable: DeliverableSchema| ProviderSchema {
            delivers: [(name.to_owned(), deliverable)].into_iter().collect(),
            projections: Vec::new(),
        };

        if !omm.chromatin.is_empty() {
            schema.providers.insert("fantom5".to_owned(), provider("ChromatinElement", genomic("ChromatinElement")));
        }
        if !omm.protein_binding.regions.is_empty() {
            schema.providers.insert("encode4_protein_binding".to_owned(), provider("ProteinBindingRegion", genomic("ProteinBindingRegion")));
        }
        if !omm.ctcf.anchors.is_empty() || !omm.ctcf.domains.is_empty() {
            let mut delivers = BTreeMap::new();
            delivers.insert("CtcfAnchor".to_owned(), genomic("CtcfAnchor"));
            delivers.insert("CtcfDomain".to_owned(), genomic("CtcfDomain"));
            schema.providers.insert("encode4_ctcf".to_owned(), ProviderSchema { delivers, projections: Vec::new() });
        }
        if !omm.experimental_loops.anchors.is_empty() || !omm.experimental_loops.loops.is_empty() {
            let mut delivers = BTreeMap::new();
            delivers.insert("ExperimentalLoopAnchor".to_owned(), genomic("ExperimentalLoopAnchor"));
            delivers.insert("ExperimentalLoop".to_owned(), DeliverableSchema { identities: Vec::new(), axes: Vec::new() });
            schema.providers.insert("experimental_loops".to_owned(), ProviderSchema { delivers, projections: Vec::new() });
        }

        // External resources remain queryable without being copied into the large
        // serialized core.  Advertise only resources recorded by this index's
        // resolved source manifest.
        if let Some(manifest) = &omm.sources {
            for resource in &manifest.resources {
                match resource.kind.as_str() {
                    "variants" => {
                        schema.providers.insert("variants".to_owned(), provider("Variant", DeliverableSchema { identities: vec!["variant".into(), "rsid".into()], axes: vec!["genomic".into()] }));
                    }
                    "clinical_effects" => {
                        schema.providers.insert("clinvar".to_owned(), provider("ClinicalEffect", DeliverableSchema { identities: vec!["variant".into()], axes: Vec::new() }));
                    }
                    "protein_interactions" | "protein_interaction_nodes" | "protein_interaction_aliases" | "interactions" | "string" => {
                        schema.providers.insert("string".to_owned(), provider("ProteinInteraction", DeliverableSchema { identities: vec!["protein".into()], axes: Vec::new() }));
                    }
                    _ => {}
                }
            }
        }
        schema
    }

    pub fn validate_query(&self, query: &Query) -> Result<()> {
        let entity = match query.entity {
            Entity::Gene => "gene",
            Entity::Transcript => "transcript",
            Entity::Protein => "protein",
            Entity::Variant => "variant",
        };
        let search = self.searches.get(entity)
            .with_context(|| format!("this Ommverse schema cannot search {entity}"))?;
        let fields: BTreeSet<_> = query.predicates.iter().map(|p| p.field()).collect();
        for field in &fields {
            if !search.fields.contains(*field) {
                bail!("{entity} field '{field}' is not searchable in this Ommverse schema");
            }
        }
        for alternatives in &search.required_any {
            if !alternatives.iter().any(|field| fields.contains(field.as_str())) {
                bail!("{entity} query requires one of: {}", alternatives.iter().cloned().collect::<Vec<_>>().join(", "));
            }
        }
        Ok(())
    }

    pub fn save_for_index(&self, index: impl AsRef<Path>) -> Result<PathBuf> {
        let path = schema_path(index);
        let file = File::create(&path).with_context(|| format!("creating schema {}", path.display()))?;
        serde_yaml::to_writer(file, self)?;
        Ok(path)
    }

    pub fn load_for_index(index: impl AsRef<Path>) -> Result<Self> {
        let path = schema_path(index);
        let file = File::open(&path).with_context(|| format!("opening Ommverse schema {}; rebuild the index to generate it", path.display()))?;
        let schema: Self = serde_yaml::from_reader(file)?;
        if schema.format != "ommverse-schema" { bail!("{} is not an Ommverse schema", path.display()); }
        if schema.version != OMMVERSE_SCHEMA_VERSION {
            bail!("unsupported Ommverse schema version {} (expected {})", schema.version, OMMVERSE_SCHEMA_VERSION);
        }
        Ok(schema)
    }
}

pub fn schema_path(index: impl AsRef<Path>) -> PathBuf {
    let index = index.as_ref();
    let mut name = index.file_name().unwrap_or_default().to_os_string();
    name.push(".schema.yaml");
    index.with_file_name(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_describes_protein_domain_projection_to_genome() {
        let schema = OmmverseSchema::current("hg38");
        let domain = &schema.providers["interpro"].delivers["ProteinDomain"];
        assert_eq!(domain.axes, vec!["proteomic"]);
        assert!(schema.providers["gtf_splice_index"].projections.iter().any(|p| p.from == "proteomic" && p.to == "genomic"));
    }

    #[test]
    fn schema_rejects_unknown_query_field_before_index_load() {
        let schema = OmmverseSchema::current("hg38");
        let query = Query::parse("SELECT gene WHERE chromosome = \"chr2\" AND position = 1:2 AND nonsense = 7").unwrap();
        assert!(schema.validate_query(&query).unwrap_err().to_string().contains("nonsense"));
    }
}
