pub mod encode4;
pub mod fantom5;
pub mod ucsc;

use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;

const HG38: &str = include_str!("../../sources/hg38.yaml");

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SourceManifest {
    pub schema_version: u32,
    pub reference: ReferenceSpec,
    #[serde(default)]
    pub resources: Vec<ResourceSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceSpec {
    pub name: String,
    pub species: String,
    pub assembly: String,
    pub assembly_accession: Option<String>,
    /// NCBI taxonomy identifier used by organism-scoped external resources.
    #[serde(default)]
    pub taxonomy_id: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResourceSpec {
    pub kind: String,
    pub source: String,
    pub release: Option<String>,
    pub format: Option<String>,
    pub file: Option<PathBuf>,
    #[serde(default)]
    pub providers: Vec<ProviderSpec>,
    #[serde(default)]
    pub resolved: Option<ResolvedSource>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResolvedSource {
    pub provider: String,
    pub url: String,
    pub bytes: Option<u64>,
    pub sha256: Option<String>,
    pub retrieved_unix_seconds: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProviderSpec {
    pub name: String,
    pub url: String,
    #[serde(default)]
    pub priority: u32,
    #[serde(default)]
    pub adapter: Option<String>,
}

impl SourceManifest {
    pub fn from_yaml(text: &str) -> Result<Self> {
        serde_yaml::from_str(text).context("parsing Ommverse sources YAML")
    }

    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        serde_yaml::from_reader(file).with_context(|| format!("parsing {}", path.display()))
    }

    pub fn load_optional(path: impl AsRef<Path>) -> Result<Option<Self>> {
        let path = path.as_ref();
        if path.is_file() { Ok(Some(Self::load(path)?)) } else { Ok(None) }
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        serde_yaml::to_writer(file, self).with_context(|| format!("writing {}", path.display()))
    }
}

fn builtin(assembly: &str) -> Result<SourceManifest> {
    match assembly {
        "hg38" => SourceManifest::from_yaml(HG38),
        _ => Ok(SourceManifest {
            schema_version: 1,
            reference: ReferenceSpec {
                name: assembly.to_owned(),
                species: "unknown".to_owned(),
                assembly: assembly.to_owned(),
                assembly_accession: None,
                taxonomy_id: None,
            },
            resources: vec![ResourceSpec {
                kind: "reference".to_owned(),
                source: "UCSC".to_owned(),
                release: None,
                format: Some("ucsc-reference-tree".to_owned()),
                file: None,
                providers: vec![ProviderSpec {
                    name: "UCSC".to_owned(),
                    url: ucsc::BASE_URL.to_owned(),
                    priority: 10,
                    adapter: Some("ucsc".to_owned()),
                }],
                resolved: None,
            }],
        }),
    }
}

fn download(url: &str, path: &Path) -> Result<()> {
    if path.is_file() { return Ok(()); }
    if let Some(parent) = path.parent() { std::fs::create_dir_all(parent)?; }
    eprintln!("[ommverse] fetch {url}");
    let status = Command::new("wget")
        .args(["--continue", "--output-document"])
        .arg(path).arg(url).status()
        .with_context(|| format!("starting wget for {url}"))?;
    if !status.success() { bail!("wget failed for {url}"); }
    Ok(())
}

fn resolved(provider: String, url: String, path: Option<&Path>) -> ResolvedSource {
    let bytes = path.and_then(|p| std::fs::metadata(p).ok().map(|m| m.len()));
    let sha256 = path.and_then(|p| {
        let out = Command::new("sha256sum").arg(p).output().ok()?;
        if !out.status.success() { return None; }
        String::from_utf8_lossy(&out.stdout).split_whitespace().next().map(str::to_owned)
    });
    let retrieved_unix_seconds = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH).ok().map(|d| d.as_secs());
    ResolvedSource { provider, url, bytes, sha256, retrieved_unix_seconds }
}

#[derive(Debug, Clone, Default)]
pub struct AnnotationOptions {
    /// None means resolve the latest GENCODE release when this assembly is supported.
    pub gencode_release: Option<String>,
    /// An explicit user GTF disables automatic GENCODE fetching.
    pub explicit_gtf: Option<PathBuf>,
}

pub fn fetch_assembly_with_annotation(assembly: &str, cache_root: &Path, annotation: &AnnotationOptions) -> Result<PathBuf> {
    let mut manifest = builtin(assembly)?;
    if manifest.reference.name != assembly {
        bail!("source manifest names {} but {assembly} was requested", manifest.reference.name);
    }
    let root = cache_root.join(assembly);
    std::fs::create_dir_all(&root)?;

    let mut ran_ucsc = false;
    for resource in &mut manifest.resources {
        let mut providers = resource.providers.clone();
        providers.sort_by_key(|p| p.priority);
        let provider = providers.first().with_context(|| format!("{} has no provider", resource.kind))?;
        let provider_name = provider.name.clone();
        let provider_url = provider.url.clone();
        let provider_adapter = provider.adapter.clone();
        if provider_adapter.as_deref() == Some("ucsc") {
            if !ran_ucsc {
                ucsc::fetch_from_base_with_gencode(
                    assembly,
                    cache_root,
                    provider_url.trim_end_matches('/'),
                    annotation.gencode_release.as_deref(),
                    annotation.explicit_gtf.is_none(),
                )?;
                ran_ucsc = true;
            }
            resource.resolved = Some(resolved(provider_name, provider_url, None));
            continue;
        }
        let Some(rel) = &resource.file else { continue; };
        let path = root.join(rel);
        download(&provider_url, &path)?;
        resource.resolved = Some(resolved(provider_name, provider_url, Some(&path)));
    }
    if let Some(gtf) = &annotation.explicit_gtf {
        let canonical = gtf.canonicalize().with_context(|| format!("GTF {}", gtf.display()))?;
        manifest.resources.push(ResourceSpec {
            kind: "gene_annotation".to_owned(),
            source: "user".to_owned(),
            release: None,
            format: Some("gtf".to_owned()),
            file: Some(canonical.clone()),
            providers: Vec::new(),
            resolved: Some(resolved("user".to_owned(), canonical.display().to_string(), Some(&canonical))),
        });
    } else if let Some(series) = match assembly { "hg38" => Some("GENCODE human"), "mm39" => Some("GENCODE mouse"), _ => None } {
        manifest.resources.push(ResourceSpec {
            kind: "gene_annotation".to_owned(),
            source: series.to_owned(),
            release: Some(annotation.gencode_release.clone().unwrap_or_else(|| "latest".to_owned())),
            format: Some("gtf.gz".to_owned()),
            file: None,
            providers: Vec::new(),
            resolved: None,
        });
    }
    manifest.save(root.join("sources.yaml"))?;
    Ok(root)
}

pub fn fetch_assembly(assembly: &str, cache_root: &Path) -> Result<PathBuf> {
    fetch_assembly_with_annotation(assembly, cache_root, &AnnotationOptions::default())
}
