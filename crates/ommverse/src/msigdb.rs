//! MSigDB GMT gene-set sidecar. Gene-set membership is association, not direct TF binding.
//! The index is independent of the OMM bincode format, so existing genome indexes stay valid.
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeneSet {
    pub name: String,
    pub description: String,
    pub source: String,
    pub release: String,
    pub genes: Vec<String>,
}

#[derive(Serialize, Deserialize)]
struct Index {
    version: u32,
    sets: Vec<GeneSet>,
    by_gene: BTreeMap<String, Vec<usize>>,
}

pub fn sidecar(index: &Path) -> PathBuf {
    PathBuf::from(format!("{}.msigdb.bin", index.display()))
}

fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}


/// Fetch the official human MSigDB GMT collections using curl. Downloads are
/// cached and atomically renamed only after validating the GMT structure.
/// The importer never redistributes MSigDB data.
pub fn download_collections(cache: &Path, release: &str) -> Result<Vec<PathBuf>> {
    use std::process::Command;
    if release.is_empty() || !release.chars().all(|c| c.is_ascii_alphanumeric() || c == '.') {
        bail!("invalid MSigDB release: {release}");
    }
    fs::create_dir_all(cache)?;
    let base = format!("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{release}");
    let mut paths = Vec::new();
    for collection in ["h", "c2", "c3"] {
        let name = format!("{collection}.all.v{release}.symbols.gmt");
        let target = cache.join(&name);
        if target.exists() {
            validate_gmt(&target)?;
        } else {
            let partial = cache.join(format!("{name}.part"));
            let url = format!("{base}/{name}");
            eprintln!("[ommverse] downloading {url}");
            let status = Command::new("curl")
                .args(["--fail", "--location", "--retry", "3", "--connect-timeout", "20", "--max-time", "600", "--silent", "--show-error", "--output"])
                .arg(&partial).arg(&url).status()
                .context("running curl (install curl to enable automatic MSigDB downloads)")?;
            if !status.success() {
                let _ = fs::remove_file(&partial);
                bail!("MSigDB download failed for {url} (curl exit: {status})");
            }
            if let Err(error) = validate_gmt(&partial) {
                let _ = fs::remove_file(&partial);
                return Err(error).with_context(|| format!("invalid download from {url}"));
            }
            fs::rename(&partial, &target)?;
        }
        paths.push(target);
    }
    Ok(paths)
}

fn validate_gmt(path: &Path) -> Result<()> {
    let mut lines = reader(path)?.lines();
    let first = lines.next().ok_or_else(|| anyhow::anyhow!("empty GMT: {}", path.display()))??;
    let fields: Vec<_> = first.split('\t').collect();
    if fields.len() < 3 || fields[0].trim().is_empty() || fields[2..].iter().all(|g| g.trim().is_empty()) {
        bail!("invalid GMT header/first row in {}", path.display());
    }
    Ok(())
}

pub fn ingest(index: &Path, paths: &[PathBuf], release: &str) -> Result<usize> {
    if !index.is_file() { bail!("Ommverse index not found: {}", index.display()); }
    if release.trim().is_empty() { bail!("--release must not be empty"); }
    let mut sets = Vec::new();
    let mut seen = BTreeSet::new();
    for path in paths {
        let source = path.file_name().unwrap_or_default().to_string_lossy().to_string();
        for (line_no, line) in reader(path)?.lines().enumerate() {
            let line = line.with_context(|| format!("{}:{}", path.display(), line_no + 1))?;
            if line.trim().is_empty() { continue; }
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() < 3 { bail!("{}:{}: expected GMT name, description, genes", path.display(), line_no + 1); }
            let name = fields[0].trim();
            if name.is_empty() { bail!("{}:{}: empty gene-set name", path.display(), line_no + 1); }
            if !seen.insert(name.to_ascii_uppercase()) {
                bail!("duplicate gene-set name {name} across input files; use non-overlapping collections");
            }
            let genes: Vec<String> = fields[2..].iter().map(|g| g.trim().to_string())
                .filter(|g| !g.is_empty()).collect::<BTreeSet<_>>().into_iter().collect();
            if genes.is_empty() { bail!("{}:{}: empty gene set {name}", path.display(), line_no + 1); }
            sets.push(GeneSet { name: name.to_string(), description: fields[1].to_string(), source: source.clone(), release: release.to_string(), genes });
        }
    }
    let mut by_gene: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, set) in sets.iter().enumerate() {
        for gene in &set.genes { by_gene.entry(gene.to_ascii_uppercase()).or_default().push(i); }
    }
    let target = sidecar(index);
    let tmp = target.with_extension("msigdb.tmp");
    let mut out = BufWriter::new(File::create(&tmp)?);
    bincode::serialize_into(&mut out, &Index { version: 1, sets, by_gene })?;
    out.flush()?;
    drop(out);
    fs::rename(&tmp, &target)?;
    Ok(load(index)?.sets.len())
}

fn load(index: &Path) -> Result<Index> {
    let path = sidecar(index);
    let data: Index = bincode::deserialize_from(BufReader::new(File::open(&path)
        .with_context(|| format!("MSigDB sidecar missing: {} (run ingest-msigdb)", path.display()))?))?;
    if data.version != 1 { bail!("unsupported MSigDB sidecar version {}", data.version); }
    Ok(data)
}

pub fn find_gene(index: &Path, gene: &str) -> Result<Vec<GeneSet>> {
    let data = load(index)?;
    Ok(data.by_gene.get(&gene.to_ascii_uppercase()).into_iter().flatten()
        .map(|&i| data.sets[i].clone()).collect())
}

pub fn find_set(index: &Path, name: &str) -> Result<Vec<GeneSet>> {
    Ok(load(index)?.sets.into_iter().filter(|s| s.name.eq_ignore_ascii_case(name)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gmt_import_and_reverse_lookup() {
        let dir = std::env::temp_dir().join(format!("omm-msig-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let index = dir.join("test.ommverse");
        fs::write(&index, "fixture").unwrap();
        let gmt = dir.join("c3.gmt");
        fs::write(&gmt, "TP53_SET\tPMID:1\tTP53\tCDKN1A\nOTHER\tNA\tMDM2\tTP53\n").unwrap();
        assert_eq!(ingest(&index, &[gmt], "test-release").unwrap(), 2);
        assert_eq!(find_gene(&index, "tp53").unwrap().len(), 2);
        assert_eq!(find_set(&index, "tp53_set").unwrap()[0].genes.len(), 2);
        fs::remove_dir_all(dir).unwrap();
    }
}

/// Read all indexed sets for regulator-to-target matching.
pub fn all_sets(index: &Path) -> Result<Vec<GeneSet>> { Ok(load(index)?.sets) }
