use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

pub const HG38_ENHANCERS_URL: &str =
    "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz";

pub fn enhancer_url(assembly: &str) -> Option<&'static str> {
    match assembly {
        "hg38" => Some(HG38_ENHANCERS_URL),
        _ => None,
    }
}

pub fn fetch(assembly: &str, root: &Path) -> Result<Option<PathBuf>> {
    let Some(url) = enhancer_url(assembly) else { return Ok(None); };
    let dir = root.join("Chromatin").join("FANTOM5");
    std::fs::create_dir_all(&dir)?;
    let out = dir.join("F5.hg38.enhancers.bed.gz");
    if !out.is_file() {
        let status = Command::new("wget")
            .args(["--continue", "--output-document"])
            .arg(&out).arg(url).status()
            .with_context(|| format!("starting wget for {url}"))?;
        if !status.success() { anyhow::bail!("wget failed for {url}"); }
    }
    Ok(Some(out))
}
