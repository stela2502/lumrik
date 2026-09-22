use anyhow::{Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

/// ENCODE4 representative DNA-associated-protein peaks.  This is the UCSC
/// integrated rPeak backbone: 912 factors aggregated across 1,152 biosamples.
pub const HG38_TF_RPEAKS_URL: &str =
    "https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode4/regulation/tfRpeak/TFrPeakClusters.bb";

pub fn tf_rpeaks_url(assembly: &str) -> Option<&'static str> {
    match assembly {
        "hg38" => Some(HG38_TF_RPEAKS_URL),
        _ => None,
    }
}

pub fn fetch(assembly: &str, root: &Path) -> Result<Option<PathBuf>> {
    let Some(url) = tf_rpeaks_url(assembly) else { return Ok(None); };
    let dir = root.join("Chromatin").join("ENCODE4");
    std::fs::create_dir_all(&dir)?;
    let out = dir.join("TFrPeakClusters.bb");
    if !out.is_file() {
        eprintln!("[ommverse] fetch {url}");
        let status = Command::new("wget")
            .args(["--continue", "--output-document"])
            .arg(&out).arg(url).status()
            .with_context(|| format!("starting wget for {url}"))?;
        if !status.success() { anyhow::bail!("wget failed for {url}"); }
    }
    Ok(Some(out))
}
