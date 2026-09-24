use anyhow::{Context, Result, bail};
use std::path::{Path, PathBuf};
use std::process::Command;

pub const BASE_URL: &str = "https://hgdownload.soe.ucsc.edu";

fn output(url: &str) -> Result<String> {
    let out = Command::new("wget")
        .args(["-qO-", url])
        .output()
        .with_context(|| format!("starting wget for {url}"))?;
    if !out.status.success() {
        bail!("resource unavailable: {url}");
    }
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

fn exists(url: &str) -> bool {
    Command::new("wget")
        .args(["--quiet", "--spider", url])
        .status()
        .is_ok_and(|s| s.success())
}

fn download(url: &str, path: &Path, required: bool) -> Result<bool> {
    if path.is_file() {
        return Ok(true);
    }
    if !exists(url) {
        if required {
            bail!("required UCSC resource unavailable: {url}");
        }
        return Ok(false);
    }
    eprintln!("[ommverse] fetch {url}");
    let status = Command::new("wget")
        .args(["--continue", "--output-document"])
        .arg(path)
        .arg(url)
        .status()
        .with_context(|| format!("starting wget for {url}"))?;
    if !status.success() {
        bail!("wget failed for {url}");
    }
    Ok(true)
}

fn hrefs(index: &str) -> impl Iterator<Item = &str> {
    index
        .split("href=\"")
        .skip(1)
        .filter_map(|s| s.split('"').next())
}

pub fn fetch(assembly: &str, cache_root: &Path) -> Result<PathBuf> {
    let root = cache_root.join(assembly);
    let genome_dir = root.join("Genome");
    let genes_dir = root.join("Genes");
    let protein_dir = root.join("Protein");
    std::fs::create_dir_all(&genome_dir)?;
    std::fs::create_dir_all(&genes_dir)?;
    std::fs::create_dir_all(&protein_dir)?;

    let bigzips = format!("{BASE_URL}/goldenPath/{assembly}/bigZips");
    let latest = format!("{bigzips}/latest");
    for (name, required) in [
        (format!("{assembly}.2bit"), true),
        (format!("{assembly}.chrom.sizes"), true),
        (format!("{assembly}.chromAlias.txt"), false),
    ] {
        let latest_url = format!("{latest}/{name}");
        let base_url = format!("{bigzips}/{name}");
        let url = if exists(&latest_url) {
            latest_url
        } else {
            base_url
        };
        download(&url, &genome_dir.join(&name), required)?;
    }

    let genes = format!("{bigzips}/genes");
    // UCSC does not publish knownGene for every assembly.  The historical
    // downloader therefore probes all supported GTF families and lets the
    // Ommverse builder consume whichever annotation is available.  Keep the
    // same behaviour here: the *set* of gene annotations is required, not any
    // one particular UCSC track.
    let mut have_gtf = false;
    for name in [
        format!("{assembly}.knownGene.gtf.gz"),
        format!("{assembly}.ncbiRefSeq.gtf.gz"),
        "refGene.gtf.gz".to_owned(),
    ] {
        have_gtf |= download(&format!("{genes}/{name}"), &genes_dir.join(&name), false)?;
    }
    if !have_gtf {
        bail!("no supported UCSC GTF annotation available for {assembly} under {genes}");
    }

    let uniprot_root = format!("{BASE_URL}/goldenPath/archive/{assembly}/uniprot");
    let index = output(&format!("{uniprot_root}/"))?;
    let mut releases: Vec<_> = hrefs(&index)
        .filter_map(|h| h.strip_suffix('/'))
        .filter(|h| {
            h.len() == 7
                && h.as_bytes()[4] == b'_'
                && h[..4].bytes().all(|b| b.is_ascii_digit())
                && h[5..].bytes().all(|b| b.is_ascii_digit())
        })
        .collect();
    releases.sort_unstable();
    let release = releases.last().context("no UCSC UniProt release found")?;
    let uniprot = format!("{uniprot_root}/{release}");
    for (name, required) in [
        ("version.txt", true),
        ("trackDb.txt", true),
        ("protMapInfo.tsv", true),
        ("liftInfo.json", false),
        ("unipFullSeq.bb", true),
        ("unipDomain.bb", false),
        ("unipLocTransMemb.bb", false),
        ("unipLocSignal.bb", false),
        ("unipLocCytopl.bb", false),
        ("unipLocExtra.bb", false),
        ("unipModif.bb", false),
        ("unipDisulfBond.bb", false),
        ("unipRepeat.bb", false),
        ("unipChain.bb", false),
        ("unipConflict.bb", false),
        ("unipInterest.bb", false),
        ("unipMut.bb", false),
        ("unipOther.bb", false),
        ("unipSplice.bb", false),
        ("unipStruct.bb", false),
        ("unipAliSwissprot.bb", false),
        ("unipAliTrembl.bb", false),
        ("unipToGenome.over.chain.gz", false),
        ("unipToGenomeLift.psl.gz", false),
    ] {
        download(
            &format!("{uniprot}/{name}"),
            &protein_dir.join(name),
            required,
        )?;
    }

    Ok(root)
}
