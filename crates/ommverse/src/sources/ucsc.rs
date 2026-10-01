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


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GencodeAnnotation {
    pub release: String,
    pub path: PathBuf,
    pub url: String,
}

fn gencode_series(assembly: &str) -> Option<(&'static str, &'static str)> {
    match assembly {
        "hg38" => Some(("Gencode_human", "")),
        "mm39" => Some(("Gencode_mouse", "M")),
        _ => None,
    }
}

fn latest_gencode_release(root: &str, prefix: &str) -> Result<String> {
    let index = output(root)?;
    let mut releases = hrefs(&index)
        .filter_map(|href| href.strip_prefix("release_").and_then(|x| x.strip_suffix('/')))
        .filter(|release| {
            let numeric = release.strip_prefix(prefix).unwrap_or("");
            release.starts_with(prefix) && !numeric.is_empty() && numeric.bytes().all(|b| b.is_ascii_digit())
        })
        .map(str::to_owned)
        .collect::<Vec<_>>();
    releases.sort_by_key(|release| {
        release.strip_prefix(prefix).and_then(|x| x.parse::<u32>().ok()).unwrap_or(0)
    });
    releases.pop().with_context(|| format!("no GENCODE releases found under {root}"))
}

pub fn fetch_gencode_annotation(
    assembly: &str,
    genes_dir: &Path,
    requested_release: Option<&str>,
) -> Result<Option<GencodeAnnotation>> {
    let Some((series, prefix)) = gencode_series(assembly) else {
        if requested_release.is_some() {
            bail!("GENCODE annotation is not configured for assembly {assembly}; use --gtf with an explicit annotation");
        }
        return Ok(None);
    };
    let root = format!("https://ftp.ebi.ac.uk/pub/databases/gencode/{series}");
    let release = match requested_release {
        None | Some("latest") => latest_gencode_release(&root, prefix)?,
        Some(value) => {
            let value = value.trim();
            if value.starts_with(prefix) { value.to_owned() } else { format!("{prefix}{value}") }
        }
    };
    let version = release.strip_prefix(prefix).unwrap_or(&release);
    let tag = if prefix.is_empty() { format!("v{version}") } else { format!("v{prefix}{version}") };
    let name = format!("gencode.{tag}.annotation.gtf.gz");
    let url = format!("{root}/release_{release}/{name}");
    let path = genes_dir.join(&name);
    download(&url, &path, true)?;
    Ok(Some(GencodeAnnotation { release, path, url }))
}

pub fn fetch_from_base_with_gencode(assembly: &str, cache_root: &Path, base_url: &str, gencode_release: Option<&str>, use_gencode: bool) -> Result<PathBuf> {
    let root = cache_root.join(assembly);
    let genome_dir = root.join("Genome");
    let genes_dir = root.join("Genes");
    let protein_dir = root.join("Protein");
    std::fs::create_dir_all(&genome_dir)?;
    std::fs::create_dir_all(&genes_dir)?;
    std::fs::create_dir_all(&protein_dir)?;

    let bigzips = format!("{base_url}/goldenPath/{assembly}/bigZips");
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
    // Prefer GENCODE when the assembly has a GENCODE annotation.  `latest`
    // is resolved at build time to a concrete release; unsupported organisms
    // simply keep the UCSC/RefSeq annotation downloaded above.
    let gencode = if use_gencode {
        fetch_gencode_annotation(assembly, &genes_dir, gencode_release)?
    } else {
        None
    };
    if !have_gtf && gencode.is_none() && use_gencode {
        bail!("no supported gene annotation available for {assembly}; supply --gtf explicitly");
    }

    let uniprot_root = format!("{base_url}/goldenPath/archive/{assembly}/uniprot");
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

pub fn fetch_from_base(assembly: &str, cache_root: &Path, base_url: &str) -> Result<PathBuf> {
    fetch_from_base_with_gencode(assembly, cache_root, base_url, None, true)
}

#[allow(dead_code)]
pub fn fetch(assembly: &str, cache_root: &Path) -> Result<PathBuf> {
    fetch_from_base(assembly, cache_root, BASE_URL)
}
