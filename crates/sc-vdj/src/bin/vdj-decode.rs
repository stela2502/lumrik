use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use sc_vdj::{
    DecodedNumericRecombinationId, DecodedRecombinationId, RecombinationId, VdjIndex,
};
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::str::FromStr;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "vdj-decode",
    about = "Decode and migrate compact sc-vdj recombination identifiers"
)]
struct Cli {
    /// Optional VDJ index. Required only for resolving legacy v1/v2 IDs to
    /// reference segments or for directory migration.
    #[arg(long)]
    index: Option<PathBuf>,

    /// Recombination IDs to decode. If omitted, IDs are read from stdin.
    #[arg(value_name = "CODE", num_args = 0..)]
    code: Vec<String>,

    /// After decoding, print only structural fields that differ across the supplied IDs.
    /// An index resolves segment names; index-free comparison uses numeric segment IDs.
    /// Legacy IDs with multiple numeric interpretations still require an index.
    #[arg(long)]
    compare: bool,

    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Copy an existing nelrune-vdj output directory and rewrite every
    /// HC:/LC: structural identifier to the self-describing v3 encoding.
    Migrate {
        /// VDJ index used by the original run.
        #[arg(long)]
        index: PathBuf,
        /// Existing nelrune-vdj output directory. It is never modified.
        #[arg(long)]
        input: PathBuf,
        /// New migrated output directory, e.g. vdj_v3.
        #[arg(long)]
        output: PathBuf,
    },
}

fn main() -> Result<()> {
    let c = Cli::parse();
    if let Some(Command::Migrate { index, input, output }) = c.command {
        return migrate_dir(&index, &input, &output);
    }

    let idx = c
        .index
        .as_ref()
        .map(|path| VdjIndex::load(path).context("loading VDJ index"))
        .transpose()?;
    let codes = if c.code.is_empty() {
        io::stdin()
            .lock()
            .lines()
            .collect::<std::result::Result<Vec<_>, _>>()?
    } else {
        c.code
    };

    if idx.is_none() {
        eprintln!(
            "warning: no VDJ index supplied; segment names are unavailable. \
             v3 IDs decode uniquely without an index; legacy v2 IDs may have multiple numeric candidates."
        );
    }

    for s in &codes {
        let id = RecombinationId::from_str(s.trim()).map_err(anyhow::Error::msg)?;
        if let Some(idx) = idx.as_ref() {
            let d = id.decode(idx).map_err(anyhow::Error::msg)?;
            println!(
                "{}\tchain={}\tV={}\tD={}\tJ={}\tv_del_3={}\tp_v3={}\tn1={}\td_del_5={}\td_retained={}\td_del_3={}\tp_d3={}\tn2={}\tj_del_5={}\tp_j5={}\tpn_alternative={}",
                id,
                d.chain,
                d.v,
                d.d.as_deref().unwrap_or(""),
                d.j,
                d.v_del_3,
                d.p_v3_len,
                d.n1_len,
                opt(d.d_del_5),
                opt(d.d_retained_len),
                opt(d.d_del_3),
                opt(d.p_d3_len),
                opt(d.n2_len),
                d.j_del_5,
                d.p_j5_len,
                d.pn_alternative
            );
        } else {
            let candidates = id.decode_numeric_candidates().map_err(anyhow::Error::msg)?;
            if candidates.len() > 1 {
                eprintln!(
                    "warning: {} is a legacy v2 ID with {} valid index-free V/D/J interpretations; printing all candidates",
                    id,
                    candidates.len()
                );
            }
            for (i, d) in candidates.iter().enumerate() {
                print_numeric(&id, d, i + 1, candidates.len());
            }
        }
    }

    if c.compare {
        if codes.len() < 2 {
            bail!("--compare requires at least two recombination IDs");
        }
        if let Some(idx) = idx.as_ref() {
            let decoded = codes
                .iter()
                .map(|s| {
                    let id = RecombinationId::from_str(s.trim()).map_err(anyhow::Error::msg)?;
                    let d = id.decode(idx).map_err(anyhow::Error::msg)?;
                    Ok((id, d))
                })
                .collect::<Result<Vec<_>>>()?;
            print_comparison(&decoded);
        } else {
            let decoded = codes
                .iter()
                .map(|s| {
                    let id = RecombinationId::from_str(s.trim()).map_err(anyhow::Error::msg)?;
                    let candidates = id.decode_numeric_candidates().map_err(anyhow::Error::msg)?;
                    if candidates.len() != 1 {
                        bail!(
                            "--compare without --index requires IDs with a unique numeric decode; {} has {} candidates",
                            id,
                            candidates.len()
                        );
                    }
                    Ok((id, candidates.into_iter().next().unwrap()))
                })
                .collect::<Result<Vec<_>>>()?;
            print_numeric_comparison(&decoded);
        }
    }
    Ok(())
}

fn print_numeric_comparison(decoded: &[(RecombinationId, DecodedNumericRecombinationId)]) {
    let rows: Vec<(&str, Vec<String>)> = vec![
        ("chain", decoded.iter().map(|(_, d)| d.chain.map(|x| x.to_string()).unwrap_or_default()).collect()),
        ("V_id", decoded.iter().map(|(_, d)| d.v_id.to_string()).collect()),
        ("D_id", decoded.iter().map(|(_, d)| d.d_id.map(|x| x.to_string()).unwrap_or_default()).collect()),
        ("J_id", decoded.iter().map(|(_, d)| d.j_id.to_string()).collect()),
        ("v_del_3", decoded.iter().map(|(_, d)| d.v_del_3.to_string()).collect()),
        ("p_v3", decoded.iter().map(|(_, d)| d.p_v3_len.to_string()).collect()),
        ("n1", decoded.iter().map(|(_, d)| d.n1_len.to_string()).collect()),
        ("p_d5", decoded.iter().map(|(_, d)| opt(d.p_d5_len)).collect()),
        ("d_del_5", decoded.iter().map(|(_, d)| opt(d.d_del_5)).collect()),
        ("d_retained", decoded.iter().map(|(_, d)| opt(d.d_retained_len)).collect()),
        ("d_del_3", decoded.iter().map(|(_, d)| opt(d.d_del_3)).collect()),
        ("p_d3", decoded.iter().map(|(_, d)| opt(d.p_d3_len)).collect()),
        ("n2", decoded.iter().map(|(_, d)| opt(d.n2_len)).collect()),
        ("j_del_5", decoded.iter().map(|(_, d)| d.j_del_5.to_string()).collect()),
        ("p_j5", decoded.iter().map(|(_, d)| d.p_j5_len.to_string()).collect()),
        ("pn_alternative", decoded.iter().map(|(_, d)| d.pn_alternative.to_string()).collect()),
    ];

    println!("comparison (differing fields only; numeric segment IDs)");
    print!("field");
    for (id, _) in decoded {
        print!("\t{id}");
    }
    println!();
    for (name, values) in rows {
        if values.windows(2).any(|w| w[0] != w[1]) {
            print!("{name}");
            for value in values {
                print!("\t{value}");
            }
            println!();
        }
    }
}

fn print_comparison(decoded: &[(RecombinationId, DecodedRecombinationId)]) {
    let rows: Vec<(&str, Vec<String>)> = vec![
        ("chain", decoded.iter().map(|(_, d)| d.chain.to_string()).collect()),
        ("V", decoded.iter().map(|(_, d)| d.v.clone()).collect()),
        ("D", decoded.iter().map(|(_, d)| d.d.clone().unwrap_or_default()).collect()),
        ("J", decoded.iter().map(|(_, d)| d.j.clone()).collect()),
        ("v_del_3", decoded.iter().map(|(_, d)| d.v_del_3.to_string()).collect()),
        ("p_v3", decoded.iter().map(|(_, d)| d.p_v3_len.to_string()).collect()),
        ("n1", decoded.iter().map(|(_, d)| d.n1_len.to_string()).collect()),
        ("p_d5", decoded.iter().map(|(_, d)| opt(d.p_d5_len)).collect()),
        ("d_del_5", decoded.iter().map(|(_, d)| opt(d.d_del_5)).collect()),
        ("d_retained", decoded.iter().map(|(_, d)| opt(d.d_retained_len)).collect()),
        ("d_del_3", decoded.iter().map(|(_, d)| opt(d.d_del_3)).collect()),
        ("p_d3", decoded.iter().map(|(_, d)| opt(d.p_d3_len)).collect()),
        ("n2", decoded.iter().map(|(_, d)| opt(d.n2_len)).collect()),
        ("j_del_5", decoded.iter().map(|(_, d)| d.j_del_5.to_string()).collect()),
        ("p_j5", decoded.iter().map(|(_, d)| d.p_j5_len.to_string()).collect()),
        (
            "pn_alternative",
            decoded.iter().map(|(_, d)| d.pn_alternative.to_string()).collect(),
        ),
    ];

    println!("comparison (differing fields only)");
    print!("field");
    for (id, _) in decoded {
        print!("\t{id}");
    }
    println!();
    for (name, values) in rows {
        if values.windows(2).any(|w| w[0] != w[1]) {
            print!("{name}");
            for value in values {
                print!("\t{value}");
            }
            println!();
        }
    }
}

fn migrate_dir(index_path: &Path, input: &Path, output: &Path) -> Result<()> {
    if !input.is_dir() {
        bail!("migration input is not a directory: {}", input.display());
    }
    if output.exists() {
        bail!(
            "migration output already exists: {} (refusing to overwrite the original or a previous migration)",
            output.display()
        );
    }
    let index = VdjIndex::load(index_path)
        .with_context(|| format!("loading VDJ index {}", index_path.display()))?;
    let mut stats = MigrationStats::default();
    copy_and_migrate_tree(input, output, &index, &mut stats)?;

    let provenance = format!(
        "Lumrik VDJ recombination-ID migration\nsource={}\nindex={}\ntarget_recombination_id_format=3\nfiles_copied={}\ntext_files_rewritten={}\nidentifiers_rewritten={}\n",
        input.display(),
        index_path.display(),
        stats.files_copied,
        stats.text_files_rewritten,
        stats.identifiers_rewritten
    );
    fs::write(output.join("recombination_id_migration.txt"), provenance)?;
    eprintln!(
        "migrated {} -> {}: {} files copied, {} text files changed, {} HC/LC IDs rewritten to v3",
        input.display(),
        output.display(),
        stats.files_copied,
        stats.text_files_rewritten,
        stats.identifiers_rewritten
    );
    Ok(())
}

#[derive(Default)]
struct MigrationStats {
    files_copied: usize,
    text_files_rewritten: usize,
    identifiers_rewritten: usize,
}

fn copy_and_migrate_tree(
    input: &Path,
    output: &Path,
    index: &VdjIndex,
    stats: &mut MigrationStats,
) -> Result<()> {
    fs::create_dir(output)
        .with_context(|| format!("creating migration directory {}", output.display()))?;
    for entry in fs::read_dir(input).with_context(|| format!("reading {}", input.display()))? {
        let entry = entry?;
        let src = entry.path();
        let dst = output.join(entry.file_name());
        let ty = entry.file_type()?;
        if ty.is_dir() {
            copy_and_migrate_tree(&src, &dst, index, stats)?;
        } else if ty.is_file() {
            migrate_file(&src, &dst, index, stats)?;
        } else if ty.is_symlink() {
            bail!("refusing to migrate symlink {}", src.display());
        }
    }
    Ok(())
}

fn migrate_file(src: &Path, dst: &Path, index: &VdjIndex, stats: &mut MigrationStats) -> Result<()> {
    let bytes = fs::read(src).with_context(|| format!("reading {}", src.display()))?;
    stats.files_copied += 1;
    if let Ok(text) = std::str::from_utf8(&bytes) {
        let (rewritten, count) = rewrite_ids(text, index)
            .with_context(|| format!("rewriting recombination IDs in {}", src.display()))?;
        if count > 0 {
            fs::write(dst, rewritten)?;
            stats.text_files_rewritten += 1;
            stats.identifiers_rewritten += count;
            return Ok(());
        }
    }
    fs::write(dst, bytes)?;
    Ok(())
}

fn rewrite_ids(text: &str, index: &VdjIndex) -> Result<(String, usize)> {
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    let mut pos = 0usize;
    let mut count = 0usize;

    while pos < bytes.len() {
        let is_id = (bytes[pos..].starts_with(b"HC:") || bytes[pos..].starts_with(b"LC:"))
            && pos + 3 < bytes.len()
            && bytes[pos + 3].is_ascii_hexdigit();
        if !is_id {
            let ch = text[pos..].chars().next().expect("valid UTF-8 character");
            out.push(ch);
            pos += ch.len_utf8();
            continue;
        }

        let start = pos;
        pos += 3;
        while pos < bytes.len() && bytes[pos].is_ascii_hexdigit() {
            pos += 1;
        }
        let token = &text[start..pos];
        let old = RecombinationId::from_str(token).map_err(anyhow::Error::msg)?;
        let new = old.to_v3(index).map_err(anyhow::Error::msg)?;
        out.push_str(&new.to_string());
        count += usize::from(new != old);
    }
    Ok((out, count))
}

fn print_numeric(
    id: &RecombinationId,
    d: &DecodedNumericRecombinationId,
    candidate: usize,
    candidate_count: usize,
) {
    let chain = d.chain.map(|x| x.to_string()).unwrap_or_default();
    let candidate_field = if candidate_count > 1 {
        format!("\tcandidate={candidate}/{candidate_count}")
    } else {
        String::new()
    };
    println!(
        "{}{}\tchain={}\tV_id={}\tD_id={}\tJ_id={}\tv_del_3={}\tp_v3={}\tn1={}\td_del_5={}\td_retained={}\td_del_3={}\tp_d3={}\tn2={}\tj_del_5={}\tp_j5={}\tpn_alternative={}",
        id,
        candidate_field,
        chain,
        d.v_id,
        d.d_id.map(|x| x.to_string()).unwrap_or_default(),
        d.j_id,
        d.v_del_3,
        d.p_v3_len,
        d.n1_len,
        opt(d.d_del_5),
        opt(d.d_retained_len),
        opt(d.d_del_3),
        opt(d.p_d3_len),
        opt(d.n2_len),
        d.j_del_5,
        d.p_j5_len,
        d.pn_alternative
    );
}

fn opt(x: Option<u16>) -> String {
    x.map(|v| v.to_string()).unwrap_or_default()
}
