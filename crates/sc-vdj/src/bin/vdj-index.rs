use anyhow::{bail, Context, Result};
use sc_vdj::{VdjMapper, VdjMapperConfig, VdjReferenceBuilder};
use std::env;
use std::ffi::OsString;
use std::path::PathBuf;

fn main() -> Result<()> {
    let mut gtf = None;
    let mut genome = None;
    let mut out = None;
    let mut args = env::args_os().skip(1);
    while let Some(arg) = args.next() {
        match arg.to_string_lossy().as_ref() {
            "-h" | "--help" => {
                println!("vdj-index - compile a reusable lumrik V(D)J germline + seed index\n\nUsage:\n  vdj-index --gtf <GTF> --genome <FASTA> --out <FILE.vdjidx>");
                return Ok(());
            }
            "--gtf" => gtf = Some(PathBuf::from(next_value(&mut args, "--gtf")?)),
            "--genome" => genome = Some(PathBuf::from(next_value(&mut args, "--genome")?)),
            "--out" => out = Some(PathBuf::from(next_value(&mut args, "--out")?)),
            other => bail!("unknown argument {other}"),
        }
    }
    let gtf = gtf.context("missing --gtf <GTF>")?;
    let genome = genome.context("missing --genome <FASTA>")?;
    let out = out.context("missing --out <FILE.vdjidx>")?;
    eprintln!("building V(D)J germline reference...");
    let reference = VdjReferenceBuilder::default().build(&gtf, &genome)?;
    eprintln!("building ambiguity-aware coverage + identity seed indices for {} segments...", reference.len());
    let mapper = VdjMapper::new(reference, VdjMapperConfig::default());
    mapper.save_index(&out).with_context(|| format!("writing {}", out.display()))?;
    println!("VDJ index v4 written to {}", out.display());
    Ok(())
}

fn next_value<I: Iterator<Item = OsString>>(args: &mut I, flag: &str) -> Result<OsString> {
    args.next().with_context(|| format!("{flag} requires a value"))
}
