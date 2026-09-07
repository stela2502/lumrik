use anyhow::{bail, Context, Result};
use clap::Parser;
use sc_vdj::{PackedRecombinationId, RecombinationMeasurements, VdjMapper};
use std::io::{self, BufRead};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about = "Decode Lumrik V(D)J recombination identifiers")]
struct Cli {
    /// Compiled VDJ index used to resolve the packed V/D/J segment identifiers.
    #[arg(long, value_name = "VDJ_INDEX")]
    index: PathBuf,

    /// HC:/LC: identifier, or input from stdin when omitted.
    #[arg(value_name = "CODE", num_args = 0.., allow_hyphen_values = true)]
    code: Vec<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let reference = VdjMapper::load_index(&cli.index)
        .with_context(|| format!("loading V(D)J index {}", cli.index.display()))?;

    if !cli.code.is_empty() {
        decode_line(&cli.code.join(" "), &reference)?;
        return Ok(());
    }

    let stdin = io::stdin();
    let mut seen = false;
    for line in stdin.lock().lines() {
        let line = line?;
        let text = line.trim();
        if text.is_empty() {
            continue;
        }
        if seen {
            println!();
        }
        decode_line(text, &reference)?;
        seen = true;
    }
    if !seen {
        bail!("no code supplied; pass an HC:/LC: identifier as an argument or on stdin");
    }
    Ok(())
}

fn decode_line(text: &str, mapper: &VdjMapper) -> Result<()> {
    if !(text.starts_with("HC:") || text.starts_with("LC:")) {
        bail!("unrecognized VDJ identifier: expected HC:<hex> or LC:<hex>");
    }
    let id: PackedRecombinationId = text.parse().map_err(anyhow::Error::msg)?;
    let measurements = id.decode(mapper.reference()).map_err(anyhow::Error::msg)?;
    println!("type\tpacked_recombination");
    println!("id\t{id}");
    print_measurements("", &measurements);
    Ok(())
}

fn print_measurements(prefix: &str, x: &RecombinationMeasurements) {
    field(prefix, "chain", x.chain.to_string());
    field(prefix, "v", &x.v);
    measurement(prefix, "v_del_3", x.v_del_3);
    measurement(prefix, "p_v3_len", x.p_v3_len);
    measurement(prefix, "n1_len", x.n1_len);
    if x.chain.has_d() {
        measurement(prefix, "p_d5_len", x.p_d5_len);
        field(prefix, "d", x.d.as_deref().unwrap_or("?"));
        measurement(prefix, "d_del_5", x.d_del_5);
        measurement(prefix, "d_retained_len", x.d_retained_len);
        measurement(prefix, "d_del_3", x.d_del_3);
        measurement(prefix, "p_d3_len", x.p_d3_len);
        measurement(prefix, "n2_len", x.n2_len);
    }
    measurement(prefix, "p_j5_len", x.p_j5_len);
    measurement(prefix, "j_del_5", x.j_del_5);
    field(prefix, "j", &x.j);
    field(prefix, "pn_alternative", x.pn_alternative.to_string());
    field(prefix, "complete", x.is_complete().to_string());
}

fn measurement(prefix: &str, name: &str, value: Option<u16>) {
    field(prefix, name, value.map_or_else(|| "?".to_string(), |x| x.to_string()));
}

fn field(prefix: &str, name: &str, value: impl std::fmt::Display) {
    println!("{prefix}{name}\t{value}");
}
