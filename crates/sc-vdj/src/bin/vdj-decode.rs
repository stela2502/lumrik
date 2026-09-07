use anyhow::{Context, Result};
use clap::Parser;
use sc_vdj::{RecombinationId, VdjIndex};
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::str::FromStr;
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "vdj-decode",
    about = "Decode compact sc-vdj recombination identifiers"
)]
struct Cli {
    #[arg(long)]
    index: PathBuf,
    #[arg(value_name="CODE",num_args=0..)]
    code: Vec<String>,
}
fn main() -> Result<()> {
    let c = Cli::parse();
    let idx = VdjIndex::load(&c.index).context("loading VDJ index")?;
    let codes = if c.code.is_empty() {
        io::stdin()
            .lock()
            .lines()
            .collect::<std::result::Result<Vec<_>, _>>()?
    } else {
        c.code
    };
    for s in codes {
        let id = RecombinationId::from_str(s.trim()).map_err(anyhow::Error::msg)?;
        let d = id.decode(&idx).map_err(anyhow::Error::msg)?;
        println!("{}\tchain={}\tV={}\tD={}\tJ={}\tv_del_3={}\tp_v3={}\tn1={}\td_del_5={}\td_retained={}\td_del_3={}\tp_d3={}\tn2={}\tj_del_5={}\tp_j5={}\tpn_alternative={}",id,d.chain,d.v,d.d.as_deref().unwrap_or(""),d.j,d.v_del_3,d.p_v3_len,d.n1_len,opt(d.d_del_5),opt(d.d_retained_len),opt(d.d_del_3),opt(d.p_d3_len),opt(d.n2_len),d.j_del_5,d.p_j5_len,d.pn_alternative)
    }
    Ok(())
}
fn opt(x: Option<u16>) -> String {
    x.map(|v| v.to_string()).unwrap_or_default()
}
