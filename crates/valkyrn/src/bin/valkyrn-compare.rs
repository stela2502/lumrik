use anyhow::{bail, Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(name = "valkyrn-compare", about = "Compare structural receptor identities across Valkyrn runs")]
struct Cli {
    /// Sample specification NAME=VALKYRN_DIR. Repeat for two or more samples.
    #[arg(long = "sample", required = true)]
    samples: Vec<String>,
    /// Output directory.
    #[arg(long, default_value = "valkyrn_compare")]
    out: PathBuf,
}

#[derive(Clone, Default)]
struct Cell {
    cell: String,
    productive_heavy: usize,
    productive_light: usize,
    paired_productive: bool,
    heavy_ids: Vec<String>,
    light_ids: Vec<String>,
    paired_ids: Vec<String>,
    hc_depth: Option<f64>,
    lc_depth: Option<f64>,
    paired_depth: Option<f64>,
}

#[derive(Default)]
struct Sample {
    cells: Vec<Cell>,
    heavy_cells: BTreeMap<String, BTreeSet<String>>,
    pair_cells: BTreeMap<String, BTreeSet<String>>,
}

fn column(header: &[&str], name: &str, path: &Path) -> Result<usize> {
    header.iter().position(|x| *x == name).with_context(|| {
        format!(
            "{} lacks column {name}; rerun Valkyrn with the current cross-sample output format",
            path.display()
        )
    })
}

fn split_ids(s: &str) -> Vec<String> {
    s.split(',').filter(|x| !x.is_empty()).map(str::to_string).collect()
}

fn parse_opt_f64(s: &str) -> Option<f64> {
    if s.is_empty() { None } else { s.parse().ok() }
}

fn read_sample(dir: &Path) -> Result<Sample> {
    let path = dir.join("valkyrn_cells.tsv");
    let file = File::open(&path).with_context(|| format!("opening {}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines.next().context("valkyrn_cells.tsv is empty")??;
    let h: Vec<&str> = header.split('\t').collect();
    let cell_i = column(&h, "cell", &path)?;
    let ph_i = column(&h, "productive_heavy", &path)?;
    let pl_i = column(&h, "productive_light", &path)?;
    let pp_i = column(&h, "paired_productive", &path)?;
    let heavy_i = column(&h, "heavy_ids", &path)?;
    let light_i = column(&h, "light_ids", &path)?;
    let pair_i = column(&h, "paired_ids", &path)?;
    let hd_i = column(&h, "hc_depth_median_nt", &path)?;
    let ld_i = column(&h, "lc_depth_median_nt", &path)?;
    let pd_i = column(&h, "paired_depth_median_nt", &path)?;
    let mut out = Sample::default();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        let c = Cell {
            cell: f.get(cell_i).copied().unwrap_or("").to_string(),
            productive_heavy: f.get(ph_i).copied().unwrap_or("0").parse().unwrap_or(0),
            productive_light: f.get(pl_i).copied().unwrap_or("0").parse().unwrap_or(0),
            paired_productive: f.get(pp_i).copied().unwrap_or("false") == "true",
            heavy_ids: split_ids(f.get(heavy_i).copied().unwrap_or("")),
            light_ids: split_ids(f.get(light_i).copied().unwrap_or("")),
            paired_ids: split_ids(f.get(pair_i).copied().unwrap_or("")),
            hc_depth: parse_opt_f64(f.get(hd_i).copied().unwrap_or("")),
            lc_depth: parse_opt_f64(f.get(ld_i).copied().unwrap_or("")),
            paired_depth: parse_opt_f64(f.get(pd_i).copied().unwrap_or("")),
        };
        for id in &c.heavy_ids { out.heavy_cells.entry(id.clone()).or_default().insert(c.cell.clone()); }
        for id in &c.paired_ids { out.pair_cells.entry(id.clone()).or_default().insert(c.cell.clone()); }
        out.cells.push(c);
    }
    Ok(out)
}

fn parse_spec(s: &str) -> Result<(String, PathBuf)> {
    let Some((name, path)) = s.split_once('=') else { bail!("--sample must be NAME=VALKYRN_DIR, got {s}"); };
    if name.is_empty() || path.is_empty() { bail!("invalid --sample {s}"); }
    Ok((name.to_string(), PathBuf::from(path)))
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() { return None; }
    v.sort_by(f64::total_cmp);
    let n = v.len();
    Some(if n % 2 == 1 { v[n/2] } else { (v[n/2-1] + v[n/2]) / 2.0 })
}
fn fmt(x: Option<f64>) -> String { x.map(|v| format!("{v:.2}")).unwrap_or_default() }
fn pct(n: usize, d: usize) -> f64 { if d == 0 { 0.0 } else { 100.0 * n as f64 / d as f64 } }

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.samples.len() < 2 { bail!("provide at least two --sample NAME=VALKYRN_DIR arguments"); }
    fs::create_dir_all(&cli.out)?;
    let mut samples = Vec::new();
    let mut names = BTreeSet::new();
    for spec in &cli.samples {
        let (name, dir) = parse_spec(spec)?;
        if !names.insert(name.clone()) { bail!("duplicate sample name: {name}"); }
        samples.push((name, read_sample(&dir)?));
    }

    let all_hc: BTreeSet<String> = samples.iter().flat_map(|(_,s)| s.heavy_cells.keys().cloned()).collect();
    let all_pairs: BTreeSet<String> = samples.iter().flat_map(|(_,s)| s.pair_cells.keys().cloned()).collect();

    // Pairwise table: both symmetric intersection and directional recovery.
    let mut pw = BufWriter::new(File::create(cli.out.join("pairwise.tsv"))?);
    writeln!(pw, "sample_a\tsample_b\ta_cells\tb_cells\ta_hc\tb_hc\tshared_hc\ta_hc_found_in_b_pct\tb_hc_found_in_a_pct\ta_paired\tb_paired\tshared_paired\ta_paired_found_in_b_pct\tb_paired_found_in_a_pct")?;
    eprintln!("\nPairwise structural overlap");
    eprintln!("sample A <-> sample B\tHC shared\tA in B\tB in A\tHC+LC shared\tA in B\tB in A");
    for a in 0..samples.len() {
        for b in a+1..samples.len() {
            let (an, av) = &samples[a]; let (bn, bv) = &samples[b];
            let sh = av.heavy_cells.keys().filter(|x| bv.heavy_cells.contains_key(*x)).count();
            let sp = av.pair_cells.keys().filter(|x| bv.pair_cells.contains_key(*x)).count();
            writeln!(pw, "{an}\t{bn}\t{}\t{}\t{}\t{}\t{sh}\t{:.2}\t{:.2}\t{}\t{}\t{sp}\t{:.2}\t{:.2}",
                av.cells.len(), bv.cells.len(), av.heavy_cells.len(), bv.heavy_cells.len(),
                pct(sh,av.heavy_cells.len()), pct(sh,bv.heavy_cells.len()), av.pair_cells.len(), bv.pair_cells.len(),
                pct(sp,av.pair_cells.len()), pct(sp,bv.pair_cells.len()))?;
            eprintln!("{an} <-> {bn}\t{sh}\t{:.1}%\t{:.1}%\t{sp}\t{:.1}%\t{:.1}%",
                pct(sh,av.heavy_cells.len()), pct(sh,bv.heavy_cells.len()), pct(sp,av.pair_cells.len()), pct(sp,bv.pair_cells.len()));
        }
    }

    // One rich row per structural HC or exact HC+LC identity.
    let mut rw = BufWriter::new(File::create(cli.out.join("receptor_summary.tsv"))?);
    write!(rw, "identity_type\treceptor_id\theavy_id\tlight_id\tn_samples\tsamples\ttotal_cells")?;
    for (name, _) in &samples { write!(rw, "\t{name}_cells\t{name}_fraction_of_receptor_pct\t{name}_depth_median_nt")?; }
    writeln!(rw)?;
    for (kind, ids) in [("HC", &all_hc), ("HC+LC", &all_pairs)] {
        for id in ids {
            let present: Vec<&str> = samples.iter().filter(|(_,s)| if kind=="HC" { s.heavy_cells.contains_key(id) } else { s.pair_cells.contains_key(id) }).map(|(n,_)|n.as_str()).collect();
            let total: usize = samples.iter().map(|(_,s)| if kind=="HC" { s.heavy_cells.get(id).map_or(0,BTreeSet::len) } else { s.pair_cells.get(id).map_or(0,BTreeSet::len) }).sum();
            let (heavy, light) = if kind == "HC+LC" { id.split_once("+LC:").map(|(h,l)|(h.to_string(),format!("LC:{l}"))).unwrap_or((id.clone(),String::new())) } else { (id.clone(),String::new()) };
            write!(rw, "{kind}\t{id}\t{heavy}\t{light}\t{}\t{}\t{total}", present.len(), present.join(","))?;
            for (_,s) in &samples {
                let n = if kind=="HC" { s.heavy_cells.get(id).map_or(0,BTreeSet::len) } else { s.pair_cells.get(id).map_or(0,BTreeSet::len) };
                let depths: Vec<f64> = s.cells.iter().filter(|c| if kind=="HC" { c.heavy_ids.contains(id) } else { c.paired_ids.contains(id) }).filter_map(|c| if kind=="HC" { c.hc_depth } else { c.paired_depth }).collect();
                write!(rw, "\t{n}\t{:.2}\t{}", pct(n,total), fmt(median(depths)))?;
            }
            writeln!(rw)?;
        }
    }

    // Scanpy/Seurat-ready cell metadata. Sample + cell is the globally safe key.
    let mut cw = BufWriter::new(File::create(cli.out.join("cell_annotations.tsv"))?);
    write!(cw, "cell_key\tsample\tcell\tproductive_heavy\tproductive_light\tpaired_productive\theavy_ids\tlight_ids\tpaired_ids\thc_depth_median_nt\tlc_depth_median_nt\tpaired_depth_median_nt\thc_max_n_samples\tpaired_max_n_samples\thc_cross_sample\tpaired_cross_sample")?;
    for (name, _) in &samples { write!(cw, "\thc_in_{name}\thc_cells_{name}\tpaired_in_{name}\tpaired_cells_{name}")?; }
    writeln!(cw)?;
    for (sn,s) in &samples {
        for c in &s.cells {
            let hc_ns = c.heavy_ids.iter().map(|id| samples.iter().filter(|(_,x)| x.heavy_cells.contains_key(id)).count()).max().unwrap_or(0);
            let pair_ns = c.paired_ids.iter().map(|id| samples.iter().filter(|(_,x)| x.pair_cells.contains_key(id)).count()).max().unwrap_or(0);
            write!(cw, "{}:{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                sn,c.cell,sn,c.cell,c.productive_heavy,c.productive_light,c.paired_productive,c.heavy_ids.join(","),c.light_ids.join(","),c.paired_ids.join(","),fmt(c.hc_depth),fmt(c.lc_depth),fmt(c.paired_depth),hc_ns,pair_ns,hc_ns>1,pair_ns>1)?;
            for (_,other) in &samples {
                let hcount: usize = c.heavy_ids.iter().map(|id| other.heavy_cells.get(id).map_or(0,BTreeSet::len)).sum();
                let pcount: usize = c.paired_ids.iter().map(|id| other.pair_cells.get(id).map_or(0,BTreeSet::len)).sum();
                write!(cw, "\t{}\t{}\t{}\t{}", hcount>0,hcount,pcount>0,pcount)?;
            }
            writeln!(cw)?;
        }
    }

    let shared_hc = all_hc.iter().filter(|id| samples.iter().filter(|(_,s)| s.heavy_cells.contains_key(*id)).count()>1).count();
    let shared_pairs = all_pairs.iter().filter(|id| samples.iter().filter(|(_,s)| s.pair_cells.contains_key(*id)).count()>1).count();
    let three_hc = all_hc.iter().filter(|id| samples.iter().all(|(_,s)| s.heavy_cells.contains_key(*id))).count();
    let three_pairs = all_pairs.iter().filter(|id| samples.iter().all(|(_,s)| s.pair_cells.contains_key(*id))).count();
    let mut sw = BufWriter::new(File::create(cli.out.join("summary.tsv"))?);
    writeln!(sw,"sample\tcells\tpaired_productive_cells\tunique_heavy_ids\tunique_paired_ids\thc_depth_median_nt\tlc_depth_median_nt\tpaired_depth_median_nt")?;
    for (name,s) in &samples {
        writeln!(sw,"{name}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",s.cells.len(),s.cells.iter().filter(|c|c.paired_productive).count(),s.heavy_cells.len(),s.pair_cells.len(),fmt(median(s.cells.iter().filter_map(|c|c.hc_depth).collect())),fmt(median(s.cells.iter().filter_map(|c|c.lc_depth).collect())),fmt(median(s.cells.iter().filter_map(|c|c.paired_depth).collect())))?;
    }
    writeln!(sw,"#union\t\t\t{}\t{}\t\t\t",all_hc.len(),all_pairs.len())?;
    writeln!(sw,"#shared_in_2plus_samples\t\t\t{shared_hc}\t{shared_pairs}\t\t\t")?;
    writeln!(sw,"#shared_in_all_samples\t\t\t{three_hc}\t{three_pairs}\t\t\t")?;
    eprintln!("\nUnion: {} HC ({} shared in >=2; {} in all), {} HC+LC ({} shared in >=2; {} in all)",all_hc.len(),shared_hc,three_hc,all_pairs.len(),shared_pairs,three_pairs);
    eprintln!("Wrote summary.tsv, pairwise.tsv, receptor_summary.tsv and cell_annotations.tsv to {}",cli.out.display());
    Ok(())
}
