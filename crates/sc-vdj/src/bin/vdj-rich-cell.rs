use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use std::cmp::Ordering;
use std::env;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct CellScore {
    column: usize,
    receptor_molecules: f64,
    total_molecules: f64,
}

#[derive(Debug)]
struct Cli {
    exonic: PathBuf,
    bam: PathBuf,
    out: PathBuf,
    cell: Option<String>,
    top: usize,
    genes: usize,
    cell_barcode_len: Option<usize>,
}

fn main() -> Result<()> {
    let cli = parse_cli(env::args_os().skip(1))?;
    let barcodes_path = find_file(&cli.exonic, &["barcodes.tsv.gz", "barcodes.tsv"])?;
    let features_path = find_file(
        &cli.exonic,
        &[
            "features.tsv.gz",
            "features.tsv",
            "genes.tsv.gz",
            "genes.tsv",
        ],
    )?;
    let matrix_path = find_file(&cli.exonic, &["matrix.mtx.gz", "matrix.mtx"])?;

    let barcodes = read_barcodes(&barcodes_path)?;
    let (feature_lines, feature_names, receptor_rows) = read_features(&features_path)?;
    let (scores, receptor_by_cell, matrix_entries) =
        score_matrix(&matrix_path, &receptor_rows, barcodes.len())?;

    let mut ranked: Vec<CellScore> = scores
        .into_iter()
        .enumerate()
        .filter_map(|(column, (receptor_molecules, total_molecules))| {
            (receptor_molecules > 0.0).then_some(CellScore {
                column,
                receptor_molecules,
                total_molecules,
            })
        })
        .collect();
    ranked.sort_by(|a, b| {
        b.receptor_molecules
            .partial_cmp(&a.receptor_molecules)
            .unwrap_or(Ordering::Equal)
            .then_with(|| {
                b.total_molecules
                    .partial_cmp(&a.total_molecules)
                    .unwrap_or(Ordering::Equal)
            })
    });
    if ranked.is_empty() {
        bail!("no cell has IG/TR expression in the exonic matrix");
    }

    println!(
        "matrix: {} features, {} cells, {} IG/TR features",
        feature_names.len(),
        barcodes.len(),
        receptor_rows.iter().filter(|&&x| x).count()
    );
    println!("\nTop receptor-rich cells:\nrank\tcell\tIG/TR molecules\ttotal molecules\tfraction");
    for (rank, cell) in ranked.iter().take(cli.top).enumerate() {
        let fraction = if cell.total_molecules > 0.0 {
            cell.receptor_molecules / cell.total_molecules
        } else {
            0.0
        };
        println!(
            "{}\t{}\t{:.0}\t{:.0}\t{:.5}",
            rank + 1,
            barcodes[cell.column],
            cell.receptor_molecules,
            cell.total_molecules,
            fraction
        );
    }

    let winner = if let Some(requested) = cli.cell.as_deref() {
        barcodes
            .iter()
            .position(|x| x == requested)
            .with_context(|| format!("requested cell {requested} is not present in matrix"))?
    } else {
        ranked[0].column
    };
    let matrix_barcode = &barcodes[winner];
    let bam_barcode = resolve_bam_barcode(&cli.bam, matrix_barcode, cli.cell_barcode_len)?;

    println!("\nMATRIX_CELL={matrix_barcode}");
    println!("BAM_CELL={bam_barcode}");
    if matrix_barcode != &bam_barcode {
        println!("NOTE: matrix barcode contains legacy u64 padding; fixture uses BAM CB barcode");
    }
    println!("\nReceptor genes:");
    let mut genes = receptor_by_cell[winner].clone();
    genes.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));
    for (row, value) in genes.iter().take(cli.genes) {
        println!("{}\t{:.0}", feature_names[*row], value);
    }

    fs::create_dir_all(&cli.out).with_context(|| format!("creating {}", cli.out.display()))?;
    let exonic_out = cli.out.join("exonic");
    fs::create_dir_all(&exonic_out)?;
    write_single_cell_mex(
        &exonic_out,
        &bam_barcode,
        &feature_lines,
        &matrix_entries[winner],
    )?;
    let bam_out = cli.out.join("cell.bam");
    let n_bam = write_cell_bam(&cli.bam, &bam_out, &bam_barcode)?;

    println!("\nfixture: {}", cli.out.display());
    println!("  BAM: {} records -> {}", n_bam, bam_out.display());
    println!(
        "  MEX: 1 cell, {} features -> {}",
        feature_lines.len(),
        exonic_out.display()
    );
    Ok(())
}

fn parse_cli<I: IntoIterator<Item = OsString>>(args: I) -> Result<Cli> {
    let mut exonic = None;
    let mut bam = None;
    let mut out = None;
    let mut cell = None;
    let mut top = 20;
    let mut genes = 40;
    let mut cell_barcode_len = None;
    let mut args = args.into_iter();
    while let Some(arg) = args.next() {
        match arg.to_string_lossy().as_ref() {
            "-h" | "--help" => {
                print_help();
                std::process::exit(0)
            }
            "--exonic" => exonic = Some(PathBuf::from(next_value(&mut args, "--exonic")?)),
            "--bam" => bam = Some(PathBuf::from(next_value(&mut args, "--bam")?)),
            "--out" => out = Some(PathBuf::from(next_value(&mut args, "--out")?)),
            "--cell" => {
                cell = Some(
                    next_value(&mut args, "--cell")?
                        .to_string_lossy()
                        .into_owned(),
                )
            }
            "--top" => top = parse_usize(next_value(&mut args, "--top")?, "--top")?,
            "--genes" => genes = parse_usize(next_value(&mut args, "--genes")?, "--genes")?,
            "--cell-barcode-len" => {
                let len = parse_usize(
                    next_value(&mut args, "--cell-barcode-len")?,
                    "--cell-barcode-len",
                )?;
                if !(1..=32).contains(&len) {
                    bail!("--cell-barcode-len must be between 1 and 32");
                }
                cell_barcode_len = Some(len);
            }
            other => bail!("unknown argument {other}\n\nRun vdj-rich-cell --help for usage."),
        }
    }
    Ok(Cli {
        exonic: exonic.context("missing --exonic <DIR>")?,
        bam: bam.context("missing --bam <BAM>")?,
        out: out.context("missing --out <DIR>")?,
        cell,
        top,
        genes,
        cell_barcode_len,
    })
}
fn next_value<I: Iterator<Item = OsString>>(args: &mut I, flag: &str) -> Result<OsString> {
    args.next()
        .with_context(|| format!("{flag} requires a value"))
}
fn parse_usize(v: OsString, flag: &str) -> Result<usize> {
    v.to_string_lossy()
        .parse()
        .with_context(|| format!("invalid value for {flag}"))
}
fn print_help() {
    println!("vdj-rich-cell - build a real one-cell Nelrune VDJ integration fixture\n\nUsage:\n  vdj-rich-cell --exonic <MEX_DIR> --bam <BAM> --out <DIR> [options]\n\nRequired:\n  --exonic <DIR>   Nelrune exonic MEX\n  --bam <FILE>     retained Nelrune mapper BAM\n  --out <DIR>      output fixture directory\n\nOptions:\n  --cell <BARCODE> use this matrix cell instead of richest IG/TR cell\n  --cell-barcode-len <N> use first N matrix-barcode bases for exact BAM CB matching\n  --top <N>        show top N cells [20]\n  --genes <N>      show top N receptor genes [40]\n  -h, --help");
}

fn find_file(dir: &Path, names: &[&str]) -> Result<PathBuf> {
    for n in names {
        let p = dir.join(n);
        if p.is_file() {
            return Ok(p);
        }
    }
    bail!("none of {} found in {}", names.join(", "), dir.display())
}
fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path)?;
    let r: Box<dyn IoRead> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(f))
    } else {
        Box::new(f)
    };
    Ok(Box::new(BufReader::new(r)))
}
fn read_barcodes(path: &Path) -> Result<Vec<String>> {
    Ok(reader(path)?
        .lines()
        .filter_map(|x| x.ok())
        .filter_map(|x| x.split('\t').next().map(str::to_string))
        .filter(|x| !x.is_empty())
        .collect())
}
fn is_receptor_gene(name: &str) -> bool {
    let u = name.to_ascii_uppercase();
    ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
        .iter()
        .any(|p| u.starts_with(p))
}
fn read_features(path: &Path) -> Result<(Vec<String>, Vec<String>, Vec<bool>)> {
    let mut lines = Vec::new();
    let mut names = Vec::new();
    let mut receptor = Vec::new();
    for line in reader(path)?.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let (name, gex) = {
            let f: Vec<&str> = line.split('\t').collect();
            let name = f.get(1).copied().unwrap_or(f[0]).to_string();
            let gex = f.get(2).map(|x| *x == "Gene Expression").unwrap_or(true);
            (name, gex)
        };
        receptor.push(gex && is_receptor_gene(&name));
        names.push(name);
        lines.push(line);
    }
    Ok((lines, names, receptor))
}

type CellGeneValues = Vec<Vec<(usize, f64)>>;
fn score_matrix(
    path: &Path,
    receptor_rows: &[bool],
    n_cells: usize,
) -> Result<(Vec<(f64, f64)>, CellGeneValues, CellGeneValues)> {
    let mut input = reader(path)?;
    let mut line = String::new();
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            bail!("matrix ended before dimensions")
        }
        if !line.starts_with('%') {
            break;
        }
    }
    let d: Vec<&str> = line.split_whitespace().collect();
    let rows: usize = d[0].parse()?;
    let cols: usize = d[1].parse()?;
    if rows != receptor_rows.len() || cols != n_cells {
        bail!("MEX dimensions do not match features/barcodes")
    }
    let mut scores = vec![(0.0, 0.0); n_cells];
    let mut receptor = vec![Vec::new(); n_cells];
    let mut entries = vec![Vec::new(); n_cells];
    for line in input.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut f = line.split_whitespace();
        let r: usize = f.next().context("matrix row")?.parse()?;
        let c: usize = f.next().context("matrix col")?.parse()?;
        let v: f64 = f.next().context("matrix value")?.parse()?;
        let r0 = r - 1;
        let c0 = c - 1;
        scores[c0].1 += v;
        entries[c0].push((r0, v));
        if receptor_rows[r0] {
            scores[c0].0 += v;
            receptor[c0].push((r0, v));
        }
    }
    Ok((scores, receptor, entries))
}

fn bam_cb(record: &bam::Record) -> Option<&str> {
    match record.aux(b"CB").ok()? {
        Aux::String(s) => Some(s),
        _ => None,
    }
}
fn bam_ub(record: &bam::Record) -> Option<&str> {
    match record.aux(b"UB").ok()? {
        Aux::String(s) => Some(s),
        _ => None,
    }
}

fn decode_hex_ascii(value: &str) -> Option<String> {
    if value.len() % 2 != 0 {
        return None;
    }
    let mut out = Vec::with_capacity(value.len() / 2);
    let bytes = value.as_bytes();
    for i in (0..bytes.len()).step_by(2) {
        let hi = (bytes[i] as char).to_digit(16)?;
        let lo = (bytes[i + 1] as char).to_digit(16)?;
        out.push(((hi << 4) | lo) as u8);
    }
    String::from_utf8(out).ok()
}

fn qname_identity(record: &bam::Record) -> Option<(String, String)> {
    let qname = std::str::from_utf8(record.qname()).ok()?;
    let fields: Vec<&str> = qname.split('|').collect();
    if fields.len() < 4 {
        return None;
    }
    let cell = decode_hex_ascii(fields[1])?;
    let umi = decode_hex_ascii(fields[3])?;
    Some((cell, umi))
}

fn record_identity(record: &bam::Record) -> Option<(String, String)> {
    if let Some(cb) = bam_cb(record) {
        return Some((cb.to_string(), bam_ub(record).unwrap_or("").to_string()));
    }
    qname_identity(record)
}

fn padded_matrix_matches_cb(matrix: &str, cb: &str) -> bool {
    matrix == cb
        || (matrix.len() > cb.len()
            && matrix.starts_with(cb)
            && matrix.as_bytes()[cb.len()..].iter().all(|&b| b == b'A'))
}
fn resolve_bam_barcode(
    path: &Path,
    matrix_barcode: &str,
    cell_barcode_len: Option<usize>,
) -> Result<String> {
    let requested = match cell_barcode_len {
        Some(len) => matrix_barcode.get(..len).with_context(|| {
            format!("matrix barcode {matrix_barcode} is shorter than --cell-barcode-len {len}")
        })?,
        None => matrix_barcode,
    };
    let mut r =
        bam::Reader::from_path(path).with_context(|| format!("opening {}", path.display()))?;
    for rec in r.records() {
        let rec = rec?;
        if let Some((cell, _)) = record_identity(&rec) {
            let matches = if cell_barcode_len.is_some() {
                cell == requested
            } else {
                padded_matrix_matches_cb(matrix_barcode, &cell)
            };
            if matches {
                return Ok(cell);
            }
        }
    }
    if let Some(len) = cell_barcode_len {
        bail!("could not resolve matrix barcode {matrix_barcode} using first {len} bases ({requested}) against BAM CB tags or encoded mapper QNAMEs")
    }
    bail!("could not resolve matrix barcode {matrix_barcode} against BAM CB tags or encoded mapper QNAMEs")
}

fn write_cell_bam(input: &Path, out: &Path, cell: &str) -> Result<usize> {
    let mut r = bam::Reader::from_path(input)?;
    let header = bam::Header::from_template(r.header());
    let mut w = bam::Writer::from_path(out, &header, bam::Format::Bam)?;
    let mut n = 0;
    for rec in r.records() {
        let mut rec = rec?;
        let Some((record_cell, umi)) = record_identity(&rec) else {
            continue;
        };
        if record_cell != cell {
            continue;
        }
        if bam_cb(&rec).is_none() {
            rec.push_aux(b"CB", Aux::String(&record_cell))?
        }
        if bam_ub(&rec).is_none() && !umi.is_empty() {
            rec.push_aux(b"UB", Aux::String(&umi))?
        }
        w.write(&rec)?;
        n += 1;
    }
    if n == 0 {
        bail!("no BAM records found for cell={cell}")
    }
    Ok(n)
}
fn gz_writer(path: &Path) -> Result<BufWriter<GzEncoder<File>>> {
    Ok(BufWriter::new(GzEncoder::new(
        File::create(path)?,
        Compression::default(),
    )))
}
fn write_single_cell_mex(
    out: &Path,
    barcode: &str,
    features: &[String],
    entries: &[(usize, f64)],
) -> Result<()> {
    let mut b = gz_writer(&out.join("barcodes.tsv.gz"))?;
    writeln!(b, "{barcode}")?;
    let mut f = gz_writer(&out.join("features.tsv.gz"))?;
    for line in features {
        writeln!(f, "{line}")?
    }
    let mut m = gz_writer(&out.join("matrix.mtx.gz"))?;
    writeln!(m, "%%MatrixMarket matrix coordinate real general")?;
    writeln!(m, "{} 1 {}", features.len(), entries.len())?;
    for (row, v) in entries {
        writeln!(m, "{} 1 {}", row + 1, v)?
    }
    Ok(())
}
