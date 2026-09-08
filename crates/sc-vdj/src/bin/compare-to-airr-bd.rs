use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use sc_primer::{Chemistry, PrimerDetector};

const DEFAULT_LUMRIK: &str =
    "target/sc-vdj-test-output/server-ZD-4631-BcellsLaneF/vdj_calls.tsv";

const DEFAULT_BD_AIRR: &str =
    "/home/med-sal/sens05_shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/db_pipeline_B_cell_run/B_cells_VDJ_Dominant_Contigs_AIRR.tsv";

#[derive(Debug, Parser)]
#[command(
    name = "compare-to-airr-bd",
    about = "Disposable development comparator for Lumrik sc-vdj versus BD Rhapsody AIRR"
)]
struct Args {
    /// 	Lumrik vdj_calls.tsv.
    #[arg(long, default_value = DEFAULT_LUMRIK)]
    lumrik: PathBuf,

    /// 	BD Rhapsody dominant AIRR TSV.
    #[arg(long, default_value = DEFAULT_BD_AIRR)]
    bd_airr: PathBuf,

    /// 	Single-cell chemistry used to translate Lumrik barcode sequences
    /// 	to the numeric BD cell IDs used in the AIRR file.
    #[arg(long, value_enum, default_value = "bd-v2-384")]
    chemistry: Chemistry,
}

#[derive(Debug, Clone)]
struct LumrikCall {
    cell_id: u64,
    chain: String,
    v: String,
    d: String,
    j: String,
    c: String,
}

#[derive(Debug, Clone)]
struct BdCall {
    cell_id: u64,
    locus: String,
    v: String,
    d: String,
    j: String,
    c: String,
    junction: String,
    junction_aa: String,
    productive: String,
    umi_count: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let detector = PrimerDetector::from_chemistry(args.chemistry)
        .map_err(|e| anyhow::anyhow!("creating sc-primer detector: {e}"))?;

    let lumrik = read_lumrik(&args.lumrik, &detector)?;
    let bd = read_bd_airr(&args.bd_airr)?;

    println!("Input");
    println!("  Lumrik: {}", args.lumrik.display());
    println!("  BD AIRR: {}", args.bd_airr.display());
    println!("  Chemistry: {:?}", args.chemistry);
    println!();

    let lumrik_cells: HashSet<u64> =
        lumrik.iter().map(|x| x.cell_id).collect();

    let bd_cells: HashSet<u64> =
        bd.iter().map(|x| x.cell_id).collect();

    let shared_cells =
        lumrik_cells.intersection(&bd_cells).count();

    let lumrik_only_cells =
        lumrik_cells.difference(&bd_cells).count();

    let bd_only_cells =
        bd_cells.difference(&lumrik_cells).count();

    println!("Calls");
    println!("  Lumrik: {}", lumrik.len());
    println!("  BD AIRR: {}", bd.len());
    println!();

    println!("Cells");
    println!("  Lumrik: {}", lumrik_cells.len());
    println!("  BD AIRR: {}", bd_cells.len());
    println!("  Shared: {}", shared_cells);
    println!("  Lumrik only: {}", lumrik_only_cells);
    println!("  BD only: {}", bd_only_cells);
    println!();

    report_lumrik_loci(&lumrik);
    report_bd_loci(&bd);

    compare_calls(&lumrik, &bd);

    Ok(())
}

fn read_lumrik(
    path: &Path,
    detector: &PrimerDetector,
) -> Result<Vec<LumrikCall>> {
    let file = File::open(path)
        .with_context(|| format!("opening Lumrik TSV {}", path.display()))?;

    let mut lines = BufReader::new(file).lines();

    let header = lines
        .next()
        .context("Lumrik TSV is empty")??;

    let columns: Vec<&str> = header.split('\t').collect();

    let cell_idx = column_index(&columns, "cell")?;
    let chain_idx = column_index(&columns, "chain")?;
    let v_idx = column_index(&columns, "v")?;
    let d_idx = column_index(&columns, "d")?;
    let j_idx = column_index(&columns, "j")?;
    let c_idx = column_index(&columns, "c")?;

    let mut calls = Vec::new();
    let mut unmapped_barcodes = HashSet::new();

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        let cell = field_ref(&fields, cell_idx);

        let Some(cell_id) =
            detector.cell_id_for_seq(cell.as_bytes())
        else {
            unmapped_barcodes.insert(cell.to_string());
            continue;
        };

        calls.push(LumrikCall {
            cell_id,
            chain: field(&fields, chain_idx),
            v: normalize_gene(field_ref(&fields, v_idx)),
            d: normalize_gene(field_ref(&fields, d_idx)),
            j: normalize_gene(field_ref(&fields, j_idx)),
            c: normalize_gene(field_ref(&fields, c_idx)),
        });
    }

    if !unmapped_barcodes.is_empty() {
        eprintln!(
            "Warning: {} Lumrik cell barcodes could not be translated to numeric BD cell IDs",
            unmapped_barcodes.len()
        );
    }

    Ok(calls)
}

fn read_bd_airr(path: &Path) -> Result<Vec<BdCall>> {
    let file = File::open(path)
        .with_context(|| format!("opening BD AIRR TSV {}", path.display()))?;

    let mut lines = BufReader::new(file).lines();

    let header = lines
        .next()
        .context("BD AIRR TSV is empty")??;

    let columns: Vec<&str> = header.split('\t').collect();

    let cell_idx = column_index(&columns, "cell_id")?;
    let locus_idx = column_index(&columns, "locus")?;
    let v_idx = column_index(&columns, "v_call")?;
    let d_idx = column_index(&columns, "d_call")?;
    let j_idx = column_index(&columns, "j_call")?;
    let c_idx = column_index(&columns, "c_call")?;
    let junction_idx = column_index(&columns, "junction")?;
    let junction_aa_idx = column_index(&columns, "junction_aa")?;
    let productive_idx = column_index(&columns, "productive")?;
    let umi_count_idx = column_index(&columns, "umi_count")?;

    let mut calls = Vec::new();

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();

        let raw_cell_id = field_ref(&fields, cell_idx);

        let cell_id: u64 = raw_cell_id
            .parse()
            .with_context(|| {
                format!("invalid BD cell_id {raw_cell_id:?}")
            })?;

        calls.push(BdCall {
            cell_id,
            locus: field(&fields, locus_idx),
            v: normalize_gene(field_ref(&fields, v_idx)),
            d: normalize_gene(field_ref(&fields, d_idx)),
            j: normalize_gene(field_ref(&fields, j_idx)),
            c: normalize_gene(field_ref(&fields, c_idx)),
            junction: field(&fields, junction_idx),
            junction_aa: field(&fields, junction_aa_idx),
            productive: field(&fields, productive_idx),
            umi_count: field(&fields, umi_count_idx),
        });
    }

    Ok(calls)
}

fn compare_calls(lumrik: &[LumrikCall], bd: &[BdCall]) {
    let mut bd_by_cell_locus: HashMap<(u64, String), Vec<&BdCall>> =
        HashMap::new();

    for call in bd {
        bd_by_cell_locus
            .entry((call.cell_id, call.locus.clone()))
            .or_default()
            .push(call);
    }

    let mut exact_vdj = 0usize;
    let mut exact_vj = 0usize;
    let mut same_cell_locus = 0usize;
    let mut no_cell_locus = 0usize;

    for call in lumrik {
        let key = (call.cell_id, call.chain.clone());

        let Some(candidates) = bd_by_cell_locus.get(&key) else {
            no_cell_locus += 1;
            continue;
        };

        if candidates.iter().any(|bd| {
            bd.v == call.v
                && bd.d == call.d
                && bd.j == call.j
                && bd.c == call.c
        }) {
            exact_vdj += 1;
            continue;
        }

        if candidates.iter().any(|bd| {
            bd.v == call.v && bd.j == call.j
        }) {
            exact_vj += 1;
            continue;
        }

        same_cell_locus += 1;
    }

    println!();
    println!("Lumrik -> BD comparison");
    println!("  Same cell/locus + V/D/J/C: {}", exact_vdj);
    println!("  Same cell/locus + V/J: {}", exact_vj);
    println!("  Same cell/locus only: {}", same_cell_locus);
    println!("  No BD call for cell/locus: {}", no_cell_locus);
}

fn report_lumrik_loci(calls: &[LumrikCall]) {
    let mut counts: HashMap<&str, usize> = HashMap::new();

    for call in calls {
        *counts.entry(call.chain.as_str()).or_default() += 1;
    }

    println!("Lumrik chains");

    let mut values: Vec<_> = counts.into_iter().collect();
    values.sort_by_key(|(chain, _)| *chain);

    for (chain, count) in values {
        println!("  {chain}: {count}");
    }

    println!();
}

fn report_bd_loci(calls: &[BdCall]) {
    let mut counts: HashMap<&str, usize> = HashMap::new();

    for call in calls {
        *counts.entry(call.locus.as_str()).or_default() += 1;
    }

    println!("BD AIRR loci");

    let mut values: Vec<_> = counts.into_iter().collect();
    values.sort_by_key(|(locus, _)| *locus);

    for (locus, count) in values {
        println!("  {locus}: {count}");
    }
}

fn column_index(columns: &[&str], name: &str) -> Result<usize> {
    columns
        .iter()
        .position(|x| *x == name)
        .with_context(|| format!("missing required column {name:?}"))
}

fn field(fields: &[&str], index: usize) -> String {
    field_ref(fields, index).to_string()
}

fn field_ref<'a>(fields: &'a [&str], index: usize) -> &'a str {
    fields.get(index).copied().unwrap_or_default()
}

fn normalize_gene(value: &str) -> String {
    value
        .split(',')
        .next()
        .unwrap_or(value)
        .split('*')
        .next()
        .unwrap_or(value)
        .trim()
        .to_ascii_uppercase()
}
