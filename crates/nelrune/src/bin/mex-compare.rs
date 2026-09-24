use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};
use int_to_dna::IntToDna;
use sc_primer::{BdCellVersion, RhapsodyWhitelist};
use scdata::{MexFeatureIndex, load_mtx_feature_matrix, read_mtx_barcodes};
use sprs::CsMat;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OtherProducer {
    Bd,
    Tenx,
}

#[derive(Debug, Parser)]
#[command(
    name = "mex-compare",
    about = "Compare a Nelrune Gene Expression MEX matrix with BD Rhapsody or 10x output"
)]
struct Cli {
    /// Nelrune Gene Expression MEX directory. Cell labels must be DNA barcodes.
    #[arg(long)]
    nelrune: PathBuf,

    /// External Gene Expression MEX directory to compare with Nelrune.
    #[arg(long)]
    other: PathBuf,

    /// Producer of the second matrix. BD numeric ids are expanded through the
    /// Rhapsody whitelist; 10x barcode labels are packed directly.
    #[arg(long, value_enum)]
    other_producer: OtherProducer,

    /// BD barcode layout used to expand numeric BD cell ids.
    #[arg(long, default_value = "v2.384")]
    bd_version: String,

    #[arg(long, default_value = "mex_compare")]
    out: PathBuf,

    #[arg(long, default_value_t = 4)]
    threads: usize,
}

#[derive(Debug)]
struct MatrixInput {
    matrix: CsMat<f32>,
    features: MexFeatureIndex,
    canonical_cells: Vec<u64>,
    canonical_labels: Vec<String>,
    labels: Vec<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    fs::create_dir_all(&cli.out)
        .with_context(|| format!("creating {}", cli.out.display()))?;

    let nelrune = load_nelrune(&cli.nelrune, cli.threads)?;
    let other = load_other(
        &cli.other,
        cli.other_producer,
        &cli.bd_version,
        cli.threads,
    )?;

    compare(&nelrune, &other, &cli.out)
}

fn load_nelrune(dir: &Path, threads: usize) -> Result<MatrixInput> {
    let (data, features, _) = load_mtx_feature_matrix(dir, "Gene Expression", threads)
        .with_context(|| format!("loading Nelrune MEX {}", dir.display()))?;
    let barcodes = read_mtx_barcodes(dir)?;
    let canonical_labels = barcodes
        .iter()
        .map(|(label, _)| canonical_dna_label(label))
        .collect::<Result<Vec<_>>>()?;
    let canonical_cells = canonical_labels
        .iter()
        .map(|label| Ok(IntToDna::new(label.as_bytes()).into_u64()))
        .collect::<Result<Vec<_>>>()?;
    ensure_unique_cells(&canonical_cells, &barcodes, "Nelrune")?;

    Ok(MatrixInput {
        matrix: data.as_sprs().map_err(anyhow::Error::msg)?,
        features,
        canonical_cells,
        canonical_labels,
        labels: barcodes.into_iter().map(|(label, _)| label).collect(),
    })
}

fn load_other(
    dir: &Path,
    producer: OtherProducer,
    bd_version: &str,
    threads: usize,
) -> Result<MatrixInput> {
    let (data, features, _) = load_mtx_feature_matrix(dir, "Gene Expression", threads)
        .with_context(|| format!("loading comparison MEX {}", dir.display()))?;
    let barcodes = read_mtx_barcodes(dir)?;

    let canonical_labels = match producer {
        OtherProducer::Tenx => barcodes
            .iter()
            .map(|(label, _)| canonical_dna_label(label))
            .collect::<Result<Vec<_>>>()?,
        OtherProducer::Bd => {
            let version = BdCellVersion::parse(bd_version)
                .map_err(|e| anyhow::anyhow!("invalid BD version {bd_version:?}: {e}"))?;
            let whitelist = RhapsodyWhitelist::builtin(version);
            barcodes
                .iter()
                .map(|(label, _)| {
                    let numeric = label.parse::<u64>().with_context(|| {
                        format!("BD barcode label {label:?} is not a numeric BD cell id")
                    })?;
                    let seq = whitelist.cell_id_to_seq(numeric).with_context(|| {
                        format!(
                            "BD cell id {numeric} is outside the {:?} Rhapsody whitelist",
                            version
                        )
                    })?;
                    String::from_utf8(seq).context("BD whitelist returned a non-UTF8 barcode")
                })
                .collect::<Result<Vec<_>>>()?
        }
    };
    let canonical_cells = canonical_labels
        .iter()
        .map(|label| Ok(IntToDna::new(label.as_bytes()).into_u64()))
        .collect::<Result<Vec<_>>>()?;
    ensure_unique_cells(&canonical_cells, &barcodes, "comparison")?;

    Ok(MatrixInput {
        matrix: data.as_sprs().map_err(anyhow::Error::msg)?,
        features,
        canonical_cells,
        canonical_labels,
        labels: barcodes.into_iter().map(|(label, _)| label).collect(),
    })
}

fn canonical_dna_label(label: &str) -> Result<String> {
    let seq = label.split_once('-').map(|(seq, _)| seq).unwrap_or(label);
    if seq.is_empty()
        || !seq
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    {
        bail!("cell label {label:?} is not a DNA barcode");
    }
    Ok(seq.to_ascii_uppercase())
}

fn ensure_unique_cells(
    ids: &[u64],
    barcodes: &[(String, u64)],
    source: &str,
) -> Result<()> {
    let mut seen = HashMap::<u64, &str>::new();
    for (id, (label, _)) in ids.iter().zip(barcodes) {
        if let Some(previous) = seen.insert(*id, label) {
            bail!(
                "{source} cell labels {previous:?} and {label:?} collapse to the same canonical cell id"
            );
        }
    }
    Ok(())
}

fn compare(a: &MatrixInput, b: &MatrixInput, out: &Path) -> Result<()> {
    let a_cells: HashMap<u64, usize> = a
        .canonical_cells
        .iter()
        .copied()
        .enumerate()
        .map(|(col, id)| (id, col))
        .collect();
    let b_cells: HashMap<u64, usize> = b
        .canonical_cells
        .iter()
        .copied()
        .enumerate()
        .map(|(col, id)| (id, col))
        .collect();

    let mut shared_cells: Vec<(u64, usize, usize)> = a_cells
        .iter()
        .filter_map(|(id, a_col)| b_cells.get(id).map(|b_col| (*id, *a_col, *b_col)))
        .collect();
    shared_cells.sort_unstable_by_key(|(_, a_col, _)| *a_col);

    let a_gene_map = unique_gene_names(&a.features);
    let b_gene_map = unique_gene_names(&b.features);
    let mut shared_genes: Vec<(String, usize, usize)> = a_gene_map
        .iter()
        .filter_map(|(name, a_row)| {
            b_gene_map
                .get(name)
                .map(|b_row| (name.clone(), *a_row, *b_row))
        })
        .collect();
    shared_genes.sort_unstable_by_key(|(_, a_row, _)| *a_row);

    let a_t = a.matrix.transpose_view().to_csr();
    let b_t = b.matrix.transpose_view().to_csr();

    write_cell_table(out, a, b, &a_t, &b_t, &shared_cells, &shared_genes)?;
    write_other_only_cells(out, a, b, &b_t)?;
    let gene_corr = write_gene_table(out, a, b, &shared_cells, &shared_genes)?;

    let a_total: f64 = a.matrix.data().iter().map(|x| f64::from(*x)).sum();
    let b_total: f64 = b.matrix.data().iter().map(|x| f64::from(*x)).sum();
    let a_set: HashSet<_> = a.canonical_cells.iter().copied().collect();
    let b_set: HashSet<_> = b.canonical_cells.iter().copied().collect();

    let mut summary = BufWriter::new(File::create(out.join("summary.txt"))?);
    writeln!(summary, "Nelrune vs external MEX comparison")?;
    writeln!(summary, "==================================")?;
    writeln!(summary, "Nelrune features: {}", a.matrix.rows())?;
    writeln!(summary, "Other features: {}", b.matrix.rows())?;
    writeln!(summary, "Shared gene names: {}", shared_genes.len())?;
    writeln!(summary, "Nelrune cells: {}", a.matrix.cols())?;
    writeln!(summary, "Other cells: {}", b.matrix.cols())?;
    writeln!(summary, "Shared canonical cells: {}", shared_cells.len())?;
    writeln!(summary, "Nelrune-only cells: {}", a_set.difference(&b_set).count())?;
    writeln!(summary, "Other-only cells: {}", b_set.difference(&a_set).count())?;
    writeln!(summary, "Nelrune total counts: {:.0}", a_total)?;
    writeln!(summary, "Other total counts: {:.0}", b_total)?;
    writeln!(
        summary,
        "Shared-gene total-count Pearson (shared cells): {}",
        format_corr(gene_corr)
    )?;

    println!("Nelrune cells: {}", a.matrix.cols());
    println!("Other cells: {}", b.matrix.cols());
    println!("Shared cells: {}", shared_cells.len());
    println!("Shared genes: {}", shared_genes.len());
    println!("Gene-total Pearson: {}", format_corr(gene_corr));
    println!("Wrote {}", out.display());
    Ok(())
}

fn unique_gene_names(index: &MexFeatureIndex) -> HashMap<String, usize> {
    let mut counts = HashMap::<&str, usize>::new();
    for feature in index.features() {
        *counts.entry(feature.name.as_str()).or_default() += 1;
    }
    index
        .features()
        .iter()
        .enumerate()
        .filter(|(_, feature)| counts.get(feature.name.as_str()) == Some(&1))
        .map(|(row, feature)| (feature.name.clone(), row))
        .collect()
}

fn write_cell_table(
    out: &Path,
    a: &MatrixInput,
    b: &MatrixInput,
    a_t: &CsMat<f32>,
    b_t: &CsMat<f32>,
    shared_cells: &[(u64, usize, usize)],
    shared_genes: &[(String, usize, usize)],
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(out.join("cells.tsv"))?);
    writeln!(
        writer,
        "cell_id\tnelrune_label\tother_label\tnelrune_counts\tother_counts\tnelrune_genes\tother_genes\tshared_gene_pearson"
    )?;

    let b_to_a: HashMap<usize, usize> = shared_genes
        .iter()
        .map(|(_, a_row, b_row)| (*b_row, *a_row))
        .collect();
    let shared_a: HashSet<usize> = shared_genes.iter().map(|(_, a_row, _)| *a_row).collect();

    for (cell_id, a_col, b_col) in shared_cells {
        let av = a_t.outer_view(*a_col).context("missing Nelrune cell column")?;
        let bv = b_t.outer_view(*b_col).context("missing comparison cell column")?;
        let a_counts: f64 = av.data().iter().map(|x| f64::from(*x)).sum();
        let b_counts: f64 = bv.data().iter().map(|x| f64::from(*x)).sum();

        let mut ax = HashMap::<usize, f64>::new();
        for (row, value) in av.iter() {
            if shared_a.contains(&row) {
                ax.insert(row, f64::from(*value));
            }
        }
        let mut pairs = Vec::new();
        for (b_row, value) in bv.iter() {
            if let Some(a_row) = b_to_a.get(&b_row) {
                pairs.push((*a_row, f64::from(*value)));
            }
        }
        let corr = sparse_pearson(&ax, &pairs, shared_genes.len());

        writeln!(
            writer,
            "{cell_id}\t{}\t{}\t{a_counts:.0}\t{b_counts:.0}\t{}\t{}\t{}",
            a.labels[*a_col],
            b.labels[*b_col],
            av.nnz(),
            bv.nnz(),
            format_corr(corr)
        )?;
    }
    Ok(())
}


fn write_other_only_cells(
    out: &Path,
    a: &MatrixInput,
    b: &MatrixInput,
    b_t: &CsMat<f32>,
) -> Result<()> {
    let a_cells: HashSet<u64> = a.canonical_cells.iter().copied().collect();
    let mut rows = Vec::<(f64, usize, u64)>::new();

    for (b_col, cell_id) in b.canonical_cells.iter().copied().enumerate() {
        if a_cells.contains(&cell_id) {
            continue;
        }
        let bv = b_t.outer_view(b_col).context("missing comparison cell column")?;
        let counts: f64 = bv.data().iter().map(|x| f64::from(*x)).sum();
        rows.push((counts, b_col, cell_id));
    }

    rows.sort_unstable_by(|left, right| {
        right
            .0
            .total_cmp(&left.0)
            .then_with(|| left.1.cmp(&right.1))
    });

    let mut writer = BufWriter::new(File::create(out.join("other_only_cells.tsv"))?);
    writeln!(
        writer,
        "canonical_cell_id\tcanonical_barcode\tother_label\tother_counts\tother_genes"
    )?;

    for (counts, b_col, cell_id) in rows {
        let bv = b_t.outer_view(b_col).context("missing comparison cell column")?;
        writeln!(
            writer,
            "{cell_id}\t{}\t{}\t{counts:.0}\t{}",
            b.canonical_labels[b_col],
            b.labels[b_col],
            bv.nnz()
        )?;
    }

    Ok(())
}

fn write_gene_table(
    out: &Path,
    a: &MatrixInput,
    b: &MatrixInput,
    shared_cells: &[(u64, usize, usize)],
    shared_genes: &[(String, usize, usize)],
) -> Result<Option<f64>> {
    let mut a_totals = vec![0.0_f64; shared_genes.len()];
    let mut b_totals = vec![0.0_f64; shared_genes.len()];
    let a_lookup: HashMap<usize, usize> = shared_genes
        .iter()
        .enumerate()
        .map(|(i, (_, row, _))| (*row, i))
        .collect();
    let b_lookup: HashMap<usize, usize> = shared_genes
        .iter()
        .enumerate()
        .map(|(i, (_, _, row))| (*row, i))
        .collect();
    let shared_a_cols: HashSet<usize> = shared_cells.iter().map(|(_, col, _)| *col).collect();
    let shared_b_cols: HashSet<usize> = shared_cells.iter().map(|(_, _, col)| *col).collect();

    for (row, vec) in a.matrix.outer_iterator().enumerate() {
        let Some(i) = a_lookup.get(&row) else { continue };
        a_totals[*i] = vec
            .iter()
            .filter(|(col, _)| shared_a_cols.contains(col))
            .map(|(_, value)| f64::from(*value))
            .sum();
    }
    for (row, vec) in b.matrix.outer_iterator().enumerate() {
        let Some(i) = b_lookup.get(&row) else { continue };
        b_totals[*i] = vec
            .iter()
            .filter(|(col, _)| shared_b_cols.contains(col))
            .map(|(_, value)| f64::from(*value))
            .sum();
    }

    let corr = dense_pearson(&a_totals, &b_totals);
    let mut writer = BufWriter::new(File::create(out.join("genes.tsv"))?);
    writeln!(writer, "gene\tnelrune_counts\tother_counts\tlog2_ratio_nelrune_over_other")?;
    for (i, (name, _, _)) in shared_genes.iter().enumerate() {
        let ratio = ((a_totals[i] + 1.0) / (b_totals[i] + 1.0)).log2();
        writeln!(writer, "{name}\t{:.0}\t{:.0}\t{ratio:.6}", a_totals[i], b_totals[i])?;
    }
    Ok(corr)
}

fn sparse_pearson(
    a: &HashMap<usize, f64>,
    b: &[(usize, f64)],
    n: usize,
) -> Option<f64> {
    if n < 2 {
        return None;
    }
    let sum_a: f64 = a.values().sum();
    let sumsq_a: f64 = a.values().map(|x| x * x).sum();
    let sum_b: f64 = b.iter().map(|(_, x)| *x).sum();
    let sumsq_b: f64 = b.iter().map(|(_, x)| x * x).sum();
    let dot: f64 = b
        .iter()
        .filter_map(|(row, y)| a.get(row).map(|x| x * y))
        .sum();
    pearson_from_sums(n, sum_a, sum_b, sumsq_a, sumsq_b, dot)
}

fn dense_pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.len() < 2 {
        return None;
    }
    let sum_a: f64 = a.iter().sum();
    let sum_b: f64 = b.iter().sum();
    let sumsq_a: f64 = a.iter().map(|x| x * x).sum();
    let sumsq_b: f64 = b.iter().map(|x| x * x).sum();
    let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
    pearson_from_sums(a.len(), sum_a, sum_b, sumsq_a, sumsq_b, dot)
}

fn pearson_from_sums(
    n: usize,
    sum_x: f64,
    sum_y: f64,
    sumsq_x: f64,
    sumsq_y: f64,
    dot: f64,
) -> Option<f64> {
    let n = n as f64;
    let numerator = dot - sum_x * sum_y / n;
    let vx = sumsq_x - sum_x * sum_x / n;
    let vy = sumsq_y - sum_y * sum_y / n;
    if vx <= 0.0 || vy <= 0.0 {
        None
    } else {
        Some(numerator / (vx * vy).sqrt())
    }
}

fn format_corr(value: Option<f64>) -> String {
    value
        .map(|x| format!("{x:.6}"))
        .unwrap_or_else(|| "NA".to_string())
}
