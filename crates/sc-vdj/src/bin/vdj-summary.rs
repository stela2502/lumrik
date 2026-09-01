use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use sc_vdj::audit::write_vdj_audit_fastas;
use sc_vdj::output::write_reports;
use sc_vdj::{
    read_bam_filtered, Chain, ExpressionMatrix, NelruneIdentityResolver, PosteriorAnalyzer,
    PosteriorConfig, VdjMapper, VdjMapperConfig, VdjReferenceBuilder,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufRead, BufReader, Read as IoRead};
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct Cli {
    exonic: PathBuf,
    bam: PathBuf,
    index: Option<PathBuf>,
    gtf: Option<PathBuf>,
    genome: Option<PathBuf>,
    out: PathBuf,
    sterile_bins: usize,
    cell_barcode_len: Option<usize>,
}

#[derive(Debug, Default)]
struct MexExpression {
    cells: Vec<String>,
    values: HashMap<(String, String), f64>,
}

impl ExpressionMatrix for MexExpression {
    fn expression(&self, cell: &str, gene: &str) -> f64 {
        *self
            .values
            .get(&(cell.to_string(), gene.to_ascii_uppercase()))
            .unwrap_or(&0.0)
    }
}

fn main() -> Result<()> {
    let cli = parse_cli(env::args_os().skip(1))?;

    let mapper = if let Some(index) = &cli.index {
        VdjMapper::load_index(index).with_context(|| format!("loading V(D)J index {}", index.display()))?
    } else {
        eprintln!("building V(D)J germline reference...");
        let reference = VdjReferenceBuilder::default()
            .build(cli.gtf.as_ref().expect("validated --gtf"), cli.genome.as_ref().expect("validated --genome"))
            .with_context(|| "building V(D)J reference")?;
        VdjMapper::new(reference, VdjMapperConfig::default())
    };
    let reference = mapper.reference();

    println!("Reference segments: {}", reference.len());
    for chain in Chain::ALL {
        let mut parts = Vec::new();
        for kind in [
            sc_vdj::SegmentKind::V,
            sc_vdj::SegmentKind::D,
            sc_vdj::SegmentKind::J,
            sc_vdj::SegmentKind::C,
        ] {
            let n = reference.segments_for(chain, kind).count();
            if n > 0 {
                parts.push(format!("{:?}={n}", kind));
            }
        }
        if !parts.is_empty() {
            println!("  {chain}: {}", parts.join(" "));
        }
    }

    let gex = read_mex(&cli.exonic, cli.cell_barcode_len)?;
    if gex.cells.is_empty() {
        bail!("MEX contains no cells");
    }
    println!("Expression matrix: {} cell(s)", gex.cells.len());

    let called_cells: HashSet<String> = gex.cells.iter().cloned().collect();
    eprintln!(
        "reading BAM evidence for {} called cell(s)...",
        called_cells.len()
    );
    let reads = read_bam_filtered(&cli.bam, &NelruneIdentityResolver, |cell| {
        called_cells.contains(cell)
    })?;
    println!("BAM evidence: {} records for called MEX cells", reads.len());

    let analyzer = PosteriorAnalyzer::new(
        reference,
        &mapper,
        PosteriorConfig {
            sterile_bins: cli.sterile_bins,
            ..PosteriorConfig::default()
        },
    );
    let summaries = analyzer.analyze(reads, &gex);
    if summaries.is_empty() {
        bail!("posterior analyzer produced no cell summaries");
    }

    write_reports(&cli.out, &summaries)?;
    write_vdj_audit_fastas(&cli.out, &summaries, reference)?;

    println!("\nPosterior V(D)J summaries");
    for cell in &summaries {
        println!("\nCELL {}", cell.cell);
        println!(
            "  development: {:?} score={:.4} RAG={:.4} evidence_code={}",
            cell.development.program,
            cell.development.probability_like_score,
            cell.development.rag_activity,
            cell.development.evidence_code
        );

        if cell.rearrangements.is_empty() {
            println!("  rearrangements: none called");
        } else {
            println!("  rearrangements:");
            for r in &cell.rearrangements {
                println!(
                    "    {} {:?}: {}  support={} UMI(s)",
                    r.chain, r.stage, r.notation, r.total_supporting_umis
                );
                print_segment("V", r.v.as_ref());
                print_segment("D", r.d.as_ref());
                print_segment("J", r.j.as_ref());
                print_segment("C", r.c.as_ref());
            }
        }

        println!("  sterile/germline locus evidence:");
        for p in &cell.sterile {
            if p.total_unique_umis == 0 {
                continue;
            }
            println!(
                "    {}: {} UMI(s), {} aligned block(s), breadth={:.4}, centroid={:.4}, proximal={:.4}, distal={:.4}",
                p.chain,
                p.total_unique_umis,
                p.supported_intervals.len(),
                p.breadth,
                p.centroid,
                p.proximal_fraction,
                p.distal_fraction
            );
        }

        if cell.development.contributions.is_empty() {
            println!("  GEX rationale: no configured developmental markers expressed");
        } else {
            println!("  GEX rationale:");
            for x in &cell.development.contributions {
                if x.expression > 0.0 {
                    println!(
                        "    {} expression={:.3} weight={:.2} contribution={:.4}",
                        x.gene, x.expression, x.weight, x.contribution
                    );
                }
            }
        }
    }

    println!("\nReports written to {}", cli.out.display());
    Ok(())
}

fn print_segment(label: &str, x: Option<&sc_vdj::GermlineSegmentSupport>) {
    if let Some(x) = x {
        println!(
            "      {label}: {} score={} reads={} UMIs={} locus={:.4} distance={} bp",
            x.id,
            x.local_alignment_score,
            x.supporting_reads,
            x.supporting_umis,
            x.locus_fraction,
            x.distance_to_recombination_center
        );
    }
}

fn parse_cli<I: IntoIterator<Item = OsString>>(args: I) -> Result<Cli> {
    let mut exonic = None;
    let mut bam = None;
    let mut index = None;
    let mut gtf = None;
    let mut genome = None;
    let mut out = None;
    let mut sterile_bins = 64usize;
    let mut cell_barcode_len = None;
    let mut args = args.into_iter();
    while let Some(arg) = args.next() {
        match arg.to_string_lossy().as_ref() {
            "-h" | "--help" => {
                print_help();
                std::process::exit(0);
            }
            "--exonic" => exonic = Some(PathBuf::from(next_value(&mut args, "--exonic")?)),
            "--bam" => bam = Some(PathBuf::from(next_value(&mut args, "--bam")?)),
            "--index" => index = Some(PathBuf::from(next_value(&mut args, "--index")?)),
            "--gtf" => gtf = Some(PathBuf::from(next_value(&mut args, "--gtf")?)),
            "--genome" => genome = Some(PathBuf::from(next_value(&mut args, "--genome")?)),
            "--out" => out = Some(PathBuf::from(next_value(&mut args, "--out")?)),
            "--sterile-bins" => {
                sterile_bins = next_value(&mut args, "--sterile-bins")?
                    .to_string_lossy()
                    .parse()
                    .context("invalid --sterile-bins")?
            }
            "--cell-barcode-len" => {
                let len: usize = next_value(&mut args, "--cell-barcode-len")?
                    .to_string_lossy()
                    .parse()
                    .context("invalid --cell-barcode-len")?;
                if len == 0 || len > 32 {
                    bail!("--cell-barcode-len must be between 1 and 32");
                }
                cell_barcode_len = Some(len);
            }
            other => bail!("unknown argument {other}\n\nRun vdj-summary --help for usage."),
        }
    }
    if index.is_some() && (gtf.is_some() || genome.is_some()) {
        bail!("use either --index <VDJIDX> or --gtf <GTF> --genome <FASTA>, not both");
    }
    if index.is_none() && (gtf.is_none() || genome.is_none()) {
        bail!("missing reference: provide --index <VDJIDX> or both --gtf <GTF> and --genome <FASTA>");
    }
    Ok(Cli {
        exonic: exonic.context("missing --exonic <MEX_DIR>")?,
        bam: bam.context("missing --bam <BAM>")?,
        index,
        gtf,
        genome,
        out: out.context("missing --out <DIR>")?,
        sterile_bins: sterile_bins.max(8),
        cell_barcode_len,
    })
}

fn next_value<I: Iterator<Item = OsString>>(args: &mut I, flag: &str) -> Result<OsString> {
    args.next()
        .with_context(|| format!("{flag} requires a value"))
}

fn print_help() {
    println!("vdj-summary - run posterior sc-vdj analysis on Nelrune BAM + exonic MEX\n\nUsage:\n  vdj-summary --exonic <MEX_DIR> --bam <BAM> --index <VDJIDX> --out <DIR> [options]\n  vdj-summary --exonic <MEX_DIR> --bam <BAM> --gtf <GTF> --genome <FASTA> --out <DIR> [options]\n\nRequired:\n  --exonic <DIR>     Nelrune exonic MEX (one or more cells)\n  --bam <FILE>       retained/extracted Nelrune BAM (CB/UB tags or mapper QNAME metadata)\n  --index <FILE>     compiled vdj-index (preferred; replaces --gtf/--genome)\n  --gtf <FILE>       matching genome annotation (fallback with --genome)\n  --genome <FASTA>   matching genome FASTA (fallback with --gtf)\n  --out <DIR>        report directory\n\nOptions:\n  --sterile-bins <N>      spatial bins per receptor locus [64]\n  --cell-barcode-len <N>   normalize MEX barcodes to first N bases (legacy BD: 27)\n  -h, --help");
}

fn find_file(dir: &Path, names: &[&str]) -> Result<PathBuf> {
    for name in names {
        let p = dir.join(name);
        if p.is_file() {
            return Ok(p);
        }
    }
    bail!("none of {} found in {}", names.join(", "), dir.display())
}

fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let input: Box<dyn IoRead> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(input)))
}

fn read_mex(dir: &Path, cell_barcode_len: Option<usize>) -> Result<MexExpression> {
    let barcode_path = find_file(dir, &["barcodes.tsv.gz", "barcodes.tsv"])?;
    let feature_path = find_file(
        dir,
        &[
            "features.tsv.gz",
            "features.tsv",
            "genes.tsv.gz",
            "genes.tsv",
        ],
    )?;
    let matrix_path = find_file(dir, &["matrix.mtx.gz", "matrix.mtx"])?;

    let raw_cells: Vec<String> = reader(&barcode_path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .filter_map(|line| line.split('\t').next().map(str::to_string))
        .filter(|x| !x.is_empty())
        .collect();

    let mut cells = Vec::with_capacity(raw_cells.len());
    let mut seen = HashSet::with_capacity(raw_cells.len());
    for raw in raw_cells {
        let cell = match cell_barcode_len {
            Some(len) => raw
                .get(..len)
                .with_context(|| {
                    format!("matrix barcode {raw} is shorter than --cell-barcode-len {len}")
                })?
                .to_string(),
            None => raw,
        };
        if !seen.insert(cell.clone()) {
            bail!("cell barcode normalization produced duplicate barcode {cell}; refusing ambiguous MEX");
        }
        cells.push(cell);
    }

    let genes: Vec<String> = reader(&feature_path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            fields
                .get(1)
                .copied()
                .unwrap_or(fields[0])
                .to_ascii_uppercase()
        })
        .collect();

    let mut input = reader(&matrix_path)?;
    let mut line = String::new();
    loop {
        line.clear();
        if input.read_line(&mut line)? == 0 {
            bail!("matrix ended before dimensions");
        }
        if !line.starts_with('%') {
            break;
        }
    }
    let dims: Vec<&str> = line.split_whitespace().collect();
    if dims.len() != 3 {
        bail!("invalid MatrixMarket dimensions line");
    }
    let rows: usize = dims[0].parse()?;
    let cols: usize = dims[1].parse()?;
    if rows != genes.len() || cols != cells.len() {
        bail!(
            "MEX dimensions {}x{} do not match {} features and {} cells",
            rows,
            cols,
            genes.len(),
            cells.len()
        );
    }

    let mut values = HashMap::new();
    for line in input.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut f = line.split_whitespace();
        let row: usize = f.next().context("matrix row")?.parse()?;
        let col: usize = f.next().context("matrix col")?.parse()?;
        let value: f64 = f.next().context("matrix value")?.parse()?;
        if row == 0 || col == 0 || row > genes.len() || col > cells.len() {
            bail!("matrix entry outside declared dimensions: row={row}, col={col}");
        }
        values.insert((cells[col - 1].clone(), genes[row - 1].clone()), value);
    }
    Ok(MexExpression { cells, values })
}
