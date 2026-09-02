use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use mapping_info::MappingInfo;
use sc_vdj::audit::AuditFastaWriter;
use sc_vdj::output::ReportWriter;
use sc_vdj::{
    bam_shard_for_cell, read_bam_filtered, shard_bam_receptor_evidence, Chain, ExpressionMatrix,
    NelruneIdentityResolver, PosteriorAnalyzer, PosteriorConfig, VdjMapper, VdjMapperConfig,
    VdjReferenceBuilder,
};
use std::collections::{HashMap, HashSet};
use std::env;
use std::ffi::OsString;
use std::fs::{self, File};
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
    bam_shards: Option<usize>,
    print_cells: bool,
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
    let mut mapping_info = MappingInfo::new(None, 0.0, 0);

    mapping_info.start_timer("wall.reference_index");
    let mapper = if let Some(index) = &cli.index {
        VdjMapper::load_index(index).with_context(|| format!("loading V(D)J index {}", index.display()))?
    } else {
        eprintln!("building V(D)J germline reference...");
        let reference = VdjReferenceBuilder::default()
            .build(cli.gtf.as_ref().expect("validated --gtf"), cli.genome.as_ref().expect("validated --genome"))
            .with_context(|| "building V(D)J reference")?;
        VdjMapper::new(reference, VdjMapperConfig::default())
    };
    mapping_info.stop_timer("wall.reference_index");
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

    mapping_info.start_timer("wall.expression_matrix");
    let gex = read_mex(&cli.exonic, cli.cell_barcode_len)?;
    mapping_info.stop_timer("wall.expression_matrix");
    if gex.cells.is_empty() {
        bail!("MEX contains no cells");
    }
    println!("Expression matrix: {} cell(s)", gex.cells.len());

    let called_cells: HashSet<String> = gex.cells.iter().cloned().collect();
    let shard_count = cli
        .bam_shards
        .unwrap_or_else(|| automatic_bam_shards(gex.cells.len()));
    println!(
        "BAM processing: {} shard(s) for {} called cell(s)",
        shard_count,
        called_cells.len()
    );
    let print_cells = cli.print_cells || gex.cells.len() <= 64;
    if !print_cells {
        println!(
            "Per-cell stdout suppressed for scale; use --print-cells to show all {} cells",
            gex.cells.len()
        );
    }

    let posterior_config = PosteriorConfig {
        sterile_bins: cli.sterile_bins,
        ..PosteriorConfig::default()
    };
    let analyzer = PosteriorAnalyzer::new(reference, &mapper, posterior_config.clone());
    let mut reports = ReportWriter::create(&cli.out)?;
    let mut audit = AuditFastaWriter::create(&cli.out, reference)?;
    let mut summary_count = 0usize;

    if shard_count == 1 {
        eprintln!(
            "reading BAM evidence for {} called cell(s)...",
            called_cells.len()
        );
        mapping_info.start_timer("wall.bam_input");
        let reads = read_bam_filtered(&cli.bam, &NelruneIdentityResolver, |cell| {
            called_cells.contains(cell)
        })?;
        mapping_info.stop_timer("wall.bam_input");
        mapping_info.total = reads.len();
        println!("BAM evidence: {} records for called MEX cells", reads.len());
        mapping_info.start_timer("wall.posterior");
        let summaries = analyzer.analyze_for_cells_with_mapping_info(
            reads,
            &gex,
            gex.cells.iter().cloned(),
            &mut mapping_info,
        );
        mapping_info.stop_timer("wall.posterior");
        summary_count += summaries.len();
        mapping_info.start_timer("wall.report_write");
        reports.write_cells(&summaries)?;
        audit.write_cells(&summaries)?;
        mapping_info.stop_timer("wall.report_write");
        if print_cells {
            print_summaries(&summaries);
        }
    } else {
        fs::create_dir_all(&cli.out)?;
        let shard_dir = cli
            .out
            .join(format!(".vdj-bam-shards-{}", std::process::id()));
        if shard_dir.exists() {
            fs::remove_dir_all(&shard_dir)
                .with_context(|| format!("removing stale shard directory {}", shard_dir.display()))?;
        }
        fs::create_dir_all(&shard_dir)
            .with_context(|| format!("creating shard directory {}", shard_dir.display()))?;
        let _shard_guard = ShardDirGuard(shard_dir.clone());
        let shard_paths: Vec<_> = (0..shard_count)
            .map(|i| shard_dir.join(format!("receptor-{i:03}.bam")))
            .collect();

        eprintln!(
            "streaming BAM once and retaining only receptor-relevant evidence for {} called cell(s)...",
            called_cells.len()
        );
        mapping_info.start_timer("wall.bam_prefilter_shard_write");
        let stats = shard_bam_receptor_evidence(
            &cli.bam,
            &NelruneIdentityResolver,
            |cell| called_cells.contains(cell),
            &mapper,
            reference,
            posterior_config.min_seed_hits,
            &shard_paths,
        )?;
        mapping_info.stop_timer("wall.bam_prefilter_shard_write");
        mapping_info.total = stats.total_records;
        mapping_info.report_n("vdj.called_cell_records", stats.called_cell_records);
        mapping_info.report_n("vdj.receptor_records_retained", stats.retained_records);
        mapping_info.report_n(
            "vdj.irrelevant_records_discarded",
            stats.discarded_irrelevant_records,
        );
        println!(
            "BAM prefilter: total={} called={} retained={} discarded_irrelevant={}",
            stats.total_records,
            stats.called_cell_records,
            stats.retained_records,
            stats.discarded_irrelevant_records
        );

        let mut cells_by_shard = vec![Vec::<String>::new(); shard_count];
        for cell in &gex.cells {
            cells_by_shard[bam_shard_for_cell(cell, shard_count)].push(cell.clone());
        }

        for shard in 0..shard_count {
            mapping_info.start_timer("wall.bam_shard_read");
            let reads = if stats.shard_records[shard] == 0 {
                Vec::new()
            } else {
                read_bam_filtered(
                    &shard_paths[shard],
                    &NelruneIdentityResolver,
                    |_| true,
                )?
            };
            mapping_info.stop_timer("wall.bam_shard_read");
            eprintln!(
                "VDJ shard {}/{}: {} cell(s), {} retained BAM record(s)",
                shard + 1,
                shard_count,
                cells_by_shard[shard].len(),
                reads.len()
            );
            mapping_info.start_timer("wall.posterior");
            let summaries = analyzer.analyze_for_cells_with_mapping_info(
                reads,
                &gex,
                cells_by_shard[shard].iter().cloned(),
                &mut mapping_info,
            );
            mapping_info.stop_timer("wall.posterior");
            summary_count += summaries.len();
            mapping_info.start_timer("wall.report_write");
            reports.write_cells(&summaries)?;
            audit.write_cells(&summaries)?;
            mapping_info.stop_timer("wall.report_write");
            if print_cells {
                print_summaries(&summaries);
            }
            // `summaries` (including all supporting read sequences and sterile
            // interval state) is dropped here before the next shard is loaded.
        }
    }

    if summary_count == 0 {
        bail!("posterior analyzer produced no cell summaries");
    }
    mapping_info.start_timer("wall.report_finish");
    reports.finish()?;
    audit.finish()?;
    mapping_info.stop_timer("wall.report_finish");

    let performance_path = cli.out.join("vdj-mapping-info.txt");
    fs::write(&performance_path, mapping_info.to_string())
        .with_context(|| format!("writing {}", performance_path.display()))?;
    println!("\n{}", mapping_info);
    println!("VDJ performance report written to {}", performance_path.display());
    println!("Reports written to {}", cli.out.display());
    Ok(())
}

struct ShardDirGuard(PathBuf);

impl Drop for ShardDirGuard {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn automatic_bam_shards(cells: usize) -> usize {
    ((cells + 127) / 128).clamp(1, 64)
}

fn print_summaries(summaries: &[sc_vdj::CellVdjSummary]) {
    println!("\nPosterior V(D)J summaries");
    for cell in summaries {
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
    let mut bam_shards = None;
    let mut print_cells = false;
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
            "--bam-shards" => {
                let n: usize = next_value(&mut args, "--bam-shards")?
                    .to_string_lossy()
                    .parse()
                    .context("invalid --bam-shards")?;
                if n == 0 || n > 256 {
                    bail!("--bam-shards must be between 1 and 256");
                }
                bam_shards = Some(n);
            }
            "--print-cells" => print_cells = true,
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
        bam_shards,
        print_cells,
    })
}

fn next_value<I: Iterator<Item = OsString>>(args: &mut I, flag: &str) -> Result<OsString> {
    args.next()
        .with_context(|| format!("{flag} requires a value"))
}

fn print_help() {
    println!("vdj-summary - run posterior sc-vdj analysis on Nelrune BAM + exonic MEX\n\nUsage:\n  vdj-summary --exonic <MEX_DIR> --bam <BAM> --index <VDJIDX> --out <DIR> [options]\n  vdj-summary --exonic <MEX_DIR> --bam <BAM> --gtf <GTF> --genome <FASTA> --out <DIR> [options]\n\nRequired:\n  --exonic <DIR>     Nelrune exonic MEX (one or more cells)\n  --bam <FILE>       retained/extracted Nelrune BAM (CB/UB tags or mapper QNAME metadata)\n  --index <FILE>     compiled vdj-index (preferred; replaces --gtf/--genome)\n  --gtf <FILE>       matching genome annotation (fallback with --genome)\n  --genome <FASTA>   matching genome FASTA (fallback with --gtf)\n  --out <DIR>        report directory\n\nOptions:\n  --sterile-bins <N>      spatial bins per receptor locus [64]\n  --cell-barcode-len <N>  normalize MEX barcodes to first N bases (legacy BD: 27)\n  --bam-shards <N>        bound whole-BAM memory with N receptor shards [auto]\n                         auto targets about 128 called cells/shard, max 64\n  --print-cells           print every per-cell posterior summary to stdout\n  -h, --help");
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
