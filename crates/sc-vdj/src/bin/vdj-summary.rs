use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use mapping_info::MappingInfo;
use sc_vdj::audit::AuditFastaWriter;
use sc_vdj::output::ReportWriter;
use sc_vdj::{
    read_bam_receptor_evidence, Chain, ExpressionMatrix, NelruneIdentityResolver,
    PosteriorAnalyzer, PosteriorConfig, VdjMapper, VdjMapperConfig, VdjReferenceBuilder,
};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read as IoRead};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    name = "vdj-summary",
    about = "Run detailed posterior single-cell V(D)J analysis on a Nelrune BAM + exonic MEX",
    override_usage = "vdj-summary --exonic <MEX_DIR> --bam <BAM> (--index <VDJIDX> | --gtf <GTF> --genome <FASTA>) --out <DIR> [OPTIONS]",
    after_help = "REFERENCE
  --index is the preferred production path. --gtf + --genome rebuild the
  reference in-memory and are mainly useful for development/debugging.

LEGACY BD
  Older BD MEX barcodes may contain u64 A-padding. Use
  --cell-barcode-len 27 to match the original 27-base BAM cell barcode."
)]
struct Cli {
    /// Nelrune exonic MEX directory.
    #[arg(long, value_name = "MEX_DIR")]
    exonic: PathBuf,

    /// Retained/extracted Nelrune BAM with cell/UMI identity.
    #[arg(long, value_name = "BAM")]
    bam: PathBuf,

    /// Compiled V(D)J index (preferred reference input).
    #[arg(long, value_name = "VDJIDX")]
    index: Option<PathBuf>,

    /// Matching genome annotation; requires --genome when --index is absent.
    #[arg(long, value_name = "GTF")]
    gtf: Option<PathBuf>,

    /// Matching genome FASTA; requires --gtf when --index is absent.
    #[arg(long, value_name = "FASTA")]
    genome: Option<PathBuf>,

    /// Report directory.
    #[arg(long, value_name = "DIR")]
    out: PathBuf,

    /// Spatial bins per receptor locus.
    #[arg(long, default_value_t = 64, value_name = "N")]
    sterile_bins: usize,

    /// Normalize MEX barcodes to the first N bases (legacy BD: 27).
    #[arg(long, value_name = "N")]
    cell_barcode_len: Option<usize>,

    /// Print every per-cell posterior summary to stdout.
    #[arg(long)]
    print_cells: bool,
}

impl Cli {
    fn validate(mut self) -> Result<Self> {
        if self.index.is_some() && (self.gtf.is_some() || self.genome.is_some()) {
            bail!("use either --index <VDJIDX> or --gtf <GTF> --genome <FASTA>, not both");
        }
        if self.index.is_none() && (self.gtf.is_none() || self.genome.is_none()) {
            bail!("missing reference: provide --index <VDJIDX> or both --gtf <GTF> and --genome <FASTA>");
        }
        if let Some(len) = self.cell_barcode_len {
            if !(1..=32).contains(&len) {
                bail!("--cell-barcode-len must be between 1 and 32");
            }
        }
        self.sterile_bins = self.sterile_bins.max(8);
        Ok(self)
    }
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
    let cli = Cli::parse().validate()?;
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
    println!(
        "BAM processing: one-pass in-memory routing for {} called cell(s)",
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
    let analyzer = PosteriorAnalyzer::new(reference, &mapper, posterior_config);
    let mut reports = ReportWriter::create(&cli.out)?;
    let mut audit = AuditFastaWriter::create(&cli.out, reference)?;

    eprintln!(
        "streaming BAM once and routing indexed receptor-segment evidence for {} called cell(s)...",
        called_cells.len()
    );
    mapping_info.start_timer("wall.bam_route_in_memory");
    let (routed, stats) = read_bam_receptor_evidence(
        &cli.bam,
        &NelruneIdentityResolver,
        |cell| called_cells.contains(cell),
        reference,
    )?;
    mapping_info.stop_timer("wall.bam_route_in_memory");
    mapping_info.total = stats.total_records;
    mapping_info.report_n("vdj.receptor_chromosome_records", stats.receptor_chromosome_records);
    mapping_info.report_n("vdj.segment_overlap_records", stats.segment_overlap_records);
    mapping_info.report_n("vdj.called_cell_vdj_records", stats.called_cell_records);
    mapping_info.report_n("vdj.routed_locus_observations", stats.routed_evidence_records);
    for (chain, count) in Chain::ALL.into_iter().zip(stats.locus_records) {
        mapping_info.report_n(&format!("vdj.routed.{}", chain.as_str().to_ascii_lowercase()), count);
    }
    println!(
        "BAM routing: total={} receptor_chr={} segment_overlap={} called_cell={} routed={}",
        stats.total_records,
        stats.receptor_chromosome_records,
        stats.segment_overlap_records,
        stats.called_cell_records,
        stats.routed_evidence_records
    );

    mapping_info.start_timer("wall.posterior");
    let summaries = analyzer.analyze_routed_for_cells_with_mapping_info(
        routed,
        &gex,
        gex.cells.iter().cloned(),
        &mut mapping_info,
    );
    mapping_info.stop_timer("wall.posterior");
    let summary_count = summaries.len();
    mapping_info.start_timer("wall.report_write");
    reports.write_cells(&summaries)?;
    audit.write_cells(&summaries)?;
    mapping_info.stop_timer("wall.report_write");
    if print_cells {
        print_summaries(&summaries);
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

fn print_summaries(summaries: &[sc_vdj::CellVdjSummary]) {
    println!("\nPosterior V(D)J summaries");
    for cell in summaries {
        println!("\nCELL {}", cell.cell);
        println!(
            "  recombination activity: RAG={:.4} RAG1={:.3} RAG2={:.3} DNTT={:.3} paired_RAG={} light_chain={}",
            cell.recombination_activity.rag_activity,
            cell.recombination_activity.rag1_expression,
            cell.recombination_activity.rag2_expression,
            cell.recombination_activity.dntt_expression,
            cell.recombination_activity.rag_pair_detected,
            cell.light_chain_status(),
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
                if let Some(junction) = &r.junction {
                    if r.chain.has_d() {
                        println!(
                            "      junction: V3del={} PV3={} N1={} PD5={} D5del={} Dlen={} D3del={} PD3={} N2={} PJ5={} J5del={} alternative_pn={}",
                            junction.v_del_3,
                            junction.p_v3_len(),
                            junction.n1_len(),
                            junction.p_d5_len(),
                            junction.d_del_5.unwrap_or(0),
                            junction.d_retained_len.unwrap_or(0),
                            junction.d_del_3.unwrap_or(0),
                            junction.p_d3_len(),
                            junction.n2_len(),
                            junction.p_j5_len(),
                            junction.j_del_5,
                            junction.pn_alternative,
                        );
                    } else {
                        println!(
                            "      junction: V3del={} PV3={} N={} PJ5={} J5del={} alternative_pn={}",
                            junction.v_del_3,
                            junction.p_v3_len(),
                            junction.n1_len(),
                            junction.p_j5_len(),
                            junction.j_del_5,
                            junction.pn_alternative,
                        );
                    }
                    println!(
                        "      junction_sequence: PV3={} N1={} PD5={} PD3={} N2={} PJ5={}",
                        String::from_utf8_lossy(&junction.p_v3),
                        String::from_utf8_lossy(&junction.n1),
                        String::from_utf8_lossy(&junction.p_d5),
                        String::from_utf8_lossy(&junction.p_d3),
                        String::from_utf8_lossy(&junction.n2),
                        String::from_utf8_lossy(&junction.p_j5),
                    );
                    println!(
                        "      naive_recombination: {}",
                        String::from_utf8_lossy(&junction.inferred_naive_sequence)
                    );
                }
                if let Some(id) = sc_vdj::PackedRecombinationId::from_call(r) {
                    println!("      recombination_id: {id}");
                }

            }
        }
        if let Some(id) = cell.strongest_heavy_id() {
            println!("  heavy_recombination_id: {id}");
        }
        if let Some(id) = cell.strongest_light_id() {
            println!("  light_recombination_id: {id}");
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
        println!("  recombination GEX evidence:");
        for x in &cell.recombination_activity.contributions {
            if x.expression > 0.0 {
                println!(
                    "    {} expression={:.3} activity={:.4}",
                    x.gene, x.expression, x.activity
                );
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
