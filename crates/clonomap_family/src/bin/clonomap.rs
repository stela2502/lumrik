use anyhow::{Context, Result, bail};
use clap::Parser;
use clonomap::{
    AlignedCell, CellReceptor, ClonoMap, Family, FamilyConfig, LightClone, LightMember,
    MutationMeasurement, Receptor, ReferenceModels, align_fragment, project_observed_bases,
    reconstruct_substitutions, rooted_categorical_hex, rooted_continuous_hex,
    translate_imgt_alignment,
};
use lumrik_status::{
    ServerContent, ServerSnapshot, StatusMetric, StatusSection, memory_status, public_hostname,
    spawn_status_server,
};
use ndarray::Array2;
use ommverse::AaModelVault;
use rayon::prelude::*;
use sc_primer::{Chemistry, PrimerDetector};
use sc_vdj::RecombinationId;
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::net::{IpAddr, Ipv4Addr, SocketAddr};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::{Arc, RwLock};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Parser, Debug)]
#[command(
    name = "clonomap",
    author,
    version,
    about = "Mutation-aware receptor-family analysis for Lumrik"
)]
struct Cli {
    /// Lumrik/nelrune-vdj output directory containing airr_rearrangements.tsv and vdj_calls.tsv
    #[arg(long, num_args = 1..)]
    vdj_out: Vec<PathBuf>,
    /// Output directory for family, mutation, model-space and plot reports
    #[arg(long)]
    out: PathBuf,
    /// Single-cell chemistry/chemistries used to identify and clip known technical sequence.
    #[arg(long, value_enum, num_args = 1.., default_values_t = [Chemistry::BdV2_384, Chemistry::BdV2_384Vdj])]
    chemistry: Vec<Chemistry>,
    /// Maximum nucleotide CDR3 edit distance during initial HC family collection
    #[arg(long, default_value_t = 3)]
    max_cdr3_distance: usize,
    /// Minimum final HC family size reported in the console selection summary and plotted regardless of Pearson significance.
    #[arg(long, default_value_t = 100)]
    min_family_size: usize,
    /// Maximum nominal Pearson p-value reported in the console selection summary and plotted regardless of family size.
    #[arg(long, default_value_t = 0.05)]
    max_pearson_p: f64,
    /// Minimum identity for the first HC mutation-alignment pass
    #[arg(long, default_value_t = 0.50)]
    min_alignment_identity: f64,
    /// Hard identity required when an ejected HC tries another family
    #[arg(long, default_value_t = 0.75)]
    hard_alignment_identity: f64,
    /// Persistent reference_curator store. Discovered sequences keep stable identities across runs.
    #[arg(long)]
    reference_curator: Option<PathBuf>,
    /// Generic Ommverse amino-acid model vault. When supplied, add speculative family model-space trajectories.
    #[arg(long)]
    vault: Option<PathBuf>,
    /// Skip the default family SVG plots.
    #[arg(long)]
    no_plots: bool,
    /// Draw rooted family plots radially instead of the default layered layout.
    #[arg(long)]
    radial_layout: bool,
    /// Starting port for the ClonoMap status server; the first free port at or above it is used.
    #[arg(long, default_value_t = 8787)]
    health_port: u16,
    /// Hostname advertised by the status server. Defaults to the machine hostname.
    #[arg(long)]
    health_hostname: Option<String>,
}

struct ClonoMapStatus {
    started_unix_ms: u128,
    finished_unix_ms: Option<u128>,
    stage: String,
    detail: String,
    public_url: Option<String>,
    output: String,
    input_sources: usize,
    receptor_calls: usize,
    hc_cells: usize,
    families: usize,
    lc_clones: usize,
    hc_only_attached: usize,
    hc_only_ambiguous: usize,
    hc_only_unresolved: usize,
    lc_only_attached: usize,
    lc_only_ambiguous: usize,
    lc_only_unresolved: usize,
    top_lc_rescues: Vec<String>,
    vault: Option<String>,
    vault_models: Option<usize>,
}

impl ClonoMapStatus {
    fn new(args: &Cli) -> Self {
        Self {
            started_unix_ms: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
            finished_unix_ms: None,
            stage: "starting".to_owned(),
            detail: "Starting ClonoMap".to_owned(),
            public_url: None,
            output: args.out.display().to_string(),
            input_sources: args.vdj_out.len(),
            receptor_calls: 0,
            hc_cells: 0,
            families: 0,
            lc_clones: 0,
            hc_only_attached: 0,
            hc_only_ambiguous: 0,
            hc_only_unresolved: 0,
            lc_only_attached: 0,
            lc_only_ambiguous: 0,
            lc_only_unresolved: 0,
            top_lc_rescues: Vec::new(),
            vault: args.vault.as_ref().map(|p| p.display().to_string()),
            vault_models: None,
        }
    }

    fn stage(&mut self, stage: &str, detail: impl Into<String>) {
        self.stage = stage.to_owned();
        self.detail = detail.into();
    }

    fn finish(&mut self) {
        self.stage("complete", "ClonoMap analysis complete");
        self.finished_unix_ms = Some(
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
        );
    }
}

impl ServerContent for ClonoMapStatus {
    fn server_snapshot(&self) -> ServerSnapshot {
        let mut run = vec![
            StatusMetric::new("Detail", self.detail.clone()),
            StatusMetric::new("Input VDJ directories", self.input_sources.to_string()),
            StatusMetric::new("Receptor calls", self.receptor_calls.to_string()),
            StatusMetric::new("Qualified/final HC cells", self.hc_cells.to_string()),
            StatusMetric::new("HC families", self.families.to_string()),
            StatusMetric::new("LC clones", self.lc_clones.to_string()),
            StatusMetric::new(
                "HC-only attached to clone",
                self.hc_only_attached.to_string(),
            ),
            StatusMetric::new("HC-only ambiguous", self.hc_only_ambiguous.to_string()),
            StatusMetric::new("HC-only unresolved", self.hc_only_unresolved.to_string()),
            StatusMetric::new(
                "LC-only attached to clone",
                self.lc_only_attached.to_string(),
            ),
            StatusMetric::new("LC-only ambiguous", self.lc_only_ambiguous.to_string()),
            StatusMetric::new("LC-only unresolved", self.lc_only_unresolved.to_string()),
        ];
        for (rank, summary) in self.top_lc_rescues.iter().enumerate() {
            run.push(StatusMetric::new(
                format!("Top rescued LC #{}", rank + 1),
                summary.clone(),
            ));
        }
        if let Some(vault) = &self.vault {
            run.push(StatusMetric::new("Model vault", vault.clone()));
            run.push(StatusMetric::new(
                "Vault models",
                self.vault_models
                    .map(|n| n.to_string())
                    .unwrap_or_else(|| "not loaded".to_owned()),
            ));
        }
        ServerSnapshot {
            title: "ClonoMap".to_owned(),
            subtitle: "Mutation-aware receptor-family analysis".to_owned(),
            started_unix_ms: self.started_unix_ms,
            finished_unix_ms: self.finished_unix_ms,
            stage: self.stage.clone(),
            public_url: self.public_url.clone(),
            sections: vec![
                StatusSection::new("Current analysis", run),
                StatusSection::new(
                    "Output",
                    vec![StatusMetric::new("Directory", self.output.clone())],
                ),
                {
                    let mem = memory_status();
                    StatusSection::new(
                        "Process",
                        vec![
                            StatusMetric::new("RSS", format!("{:.1} MiB", mem.process_rss_mib)),
                            StatusMetric::new(
                                "Peak RSS",
                                format!("{:.1} MiB", mem.process_peak_rss_mib),
                            ),
                            StatusMetric::new(
                                "System available",
                                format!("{:.1} MiB", mem.system_available_mib),
                            ),
                        ],
                    )
                },
            ],
        }
    }
}

fn set_status(state: &Arc<RwLock<ClonoMapStatus>>, stage: &str, detail: impl Into<String>) {
    if let Ok(mut status) = state.write() {
        status.stage(stage, detail);
    }
}

#[derive(Clone)]
struct CallRow {
    cell: String,
    source: String,
    receptor: Receptor,
    productive: bool,
}

#[derive(Debug)]
struct UnassignedCell {
    cell: CellReceptor,
    from_family: String,
}

#[derive(Debug, Clone)]
struct LightOnlyCell {
    cell: String,
    evidence_reads: usize,
    support_umis: usize,
}

#[derive(Debug, Clone)]
struct LightOnlyGroup {
    receptor: Receptor,
    cells: Vec<LightOnlyCell>,
}

#[derive(Debug, Clone)]
struct LightOnlyAssignment {
    family: String,
    clone_id: String,
    receptor: Receptor,
    cells: Vec<LightOnlyCell>,
}

fn mean_sd(values: &[usize]) -> (f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0);
    }
    let n = values.len() as f64;
    let mean = values.iter().map(|&v| v as f64).sum::<f64>() / n;
    let variance = values
        .iter()
        .map(|&v| {
            let d = v as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n;
    (mean, variance.sqrt())
}

fn lc_rescue_summary(assignment: &LightOnlyAssignment, total: usize) -> String {
    let reads: Vec<usize> = assignment.cells.iter().map(|c| c.evidence_reads).collect();
    let umis: Vec<usize> = assignment.cells.iter().map(|c| c.support_umis).collect();
    let (read_mean, read_sd) = mean_sd(&reads);
    let (umi_mean, umi_sd) = mean_sd(&umis);
    let share = if total == 0 {
        0.0
    } else {
        100.0 * assignment.cells.len() as f64 / total as f64
    };
    format!(
        "{} cells ({share:.1}%), evidence reads {read_mean:.2} +/- {read_sd:.2}, UMIs {umi_mean:.2} +/- {umi_sd:.2}; {}",
        assignment.cells.len(),
        assignment.clone_id
    )
}

/// Deliberately lossy, mutation-agnostic chain identity used by ClonoMap clones.
/// Layout is five fixed-width u16 hex fields: V segment id, V 3' clipping,
/// reconstructed V->J intervening length, J 5' clipping, J segment id.
fn clone_chain_hex(receptor: &Receptor) -> Result<String> {
    let id = RecombinationId::from_str(&receptor.id)
        .map_err(anyhow::Error::msg)
        .with_context(|| format!("decode recombination id {}", receptor.id))?;
    let decoded = id
        .decode_numeric_candidates()
        .map_err(anyhow::Error::msg)
        .with_context(|| format!("decode recombination id {}", receptor.id))?;
    if decoded.len() != 1 {
        bail!(
            "clone id requires an unambiguous numeric recombination id, got {} candidates for {}",
            decoded.len(),
            receptor.id
        );
    }
    let d = decoded.into_iter().next().unwrap();
    let vj_distance = u32::from(d.p_v3_len)
        + u32::from(d.n1_len)
        + u32::from(d.p_d5_len.unwrap_or(0))
        + u32::from(d.d_retained_len.unwrap_or(0))
        + u32::from(d.p_d3_len.unwrap_or(0))
        + u32::from(d.n2_len.unwrap_or(0))
        + u32::from(d.p_j5_len);
    if d.v_id > u16::MAX as usize || d.j_id > u16::MAX as usize || vj_distance > u16::MAX as u32 {
        bail!(
            "clone structural fields exceed 16-bit encoding for {}",
            receptor.id
        );
    }
    Ok(format!(
        "{:04X}{:04X}{:04X}{:04X}{:04X}",
        d.v_id, d.v_del_3, vj_distance, d.j_del_5, d.j_id
    ))
}

fn clone_id(hc: &Receptor, lc: &Receptor) -> Result<String> {
    Ok(format!(
        "CLONE:{}:{}",
        clone_chain_hex(hc)?,
        clone_chain_hex(lc)?
    ))
}

struct StageLog {
    out: BufWriter<File>,
    started: Instant,
    stage_started: Instant,
}

impl StageLog {
    fn create(path: &Path) -> Result<Self> {
        let mut out = BufWriter::new(
            File::create(path).with_context(|| format!("create {}", path.display()))?,
        );
        writeln!(out, "stage\tstage_seconds\ttotal_seconds\tdetail")?;
        Ok(Self {
            out,
            started: Instant::now(),
            stage_started: Instant::now(),
        })
    }

    fn mark(&mut self, stage: &str, detail: impl AsRef<str>) -> Result<()> {
        let now = Instant::now();
        writeln!(
            self.out,
            "{}\t{:.3}\t{:.3}\t{}",
            stage,
            now.duration_since(self.stage_started).as_secs_f64(),
            now.duration_since(self.started).as_secs_f64(),
            detail
                .as_ref()
                .replace('\t', " ")
                .replace('\n', " ")
                .replace('\r', " "),
        )?;
        self.out.flush()?;
        self.stage_started = now;
        Ok(())
    }
}

fn main() -> Result<()> {
    run(Cli::parse())
}

fn run(args: Cli) -> Result<()> {
    if !(0.0..=1.0).contains(&args.max_pearson_p) {
        bail!("--max-pearson-p must be between 0 and 1");
    }
    let status = Arc::new(RwLock::new(ClonoMapStatus::new(&args)));
    let server = spawn_status_server(
        Arc::clone(&status),
        SocketAddr::new(IpAddr::V4(Ipv4Addr::UNSPECIFIED), args.health_port),
    )?;
    let url = format!(
        "http://{}:{}",
        public_hostname(args.health_hostname.as_deref()),
        server.addr().port()
    );
    if let Ok(mut s) = status.write() {
        s.public_url = Some(url.clone());
    }
    eprintln!("[clonomap] health server: {url}");

    set_status(
        &status,
        "preparing output",
        format!("Creating {}", args.out.display()),
    );
    create_dir_all(&args.out).with_context(|| format!("create {}", args.out.display()))?;
    let mut stage_log = StageLog::create(&args.out.join("clonomap.log"))?;
    let cfg = FamilyConfig {
        max_cdr3_distance: args.max_cdr3_distance,
        min_alignment_identity: args.min_alignment_identity,
        hard_alignment_identity: args.hard_alignment_identity,
        ..FamilyConfig::default()
    };

    set_status(
        &status,
        "discovering inputs",
        "Resolving VDJ output directories",
    );
    let sources = input_sources(&args.vdj_out)?;
    set_status(
        &status,
        "loading receptors",
        format!(
            "Reading receptor calls from {} VDJ source(s)",
            sources.len()
        ),
    );
    let mut rows = Vec::new();
    for (source, dir) in &sources {
        let mut input = read_vdj_output(dir, source)?;
        println!(
            "Loaded {} receptor calls from {} ({})",
            input.len(),
            source,
            dir.display()
        );
        rows.append(&mut input);
    }
    if let Ok(mut s) = status.write() {
        s.receptor_calls = rows.len();
    }
    stage_log.mark(
        "load_receptors",
        format!("sources={} receptor_calls={}", sources.len(), rows.len()),
    )?;
    set_status(
        &status,
        "opening reference curator",
        "Loading persistent reference-candidate state",
    );
    // ClonoMap discovers reference-incompatible fragments; reference_curator
    // owns their persistent identity, resolution state and external evidence.
    let curator_store = args
        .reference_curator
        .clone()
        .unwrap_or_else(|| args.out.join("reference_curator.bin"));
    let mut novel_v = ReferenceModels::open(Some(&curator_store), &args.chemistry)
        .map_err(anyhow::Error::msg)
        .context("opening reference curator")?;

    // Learn the ordinary distance-to-supplied-germline background from the full
    // repertoire, before any family-size or CDR3 grouping can bias it. HC and LC
    // have separate backgrounds. Receptors beyond mean + 2 SD receive a
    // hypothetical V model through reference_curator, then the normal ClonoMap
    // family machinery starts from scratch with those augmented assignments.
    set_status(
        &status,
        "reference refinement",
        format!(
            "Refining {} receptor calls against supplied germline",
            rows.len()
        ),
    );
    let refinement = refine_global_reference_outliers(&mut rows, &mut novel_v);
    println!("\nGlobal reference refinement (max distance from supplied germline)");
    print_refinement_line("HC", &refinement.hc);
    print_refinement_line("LC", &refinement.lc);

    let mut novel_receptors = refinement.hc.hypothesized + refinement.lc.hypothesized;
    for row in &mut rows {
        if row.productive
            && !row.receptor.v.starts_with("RC-")
            && novel_v.resolve_receptor_for_sample(
                &mut row.receptor,
                cfg.hard_alignment_identity,
                Some(&row.source),
            )
        {
            novel_receptors += 1;
        }
    }
    stage_log.mark(
        "reference_refinement",
        format!(
            "hc_n={} hc_outliers={} hc_hypothesized={} lc_n={} lc_outliers={} lc_hypothesized={} receptors_assigned={}",
            refinement.hc.n, refinement.hc.outliers, refinement.hc.hypothesized,
            refinement.lc.n, refinement.lc.outliers, refinement.lc.hypothesized, novel_receptors
        ),
    )?;

    // Keep the normalized receptor calls for the output contract. Family inference
    // may discard ambiguous/non-functional alternatives, but cells.tsv and the
    // duplicate audit must still describe what was observed in each cell.
    let receptor_audit = rows.clone();
    let mut by_cell: BTreeMap<String, Vec<CallRow>> = BTreeMap::new();
    for row in rows {
        by_cell.entry(row.cell.clone()).or_default().push(row);
    }

    set_status(
        &status,
        "HC qualification",
        format!(
            "Selecting one valid heavy chain across {} cells",
            by_cell.len()
        ),
    );
    let examined = by_cell.len();
    let mut no_hc = 0usize;
    let mut multiple_hc = 0usize;
    let mut qualified = Vec::new();
    let mut light_only_raw: Vec<(String, Receptor)> = Vec::new();

    // LC normalization is deliberately performed while the cell is assembled.
    // Obvious non-functional alternatives must never leak into Family/plot/HMM
    // logic.  Multiple *functional* LCs are retained as candidates: resolving
    // those requires family/clone context and must not be guessed here.
    let mut lc_no_call = 0usize;
    let mut lc_one_functional = 0usize;
    let mut lc_one_nonfunctional = 0usize;
    let mut lc_multiple_one_functional = 0usize;
    let mut lc_multiple_none_functional = 0usize;
    let mut lc_multiple_functional = 0usize;
    let mut lc_nonfunctional_discarded = 0usize;

    for (cell_id, rows) in by_cell {
        let hcs: Vec<&CallRow> = rows
            .iter()
            .filter(|r| r.productive && r.receptor.chain == "IGH" && receptor_complete(&r.receptor))
            .collect();
        match hcs.as_slice() {
            [] => {
                no_hc += 1;
                for row in rows.iter().filter(|r| {
                    matches!(r.receptor.chain.as_str(), "IGK" | "IGL")
                        && r.productive
                        && receptor_complete(&r.receptor)
                }) {
                    light_only_raw.push((cell_id.clone(), row.receptor.clone()));
                }
            }
            [hc] => {
                let lc_rows: Vec<&CallRow> = rows
                    .iter()
                    .filter(|r| matches!(r.receptor.chain.as_str(), "IGK" | "IGL"))
                    .collect();
                let functional_lc: Vec<&CallRow> = lc_rows
                    .iter()
                    .copied()
                    .filter(|r| r.productive && receptor_complete(&r.receptor))
                    .collect();

                let discarded_here = lc_rows.len().saturating_sub(functional_lc.len());
                lc_nonfunctional_discarded += discarded_here;

                match (lc_rows.len(), functional_lc.len()) {
                    (0, 0) => lc_no_call += 1,
                    (1, 0) => lc_one_nonfunctional += 1,
                    (1, 1) => lc_one_functional += 1,
                    (n, 0) if n > 1 => lc_multiple_none_functional += 1,
                    (n, 1) if n > 1 => lc_multiple_one_functional += 1,
                    (_, n) if n > 1 => lc_multiple_functional += 1,
                    _ => unreachable!("LC normalization count combination"),
                }

                // Zero or one functional LC is now fully resolved.  More than
                // one remains explicit candidate evidence for later contextual
                // resolution inside the finalized HC family.
                let lc = functional_lc
                    .into_iter()
                    .map(|r| r.receptor.clone())
                    .collect();
                qualified.push(CellReceptor {
                    cell_id,
                    hc: hc.receptor.clone(),
                    lc,
                });
            }
            _ => multiple_hc += 1,
        }
    }

    // Compact LC-only observations before any family lookup. Identical observed
    // receptor sequence/state is represented once and carries all source cells.
    let mut light_only_by_sequence: BTreeMap<(String, String), LightOnlyGroup> = BTreeMap::new();
    for (cell, receptor) in light_only_raw {
        let key = (receptor.id.clone(), receptor.observed.clone());
        let observation = LightOnlyCell {
            cell,
            evidence_reads: receptor.evidence_reads,
            support_umis: receptor.support_umis,
        };
        light_only_by_sequence
            .entry(key)
            .and_modify(|g| g.cells.push(observation.clone()))
            .or_insert_with(|| LightOnlyGroup {
                receptor,
                cells: vec![observation],
            });
    }
    let light_only_groups: Vec<LightOnlyGroup> = light_only_by_sequence.into_values().collect();

    println!("HC qualification");
    println!("  cells examined:                 {examined}");
    println!("  accepted: exactly one valid HC: {}", qualified.len());
    println!("  skipped: no valid HC:           {no_hc}");
    println!("  skipped: multiple valid HCs:    {multiple_hc}");
    println!("  NOTE: apparent duplicate HCs are intentionally not collapsed; investigate later.");
    println!("LC normalization (HC-qualified cells)");
    println!("  no LC call:                            {lc_no_call}");
    println!("  one LC, functional:                    {lc_one_functional}");
    println!("  one LC, nonfunctional/incomplete:      {lc_one_nonfunctional}");
    println!("  multiple LC, exactly one functional:   {lc_multiple_one_functional}");
    println!("  multiple LC, none functional:          {lc_multiple_none_functional}");
    println!("  multiple LC, >1 functional candidates: {lc_multiple_functional}");
    println!("  nonfunctional/incomplete LC discarded: {lc_nonfunctional_discarded}");

    stage_log.mark(
        "hc_qualification",
        format!(
            "cells_examined={examined} accepted={} no_valid_hc={no_hc} multiple_valid_hc={multiple_hc}",
            qualified.len()
        ),
    )?;

    if let Ok(mut s) = status.write() {
        s.hc_cells = qualified.len();
    }
    set_status(
        &status,
        "HC family construction",
        format!("Clustering {} qualified HC cells", qualified.len()),
    );
    let mut families = build_initial_families(qualified, &cfg);
    let provisional_cells: usize = families.iter().map(Family::provisional_len).sum();
    println!(
        "\nInitial HC families: {} ({provisional_cells} cells)",
        families.len()
    );
    stage_log.mark(
        "initial_hc_families",
        format!(
            "families={} provisional_cells={provisional_cells}",
            families.len()
        ),
    )?;

    if let Ok(mut s) = status.write() {
        s.families = families.len();
    }
    set_status(
        &status,
        "HC mutation validation",
        format!("Validating {} provisional HC families", families.len()),
    );
    // Family owns only accepted members. align() returns bare CellReceptors;
    // reassignment policy and unresolved cells belong to this outer layer.
    let mut rejected = Vec::<(String, CellReceptor)>::new();
    for family in &mut families {
        let from = family.name.clone();
        rejected.extend(
            family
                .align(&cfg, &novel_v)
                .into_iter()
                .map(|cell| (from.clone(), cell)),
        );
    }
    let failed_original = rejected.len();
    let unassigned = reassign_cells(&mut families, rejected, &cfg, &novel_v);
    let reassigned = failed_original.saturating_sub(unassigned.len());
    let final_cells: usize = families.iter().map(|f| f.members.len()).sum();

    println!("\nHC post-cluster mutation validation");
    println!("  provisional cells:             {provisional_cells}");
    println!("  passed/reassigned final cells: {final_cells}");
    println!("  failed original family:        {failed_original}");
    println!("  reassigned to another family:  {reassigned}");
    println!("  remained unassigned:           {}", unassigned.len());
    stage_log.mark(
        "hc_mutation_validation",
        format!(
            "provisional_cells={provisional_cells} final_cells={final_cells} failed_original={failed_original} reassigned={reassigned} unassigned={}",
            unassigned.len()
        ),
    )?;

    if let Ok(mut s) = status.write() {
        s.hc_cells = final_cells;
    }
    set_status(
        &status,
        "LC clone assignment",
        format!(
            "Assigning light-chain clones in {} finalized HC families",
            families.len()
        ),
    );
    // LC clone ownership starts only after HC membership is final.
    for family in &mut families {
        family.split_light_chains(&cfg, &novel_v);
    }
    // Freeze mutation-agnostic clone identities only after HC and LC clustering
    // is complete. The clone id is a tandem of the deliberately lossy HC and LC
    // structural hexes, never a concatenation of the full sc-vdj IDs.
    for family in &mut families {
        for lc in &mut family.light_clones {
            lc.name = clone_id(&family.root, &lc.root)?;
        }
    }

    // Reverse lookup: a compact LC structural clone may nominate an HC family
    // only when every observed HC/LC connection agrees on that HC root.
    let mut lc_to_hc: BTreeMap<String, BTreeSet<(String, String, String)>> = BTreeMap::new();
    for family in &families {
        let hc_hex = clone_chain_hex(&family.root)?;
        for lc in &family.light_clones {
            let lc_hex = clone_chain_hex(&lc.root)?;
            lc_to_hc.entry(lc_hex).or_default().insert((
                hc_hex.clone(),
                family.name.clone(),
                lc.name.clone(),
            ));
        }
    }
    let mut light_only_assignments = Vec::new();
    let mut light_only_ambiguous = 0usize;
    let mut light_only_unmatched = 0usize;
    for group in light_only_groups {
        let lc_hex = clone_chain_hex(&group.receptor)?;
        match lc_to_hc.get(&lc_hex) {
            Some(targets) if targets.len() == 1 => {
                let (_hc_hex, family, clone_id) = targets.iter().next().unwrap().clone();
                light_only_assignments.push(LightOnlyAssignment {
                    family,
                    clone_id,
                    receptor: group.receptor,
                    cells: group.cells,
                });
            }
            Some(_) => light_only_ambiguous += group.cells.len(),
            None => light_only_unmatched += group.cells.len(),
        }
    }
    let lc_only_attached: usize = light_only_assignments.iter().map(|a| a.cells.len()).sum();

    // HC-only cells already belong to a finalized HC family. They can inherit
    // the combined HC/LC clone id only when that family has exactly one LC
    // clone. Multiple LC clones make the missing light chain genuinely
    // ambiguous; zero LC clones leave it unresolved.
    let mut hc_only_attached = 0usize;
    let mut hc_only_ambiguous = 0usize;
    let mut hc_only_unresolved = 0usize;
    for family in &families {
        let paired_cells: BTreeSet<&str> = family
            .light_clones
            .iter()
            .flat_map(|lc| lc.members.iter().map(|m| m.cell.as_str()))
            .collect();
        let hc_only = family
            .members
            .iter()
            .filter(|m| !paired_cells.contains(m.cell.cell_id.as_str()))
            .count();
        match family.light_clones.len() {
            0 => hc_only_unresolved += hc_only,
            1 => hc_only_attached += hc_only,
            _ => hc_only_ambiguous += hc_only,
        }
    }

    println!("Partial-chain clone rescue");
    println!("  HC-only cells attached:   {hc_only_attached}");
    println!("  HC-only cells ambiguous:  {hc_only_ambiguous}");
    println!("  HC-only cells unresolved: {hc_only_unresolved}");
    println!("  LC-only cells attached:   {lc_only_attached}");
    println!("  LC-only cells ambiguous:  {light_only_ambiguous}");
    println!("  LC-only cells unresolved: {light_only_unmatched}");

    let mut top_lc_rescues: Vec<&LightOnlyAssignment> = light_only_assignments.iter().collect();
    top_lc_rescues.sort_by_key(|a| std::cmp::Reverse(a.cells.len()));
    top_lc_rescues.truncate(5);
    if !top_lc_rescues.is_empty() {
        println!("  Top rescued LC-only groups:");
        for (rank, assignment) in top_lc_rescues.iter().enumerate() {
            println!(
                "    {}. {}",
                rank + 1,
                lc_rescue_summary(assignment, lc_only_attached)
            );
        }
    }
    let top_lc_rescue_summaries: Vec<String> = top_lc_rescues
        .iter()
        .map(|a| lc_rescue_summary(a, lc_only_attached))
        .collect();

    let lc_clones: usize = families
        .iter()
        .map(|family| family.light_clones.len())
        .sum();
    stage_log.mark(
        "lc_clone_assignment",
        format!(
            "families={} lc_clones={lc_clones} hc_only_attached={hc_only_attached} hc_only_ambiguous={hc_only_ambiguous} hc_only_unresolved={hc_only_unresolved} lc_only_attached={lc_only_attached} lc_only_ambiguous={light_only_ambiguous} lc_only_unresolved={light_only_unmatched}",
            families.len()
        ),
    )?;

    if let Ok(mut s) = status.write() {
        s.lc_clones = lc_clones;
        s.hc_only_attached = hc_only_attached;
        s.hc_only_ambiguous = hc_only_ambiguous;
        s.hc_only_unresolved = hc_only_unresolved;
        s.lc_only_attached = lc_only_attached;
        s.lc_only_ambiguous = light_only_ambiguous;
        s.lc_only_unresolved = light_only_unmatched;
        s.top_lc_rescues = top_lc_rescue_summaries;
    }
    set_status(
        &status,
        "writing analysis reports",
        "Writing family, cell, overlap and reference reports",
    );
    write_family_report(&args.out.join("families.tsv"), &families)?;
    let cell_id_detectors: Vec<PrimerDetector> = args
        .chemistry
        .iter()
        .copied()
        .filter_map(|chemistry| PrimerDetector::from_chemistry(chemistry).ok())
        .collect();
    let receptor_rows_tmp = args.out.join(".cells.receptor_rows.tmp.tsv");
    write_cell_report(
        &receptor_rows_tmp,
        &families,
        &light_only_assignments,
        &cell_id_detectors,
    )?;
    write_unique_cells_and_duplicates(
        &receptor_rows_tmp,
        &args.out.join("cells.tsv"),
        &args.out.join("potential_duplicates.tsv"),
        &receptor_audit,
        &cell_id_detectors,
    )?;
    std::fs::remove_file(&receptor_rows_tmp).ok();
    write_overlap_reports(&args.out, &families)?;
    write_reference_candidate_report(&args.out.join("reference_candidates.tsv"), &novel_v)?;
    novel_v
        .save(&curator_store)
        .map_err(anyhow::Error::msg)
        .with_context(|| format!("saving reference curator {}", curator_store.display()))?;
    write_unassigned(&args.out.join("hc_unassigned.tsv"), &unassigned)?;
    stage_log.mark(
        "analysis_reports",
        format!(
            "families={} final_hc_cells={} reference_models={}",
            families.len(),
            final_cells,
            novel_v.total()
        ),
    )?;

    println!(
        "\nReference-incompatible models: {} ({} receptors assigned; evidence in reference_candidates.tsv)",
        novel_v.total(),
        novel_receptors
    );

    if let Some(vault_path) = &args.vault {
        set_status(
            &status,
            "loading model vault",
            format!("Loading {}", vault_path.display()),
        );
        let vault = AaModelVault::load(vault_path)
            .with_context(|| format!("loading Ommverse model vault {}", vault_path.display()))?;
        let vault_models = vault.models.len() + vault.interpro_models.len();
        if let Ok(mut s) = status.write() {
            s.vault_models = Some(vault_models);
        }
        set_status(
            &status,
            "model-vault scoring",
            format!("Scoring family receptor states against {vault_models} models"),
        );
        write_model_space(&args.out.join("model_space"), &families, &novel_v, &vault)?;
        stage_log.mark(
            "model_vault",
            format!(
                "vault={} models={}",
                vault_path.display(),
                vault.models.len() + vault.interpro_models.len()
            ),
        )?;
    }

    if !args.no_plots {
        set_status(
            &status,
            "plotting",
            format!("Rendering family plots for {} families", families.len()),
        );
        write_family_plots(
            &args.out.join("plots"),
            &families,
            args.min_family_size,
            args.max_pearson_p,
            args.radial_layout,
        )?;
        stage_log.mark("plots", format!("families={}", families.len()))?;
    }

    let mut selected: Vec<_> = families
        .iter()
        .filter(|f| family_is_selected(f, args.min_family_size, args.max_pearson_p))
        .collect();
    selected.sort_by_key(|f| std::cmp::Reverse(f.members.len()));
    println!(
        "\nFinal HC families selected by size >= {} OR Pearson p <= {}: {}",
        args.min_family_size,
        args.max_pearson_p,
        selected.len()
    );
    println!(
        "  {:<72} {:>6}  {:>6}  {:>9}  {:>6}  {:>8}  {:>8}  {:>16}",
        "HC family", "cells", "HC max", "LC clones", "LC n", "LC mean", "LC SD", "LC top 3"
    );
    println!("  {}", "-".repeat(145));
    for f in &selected {
        let r = f.mutation_report();
        println!(
            "  {:<72} {:>6}  {:>6}  {:>9}  {:>6}  {:>8}  {:>8}  {:>16}",
            r.family,
            r.cells,
            fmt_usize_option(r.max_hc_mutations),
            r.lc_clones,
            r.lc_mutations.n,
            fmt_float_option(r.lc_mutations.mean),
            fmt_float_option(r.lc_mutations.sd),
            fmt_max3(&r.lc_mutations.max3),
        );
    }
    set_status(
        &status,
        "finalizing",
        "Preparing final family selection summary",
    );
    stage_log.mark(
        "finished",
        format!(
            "families={} selected_families={}",
            families.len(),
            selected.len()
        ),
    )?;
    if let Ok(mut s) = status.write() {
        s.finish();
    }
    println!("\nReports written to {}", args.out.display());
    println!(
        "Stage timings written to {}",
        args.out.join("clonomap.log").display()
    );
    Ok(())
}

#[derive(Default)]
struct PlotState {
    cells: BTreeSet<String>,
    isotypes: BTreeMap<String, usize>,
    lc_name: String,
    hc_depth: usize,
    lc_depth: Option<usize>,
}

#[derive(Clone, Copy, Default)]
struct PearsonStat {
    n: usize,
    r: Option<f64>,
    p: Option<f64>,
}

#[derive(Default)]
struct FamilyPlotAnalysis {
    abundance_hc: PearsonStat,
    abundance_lc: PearsonStat,
    abundance_paired: PearsonStat,
    hc_lc: PearsonStat,
    isotype_hex: String,
    light_chain_hex: String,
    hc_depth_hex: String,
    lc_depth_hex: String,
    paired_depth_hex: String,
    hclc: BTreeMap<String, HclcAnalysis>,
}

#[derive(Clone, Copy, Default)]
struct HclcAnalysis {
    abundance_hc: PearsonStat,
    abundance_lc: PearsonStat,
    abundance_paired: PearsonStat,
    hc_lc: PearsonStat,
}

fn mutation_distance(m: &MutationMeasurement) -> usize {
    m.substitutions + m.indels.iter().map(|x| x.len).sum::<usize>()
}

#[derive(Clone, Copy, Default)]
struct RefinementStats {
    n: usize,
    mean: Option<f64>,
    sd: Option<f64>,
    threshold: Option<f64>,
    outliers: usize,
    hypothesized: usize,
}

#[derive(Clone, Copy, Default)]
struct GlobalRefinement {
    hc: RefinementStats,
    lc: RefinementStats,
}

fn distance_stats(values: &[usize]) -> RefinementStats {
    if values.is_empty() {
        return RefinementStats::default();
    }
    let n = values.len();
    let mean = values.iter().map(|&x| x as f64).sum::<f64>() / n as f64;
    let variance = values
        .iter()
        .map(|&x| {
            let d = x as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n as f64;
    let sd = variance.sqrt();
    RefinementStats {
        n,
        mean: Some(mean),
        sd: Some(sd),
        threshold: Some(mean + 2.0 * sd),
        ..Default::default()
    }
}

fn supplied_germline_distance(row: &CallRow) -> Option<usize> {
    if !row.productive || !receptor_complete(&row.receptor) {
        return None;
    }
    align_fragment(&row.receptor.naive, &row.receptor.observed).map(|m| mutation_distance(&m))
}

fn refine_global_reference_outliers(
    rows: &mut [CallRow],
    models: &mut ReferenceModels,
) -> GlobalRefinement {
    // Alignment against the supplied germline dominates this stage. Compute it
    // once per receptor, in parallel, then reuse the result both for the global
    // HC/LC background and for outlier selection. Reference-model registration
    // remains serial so RC IDs and model merging stay deterministic.
    let distances: Vec<Option<usize>> = rows.par_iter().map(supplied_germline_distance).collect();

    let hc_values: Vec<usize> = rows
        .iter()
        .zip(&distances)
        .filter(|(r, _)| r.receptor.chain == "IGH")
        .filter_map(|(_, distance)| *distance)
        .collect();
    let lc_values: Vec<usize> = rows
        .iter()
        .zip(&distances)
        .filter(|(r, _)| matches!(r.receptor.chain.as_str(), "IGK" | "IGL"))
        .filter_map(|(_, distance)| *distance)
        .collect();
    let mut result = GlobalRefinement {
        hc: distance_stats(&hc_values),
        lc: distance_stats(&lc_values),
    };

    for (row, distance) in rows.iter_mut().zip(distances) {
        let stats = match row.receptor.chain.as_str() {
            "IGH" => &mut result.hc,
            "IGK" | "IGL" => &mut result.lc,
            _ => continue,
        };
        let Some(threshold) = stats.threshold else {
            continue;
        };
        let Some(distance) = distance else {
            continue;
        };
        if distance as f64 <= threshold {
            continue;
        }
        stats.outliers += 1;
        if models.hypothesize_receptor_for_sample(&mut row.receptor, Some(&row.source)) {
            stats.hypothesized += 1;
        }
    }
    result
}

fn print_refinement_line(label: &str, stats: &RefinementStats) {
    match (stats.mean, stats.sd, stats.threshold) {
        (Some(mean), Some(sd), Some(threshold)) => println!(
            "  {label}: n={} mean={mean:.2} SD={sd:.2} threshold={threshold:.2} outliers={} hypothetical models assigned={}",
            stats.n, stats.outliers, stats.hypothesized
        ),
        _ => println!("  {label}: no evaluable receptors"),
    }
}

fn pearson(x: &[f64], y: &[f64]) -> PearsonStat {
    let n = x.len().min(y.len());
    if n < 3 {
        return PearsonStat {
            n,
            ..Default::default()
        };
    }
    let (x, y) = (&x[..n], &y[..n]);
    let mx = x.iter().sum::<f64>() / n as f64;
    let my = y.iter().sum::<f64>() / n as f64;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        let dx = a - mx;
        let dy = b - my;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if sxx == 0.0 || syy == 0.0 {
        return PearsonStat {
            n,
            ..Default::default()
        };
    }
    let r = (sxy / (sxx * syy).sqrt()).clamp(-1.0, 1.0);
    let p = if r.abs() >= 1.0 {
        0.0
    } else {
        let t = r.abs() * (((n - 2) as f64) / (1.0 - r * r)).sqrt();
        let dist = StudentsT::new(0.0, 1.0, (n - 2) as f64).expect("valid Student t distribution");
        2.0 * (1.0 - dist.cdf(t))
    };
    PearsonStat {
        n,
        r: Some(r),
        p: Some(p),
    }
}

fn family_plot_analysis(family: &Family) -> FamilyPlotAnalysis {
    let hc_by_cell: HashMap<&str, (&str, usize)> = family
        .members
        .iter()
        .map(|m| {
            (
                m.cell.cell_id.as_str(),
                (m.cell.hc.c.as_str(), m.hc_mutations.mutation_events()),
            )
        })
        .collect();
    let mut paired_cells = BTreeSet::new();
    let mut states: BTreeMap<(usize, String, Option<usize>), PlotState> = BTreeMap::new();
    for lc in &family.light_clones {
        for member in &lc.members {
            let Some(&(isotype, hc_depth)) = hc_by_cell.get(member.cell.as_str()) else {
                continue;
            };
            paired_cells.insert(member.cell.clone());
            let lc_depth = member.mutations.mutation_events();
            let state = states
                .entry((hc_depth, lc.name.clone(), Some(lc_depth)))
                .or_default();
            state.cells.insert(member.cell.clone());
            *state
                .isotypes
                .entry(
                    if isotype.is_empty() {
                        "unknown"
                    } else {
                        isotype
                    }
                    .to_string(),
                )
                .or_default() += 1;
            state.lc_name = lc.name.clone();
            state.hc_depth = hc_depth;
            state.lc_depth = Some(lc_depth);
        }
    }
    for member in &family.members {
        if paired_cells.contains(&member.cell.cell_id) {
            continue;
        }
        let hc_depth = member.hc_mutations.mutation_events();
        let state = states
            .entry((hc_depth, "unpaired".to_string(), None))
            .or_default();
        state.cells.insert(member.cell.cell_id.clone());
        let iso = if member.cell.hc.c.is_empty() {
            "unknown"
        } else {
            member.cell.hc.c.as_str()
        };
        *state.isotypes.entry(iso.to_string()).or_default() += 1;
        state.lc_name = "unpaired".to_string();
        state.hc_depth = hc_depth;
    }
    let sv: Vec<PlotState> = states.into_values().collect();
    let paired: Vec<&PlotState> = sv.iter().filter(|s| s.lc_depth.is_some()).collect();
    let abundance: Vec<f64> = paired.iter().map(|s| s.cells.len() as f64).collect();
    let hc: Vec<f64> = paired.iter().map(|s| s.hc_depth as f64).collect();
    let lc: Vec<f64> = paired.iter().map(|s| s.lc_depth.unwrap() as f64).collect();
    let total: Vec<f64> = paired
        .iter()
        .map(|s| (s.hc_depth + s.lc_depth.unwrap()) as f64)
        .collect();

    let mut iso_cat = vec!["HC NAIVE".to_string()];
    let mut lc_cat = vec!["HC NAIVE".to_string()];
    let mut hd = vec![None];
    let mut ld = vec![None];
    let mut pd = vec![None];
    for s in &sv {
        iso_cat.push(
            s.isotypes
                .iter()
                .max_by(|(ka, va), (kb, vb)| va.cmp(vb).then_with(|| kb.cmp(ka)))
                .map(|(k, _)| k.clone())
                .unwrap_or_else(|| "unknown".into()),
        );
        lc_cat.push(s.lc_name.clone());
        hd.push(Some(s.hc_depth as f32));
        ld.push(s.lc_depth.map(|x| x as f32));
        pd.push(Some((s.hc_depth + s.lc_depth.unwrap_or(0)) as f32));
    }
    let mut lc_order: Vec<String> = family.light_clones.iter().map(|x| x.name.clone()).collect();
    lc_order.sort_by_key(|name| {
        std::cmp::Reverse(
            family
                .light_clones
                .iter()
                .find(|x| &x.name == name)
                .map(|x| x.members.len())
                .unwrap_or(0),
        )
    });
    lc_order.dedup();
    if lc_cat.iter().any(|x| x == "unpaired") {
        lc_order.push("unpaired".into());
    }
    let mut hclc = BTreeMap::new();
    for lc_name in family.light_clones.iter().map(|x| x.name.clone()) {
        let subset: Vec<&PlotState> = paired
            .iter()
            .copied()
            .filter(|s| s.lc_name == lc_name)
            .collect();
        let abundance: Vec<f64> = subset.iter().map(|s| s.cells.len() as f64).collect();
        let hc: Vec<f64> = subset.iter().map(|s| s.hc_depth as f64).collect();
        let lc: Vec<f64> = subset.iter().map(|s| s.lc_depth.unwrap() as f64).collect();
        let total: Vec<f64> = subset
            .iter()
            .map(|s| (s.hc_depth + s.lc_depth.unwrap()) as f64)
            .collect();
        hclc.insert(
            lc_name,
            HclcAnalysis {
                abundance_hc: pearson(&abundance, &hc),
                abundance_lc: pearson(&abundance, &lc),
                abundance_paired: pearson(&abundance, &total),
                hc_lc: pearson(&hc, &lc),
            },
        );
    }
    let iso_colors = rooted_categorical_hex(&iso_cat, None, 0);
    let lc_colors = rooted_categorical_hex(&lc_cat, Some(&lc_order), 0);
    let hd_colors = rooted_continuous_hex(&hd, 0);
    let ld_colors = rooted_continuous_hex(&ld, 0);
    let pd_colors = rooted_continuous_hex(&pd, 0);
    let encode = |colors: &[String]| -> String {
        sv.iter()
            .enumerate()
            .flat_map(|(i, s)| {
                s.cells
                    .iter()
                    .map(move |cell| format!("{}={}", cell, colors[i + 1]))
            })
            .collect::<Vec<_>>()
            .join(";")
    };
    FamilyPlotAnalysis {
        abundance_hc: pearson(&abundance, &hc),
        abundance_lc: pearson(&abundance, &lc),
        abundance_paired: pearson(&abundance, &total),
        hc_lc: pearson(&hc, &lc),
        isotype_hex: encode(&iso_colors),
        light_chain_hex: encode(&lc_colors),
        hc_depth_hex: encode(&hd_colors),
        lc_depth_hex: encode(&ld_colors),
        paired_depth_hex: encode(&pd_colors),
        hclc,
    }
}

fn stats_significant(stats: &[PearsonStat], max_p: f64) -> bool {
    stats.iter().any(|s| s.p.is_some_and(|p| p <= max_p))
}

fn family_is_significant(a: &FamilyPlotAnalysis, max_p: f64) -> bool {
    stats_significant(
        &[a.abundance_hc, a.abundance_lc, a.abundance_paired, a.hc_lc],
        max_p,
    ) || a.hclc.values().any(|x| {
        stats_significant(
            &[x.abundance_hc, x.abundance_lc, x.abundance_paired, x.hc_lc],
            max_p,
        )
    })
}

fn family_is_selected(family: &Family, min_size: usize, max_p: f64) -> bool {
    family.members.len() >= min_size || family_is_significant(&family_plot_analysis(family), max_p)
}

fn write_family_plots(
    out: &Path,
    families: &[Family],
    min_size: usize,
    max_p: f64,
    radial_layout: bool,
) -> Result<()> {
    create_dir_all(out).with_context(|| format!("create {}", out.display()))?;
    let mut hc_plot_name_counts: HashMap<String, usize> = HashMap::new();
    for family in families
        .iter()
        .filter(|f| family_is_selected(f, min_size, max_p))
    {
        // The HC recombination ID is the stable biological identity of the plotted
        // family. Only add a suffix if more than one selected family has the same
        // HC root ID; the first keeps the clean unsuffixed name.
        let count = hc_plot_name_counts
            .entry(family.root.id.clone())
            .or_default();
        *count += 1;
        let plot_dir_name = if *count == 1 {
            family.root.id.clone()
        } else {
            format!("{}_{}", family.root.id, count)
        };
        write_one_family_plot(out, family, &plot_dir_name, min_size, max_p, radial_layout)?;
    }
    Ok(())
}

#[derive(Debug)]
struct ChainStateValidation {
    status: &'static str,
    receptor_ids: Vec<String>,
    representative_nt: Option<String>,
    detail: String,
}

fn indel_signature(m: &MutationMeasurement) -> String {
    m.indels
        .iter()
        .map(|x| {
            format!(
                "{}:{}:{}",
                x.naive_pos,
                if x.inserted { "ins" } else { "del" },
                x.len
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn validate_chain_state<'a, I>(reference: &str, receptors: I) -> ChainStateValidation
where
    I: IntoIterator<Item = (&'a Receptor, &'a MutationMeasurement)>,
{
    let rows: Vec<_> = receptors.into_iter().collect();
    if rows.is_empty() {
        return ChainStateValidation {
            status: "absent",
            receptor_ids: Vec::new(),
            representative_nt: None,
            detail: "no receptor assigned to graph state".into(),
        };
    }
    let receptor_ids = rows.iter().map(|(r, _)| r.id.clone()).collect::<Vec<_>>();
    let mut merged = reference
        .as_bytes()
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect::<Vec<_>>();
    let mut observed_at: BTreeMap<usize, u8> = BTreeMap::new();
    let mut conflicts = Vec::new();
    let mut projection_failures = Vec::new();
    let mut indels = BTreeSet::new();
    for (receptor, mutations) in &rows {
        indels.insert(indel_signature(mutations));
        let Some(projected) = project_observed_bases(reference, &receptor.observed) else {
            projection_failures.push(receptor.id.clone());
            continue;
        };
        for (pos, base) in projected {
            if let Some(previous) = observed_at.get(&pos) {
                if *previous != base {
                    conflicts.push(format!(
                        "{}:{}>{}/{}",
                        pos,
                        reference.as_bytes().get(pos).copied().unwrap_or(b'?') as char,
                        *previous as char,
                        base as char
                    ));
                    continue;
                }
            } else {
                observed_at.insert(pos, base);
                if pos < merged.len() {
                    merged[pos] = base;
                }
            }
        }
    }
    if !conflicts.is_empty() {
        return ChainStateValidation {
            status: "CONFLICT",
            receptor_ids,
            representative_nt: None,
            detail: format!(
                "overlapping receptor observations disagree at {}",
                conflicts.join(",")
            ),
        };
    }
    if !projection_failures.is_empty() {
        return ChainStateValidation {
            status: "UNRESOLVED",
            receptor_ids,
            representative_nt: None,
            detail: format!(
                "could not place receptor(s) on clone reference: {}",
                projection_failures.join(",")
            ),
        };
    }
    // Different non-empty indel signatures cannot yet be merged into one exact
    // nucleotide representative.  Report that explicitly instead of hiding it
    // behind the substitution-only projection used by the AA/HMM prototype.
    let nonempty_indels = indels.iter().filter(|x| !x.is_empty()).count();
    if nonempty_indels > 1 || (nonempty_indels == 1 && indels.contains("")) {
        return ChainStateValidation {
            status: "INDEL_AMBIGUOUS",
            receptor_ids,
            representative_nt: None,
            detail: format!(
                "clone-state members have different observed indel evidence: {}",
                indels.into_iter().collect::<Vec<_>>().join(" | ")
            ),
        };
    }
    ChainStateValidation {
        status: if rows.len() == 1 {
            "singleton"
        } else {
            "compatible"
        },
        receptor_ids,
        representative_nt: String::from_utf8(merged).ok(),
        detail: format!(
            "{} receptor observation(s) merge without contradictory overlapping bases",
            rows.len()
        ),
    }
}

fn write_graph_state_validation(dir: &Path, family: &Family, states: &[PlotState]) -> Result<()> {
    let hc_by_cell: HashMap<&str, &AlignedCell> = family
        .members
        .iter()
        .map(|m| (m.cell.cell_id.as_str(), m))
        .collect();
    let lc_by_name: HashMap<&str, &LightClone> = family
        .light_clones
        .iter()
        .map(|lc| (lc.name.as_str(), lc))
        .collect();
    let mut out = BufWriter::new(File::create(dir.join("graph_state_validation.tsv"))?);
    writeln!(
        out,
        "family\tstate\tlc_clone\thc_depth\tlc_depth\tn_cells\tcells\thc_status\thc_recombination_ids\thc_representative_nt\thc_detail\tlc_status\tlc_recombination_ids\tlc_representative_nt\tlc_detail"
    )?;
    for (i, state) in states.iter().enumerate() {
        let hc_rows = state
            .cells
            .iter()
            .filter_map(|cell| {
                hc_by_cell
                    .get(cell.as_str())
                    .map(|m| (&m.cell.hc, &m.hc_mutations))
            })
            .collect::<Vec<_>>();
        let hc = validate_chain_state(&family.root.naive, hc_rows);
        let lc = if state.lc_name == "unpaired" {
            ChainStateValidation {
                status: "absent",
                receptor_ids: Vec::new(),
                representative_nt: None,
                detail: "HC-only graph state".into(),
            }
        } else if let Some(clone) = lc_by_name.get(state.lc_name.as_str()) {
            let by_cell: HashMap<&str, &LightMember> =
                clone.members.iter().map(|m| (m.cell.as_str(), m)).collect();
            let rows = state
                .cells
                .iter()
                .filter_map(|cell| {
                    by_cell
                        .get(cell.as_str())
                        .map(|m| (&m.receptor, &m.mutations))
                })
                .collect::<Vec<_>>();
            validate_chain_state(&clone.root.naive, rows)
        } else {
            ChainStateValidation {
                status: "UNRESOLVED",
                receptor_ids: Vec::new(),
                representative_nt: None,
                detail: "LC clone named by plot state is absent from Family".into(),
            }
        };
        let clean = |x: &str| x.replace('\t', " ").replace('\n', " ");
        writeln!(
            out,
            "{}\tstate{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            family.name,
            i + 1,
            state.lc_name,
            state.hc_depth,
            state
                .lc_depth
                .map(|x| x.to_string())
                .unwrap_or_else(|| "NA".into()),
            state.cells.len(),
            state.cells.iter().cloned().collect::<Vec<_>>().join(","),
            hc.status,
            hc.receptor_ids.join(","),
            hc.representative_nt.as_deref().unwrap_or(""),
            clean(&hc.detail),
            lc.status,
            lc.receptor_ids.join(","),
            lc.representative_nt.as_deref().unwrap_or(""),
            clean(&lc.detail)
        )?;
    }
    Ok(())
}

fn write_one_family_plot(
    out: &Path,
    family: &Family,
    plot_dir_name: &str,
    _min_size: usize,
    _max_p: f64,
    radial_layout: bool,
) -> Result<()> {
    // Plotting consumes the frozen Family.  It deliberately uses the mutation
    // measurements already accepted by HC/LC family logic; it never realigns a
    // receptor or changes membership.  Until span-aware positional mutation
    // coordinates are available, topology uses the honest comparable quantities
    // we have for every partial receptor: HC and LC mutation-event depth.
    let hc_by_cell: HashMap<&str, (&str, usize)> = family
        .members
        .iter()
        .map(|m| {
            (
                m.cell.cell_id.as_str(),
                (m.cell.hc.c.as_str(), m.hc_mutations.mutation_events()),
            )
        })
        .collect();
    let mut paired_cells = BTreeSet::new();
    let mut states: BTreeMap<(usize, String, Option<usize>), PlotState> = BTreeMap::new();

    for lc in &family.light_clones {
        for member in &lc.members {
            let Some(&(isotype, hc_depth)) = hc_by_cell.get(member.cell.as_str()) else {
                continue;
            };
            paired_cells.insert(member.cell.clone());
            let lc_depth = member.mutations.mutation_events();
            let key = (hc_depth, lc.name.clone(), Some(lc_depth));
            let state = states.entry(key).or_default();
            state.cells.insert(member.cell.clone());
            *state
                .isotypes
                .entry(
                    if isotype.is_empty() {
                        "unknown"
                    } else {
                        isotype
                    }
                    .to_string(),
                )
                .or_default() += 1;
            state.lc_name = lc.name.clone();
            state.hc_depth = hc_depth;
            state.lc_depth = Some(lc_depth);
        }
    }
    for member in &family.members {
        if paired_cells.contains(&member.cell.cell_id) {
            continue;
        }
        let hc_depth = member.hc_mutations.mutation_events();
        let key = (hc_depth, "unpaired".to_string(), None);
        let state = states.entry(key).or_default();
        state.cells.insert(member.cell.cell_id.clone());
        let iso = if member.cell.hc.c.is_empty() {
            "unknown"
        } else {
            member.cell.hc.c.as_str()
        };
        *state.isotypes.entry(iso.to_string()).or_default() += 1;
        state.lc_name = "unpaired".to_string();
        state.hc_depth = hc_depth;
    }
    if states.len() < 2 {
        return Ok(());
    }

    const ROOT: &str = "HC NAIVE";
    let state_values: Vec<PlotState> = states.into_values().collect();
    let n = state_values.len() + 1;
    let mut labels = Vec::with_capacity(n);
    labels.push(ROOT.to_string());
    labels.extend(state_values.iter().enumerate().map(|(i, s)| {
        format!(
            "state{}|HC{}|{}|LC{}",
            i + 1,
            s.hc_depth,
            s.lc_name,
            s.lc_depth
                .map(|x| x.to_string())
                .unwrap_or_else(|| "NA".into())
        )
    }));
    let mut features = Array2::<f32>::zeros((n, 2));
    let mut groups = Vec::with_capacity(n);
    groups.push(ROOT.to_string());
    for (i, s) in state_values.iter().enumerate() {
        features[[i + 1, 0]] = s.hc_depth as f32;
        features[[i + 1, 1]] = s.lc_depth.unwrap_or(0) as f32;
        groups.push(s.lc_name.clone());
    }
    let model = ClonoMap::from_grouped_feature_matrix(labels, features, groups, 0, 1, 2)
        .map_err(|e| anyhow::anyhow!("ClonoMap plotting model failed for {}: {e}", family.name))?;

    let mut abundance = vec![1usize];
    let mut iso_cat = vec![ROOT.to_string()];
    let mut iso_mixed = vec![false];
    let mut lc_cat = vec![ROOT.to_string()];
    let mut lc_mixed = vec![false];
    let mut hc_depth = vec![None];
    let mut lc_depth = vec![None];
    let mut paired_depth = vec![None];
    for s in &state_values {
        abundance.push(s.cells.len());
        let dominant_iso = s
            .isotypes
            .iter()
            .max_by(|(ka, va), (kb, vb)| va.cmp(vb).then_with(|| kb.cmp(ka)))
            .map(|(k, _)| k.clone())
            .unwrap_or_else(|| "unknown".to_string());
        iso_cat.push(dominant_iso);
        iso_mixed.push(s.isotypes.len() > 1);
        lc_cat.push(s.lc_name.clone());
        lc_mixed.push(false);
        hc_depth.push(Some(s.hc_depth as f32));
        lc_depth.push(s.lc_depth.map(|x| x as f32));
        paired_depth.push(Some((s.hc_depth + s.lc_depth.unwrap_or(0)) as f32));
    }
    let mut lc_order: Vec<String> = family.light_clones.iter().map(|x| x.name.clone()).collect();
    lc_order.sort_by_key(|name| {
        std::cmp::Reverse(
            family
                .light_clones
                .iter()
                .find(|x| &x.name == name)
                .map(|x| x.members.len())
                .unwrap_or(0),
        )
    });
    lc_order.dedup();
    if lc_cat.iter().any(|x| x == "unpaired") {
        lc_order.push("unpaired".to_string());
    }

    // Use the HC recombination ID as the filesystem identity as well as the
    // biological identity. Unlike family.name (V/J/CDR3), the compact HC:<HEX>
    // ID remains comfortably below filesystem component limits.
    let dir = out.join(plot_dir_name);
    create_dir_all(&dir).with_context(|| format!("create {}", dir.display()))?;
    write_graph_state_validation(&dir, family, &state_values)?;
    let title = format!(
        "{} | {} cells | {} LC clones",
        family.name,
        family.members.len(),
        family.light_clones.len()
    );
    let path = |name: &str| dir.join(name).to_string_lossy().into_owned();
    model
        .tree
        .plot_rooted_annotated_cached(
            n,
            0,
            &iso_cat,
            None,
            &iso_mixed,
            &abundance,
            "IGH constant class",
            &title,
            &path("mst_rooted_isotype.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("isotype SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_annotated_cached(
            n,
            0,
            &lc_cat,
            Some(&lc_order),
            &lc_mixed,
            &abundance,
            "Light-chain clone",
            &title,
            &path("mst_rooted_light_chain.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("LC SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_continuous_cached(
            n,
            0,
            &hc_depth,
            &abundance,
            "HC mutational depth (events)",
            &title,
            &path("mst_rooted_hc_depth.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("HC-depth SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_continuous_cached(
            n,
            0,
            &lc_depth,
            &abundance,
            "LC mutational depth (events)",
            &title,
            &path("mst_rooted_lc_depth.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("LC-depth SVG failed for {}: {e}", family.name))?;
    model
        .tree
        .plot_rooted_continuous_cached(
            n,
            0,
            &paired_depth,
            &abundance,
            "Paired HC+LC mutational depth (events)",
            &title,
            &path("mst_rooted_paired_depth.svg"),
            radial_layout,
        )
        .map_err(|e| anyhow::anyhow!("paired-depth SVG failed for {}: {e}", family.name))?;
    Ok(())
}

fn receptor_complete(r: &Receptor) -> bool {
    !r.v.is_empty()
        && !r.j.is_empty()
        && !r.cdr3_nt.is_empty()
        && !r.naive.is_empty()
        && !r.observed.is_empty()
}

fn build_initial_families(mut cells: Vec<CellReceptor>, cfg: &FamilyConfig) -> Vec<Family> {
    cells.sort_by(|a, b| {
        (&a.hc.v, &a.hc.j, &a.hc.cdr3_nt, &a.cell_id).cmp(&(
            &b.hc.v,
            &b.hc.j,
            &b.hc.cdr3_nt,
            &b.cell_id,
        ))
    });
    let mut families: Vec<Family> = Vec::new();
    for cell in cells {
        let mut pending = Some(cell);
        for family in &mut families {
            let c = pending.take().unwrap();
            match family.add_candidate(c, cfg) {
                Ok(()) => {
                    pending = None;
                    break;
                }
                Err(c) => pending = Some(c),
            }
        }
        if let Some(seed) = pending {
            let name = format!("HC:{}:{}:CDR3:{}", seed.hc.v, seed.hc.j, seed.hc.cdr3_nt);
            families.push(Family::new(name, seed));
        }
    }
    families
}

fn reassign_cells(
    families: &mut [Family],
    rejected: Vec<(String, CellReceptor)>,
    cfg: &FamilyConfig,
    novel_v: &ReferenceModels,
) -> Vec<UnassignedCell> {
    let mut unassigned = Vec::new();
    for (from_family, rejected_cell) in rejected {
        let mut cell = Some(rejected_cell);
        for family in families.iter_mut() {
            if family.name == from_family {
                continue;
            }
            let candidate = cell.take().unwrap();
            match family.try_integrate(candidate, cfg, novel_v) {
                Ok(()) => {
                    cell = None;
                    break;
                }
                Err(candidate) => cell = Some(candidate),
            }
        }
        if let Some(cell) = cell {
            unassigned.push(UnassignedCell { cell, from_family });
        }
    }
    unassigned
}

const CELL_SCOPE_SEPARATOR: char = '\u{1f}';

fn scoped_cell_id(source: &str, cell: &str) -> String {
    format!("{source}{CELL_SCOPE_SEPARATOR}{cell}")
}

fn split_scoped_cell_id(cell: &str) -> (&str, &str) {
    cell.split_once(CELL_SCOPE_SEPARATOR)
        .unwrap_or(("input", cell))
}

fn source_label(dir: &Path) -> String {
    let name = dir.file_name().and_then(|x| x.to_str()).unwrap_or("input");
    if name == "vdj_out" {
        dir.parent()
            .and_then(Path::parent)
            .and_then(Path::file_name)
            .and_then(|x| x.to_str())
            .or_else(|| {
                dir.parent()
                    .and_then(Path::file_name)
                    .and_then(|x| x.to_str())
            })
            .unwrap_or("input")
            .to_string()
    } else {
        name.to_string()
    }
}

fn input_sources(paths: &[PathBuf]) -> Result<Vec<(String, PathBuf)>> {
    let mut seen = BTreeSet::new();
    let mut out = Vec::with_capacity(paths.len());
    for path in paths {
        let source = source_label(path);
        if !seen.insert(source.clone()) {
            bail!(
                "multiple --vdj-out inputs resolve to source label {source:?}; use distinct sample directories"
            );
        }
        out.push((source, path.clone()));
    }
    Ok(out)
}

fn read_vdj_output(dir: &Path, source: &str) -> Result<Vec<CallRow>> {
    let airr_path = dir.join("airr_rearrangements.tsv");
    let calls_path = dir.join("vdj_calls.tsv");
    if !airr_path.is_file() {
        bail!("missing {}", airr_path.display());
    }
    if !calls_path.is_file() {
        bail!("missing {}", calls_path.display());
    }

    // naive_recombination is a per-cell reconstruction span. The compact
    // recombination ID may legitimately be shared by hundreds of cells, so the
    // join key is (cell, recombination_id), never recombination_id alone.
    let reconstruction_by_cell_and_id = read_reconstructions(&calls_path)?;

    let f = File::open(&airr_path).with_context(|| format!("open {}", airr_path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().context("empty airr_rearrangements.tsv")??;
    let names: Vec<&str> = header.split('\t').collect();
    let index: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let req = |name: &str| -> Result<usize> {
        index.get(name).copied().with_context(|| {
            format!(
                "{}: missing required AIRR column {name}",
                airr_path.display()
            )
        })
    };
    let cell_i = req("cell_id")?;
    let id_i = req("lumrik_recombination_id")?;
    let chain_i = req("locus")?;
    let v_i = req("v_call")?;
    let d_i = req("d_call")?;
    let j_i = req("j_call")?;
    let c_i = req("c_call")?;
    let cdr3_i = index
        .get("cdr3")
        .copied()
        .or_else(|| index.get("junction").copied())
        .with_context(|| format!("{}: missing cdr3/junction", airr_path.display()))?;
    let observed_i = req("sequence")?;
    let productive_i = req("productive")?;
    let max_i = [
        cell_i,
        id_i,
        chain_i,
        v_i,
        d_i,
        j_i,
        c_i,
        cdr3_i,
        observed_i,
        productive_i,
    ]
    .into_iter()
    .max()
    .unwrap();

    let mut out = Vec::new();
    let mut missing_naive = 0usize;
    for (line_no, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= max_i {
            bail!(
                "{}:{} has {} fields, need at least {}",
                airr_path.display(),
                line_no + 2,
                fields.len(),
                max_i + 1
            );
        }
        let cell = fields[cell_i].trim();
        let id = fields[id_i].trim();
        let key = (cell.to_string(), id.to_string());
        let Some(reconstruction) = reconstruction_by_cell_and_id.get(&key) else {
            missing_naive += 1;
            continue;
        };
        out.push(CallRow {
            cell: scoped_cell_id(source, cell),
            source: source.to_string(),
            productive: productive_value(fields[productive_i]),
            receptor: Receptor {
                id: id.to_string(),
                chain: fields[chain_i].to_string(),
                v: fields[v_i].to_string(),
                d: fields[d_i].to_string(),
                j: fields[j_i].to_string(),
                c: fields[c_i].to_string(),
                cdr3_nt: fields[cdr3_i].to_string(),
                naive: reconstruction.naive.clone(),
                observed: reconstruction.observed_receptor.clone(),
                alternative_reconstructions: reconstruction.alternatives,
                evidence_reads: reconstruction.evidence_reads,
                support_umis: reconstruction.support_umis,
            },
        });
    }
    if missing_naive > 0 {
        eprintln!(
            "ClonoMap: skipped {missing_naive} AIRR rows without matching (cell, recombination_id) naive_recombination in vdj_calls.tsv"
        );
    }
    Ok(out)
}

#[derive(Debug, Clone)]
struct SelectedReconstruction {
    naive: String,
    observed_receptor: String,
    alternatives: usize,
    evidence_reads: usize,
    support_umis: usize,
}

#[derive(Debug, Clone)]
struct ReconstructionCandidate {
    naive: String,
    observed_rearrangement: String,
    observed_receptor: String,
    productive: bool,
    mutation_events: usize,
    informative_pairs: usize,
    support_features: usize,
    support_umis: usize,
    junction_support_reads: usize,
    junction_spanning_reads: usize,
    receptor_rediscovery_reads: usize,
}

/// Collapse raw sc-vdj rows to one receptor reconstruction per
/// (cell, recombination_id) before CellReceptor construction.
///
/// A compact recombination ID describes the recombination structure, but the
/// same cell can carry competing nucleotide reconstructions for that structure.
/// We deliberately resolve that ambiguity here, outside the Family model:
///   1. productive beats non-productive;
///   2. lower observed-vs-naive mutation burden wins;
///   3. stronger junction/feature/rediscovery evidence breaks ties;
///   4. longer informative alignment, then lexical sequence order make the
///      result deterministic.
fn read_reconstructions(path: &Path) -> Result<HashMap<(String, String), SelectedReconstruction>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().context("empty vdj_calls.tsv")??;
    let names: Vec<&str> = header.split('\t').collect();
    let index: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let req = |name: &str| -> Result<usize> {
        index
            .get(name)
            .copied()
            .with_context(|| format!("{}: missing required column {name}", path.display()))
    };
    let cell_i = req("cell")?;
    let id_i = req("recombination_id")?;
    let naive_i = req("naive_recombination")?;
    let rearr_i = req("observed_rearrangement")?;
    let receptor_i = req("observed_receptor_sequence")?;
    let productivity_i = req("productivity_status")?;
    let support_i = req("support_features")?;
    let support_umis_i = req("support_umis")?;
    let junction_support_i = req("junction_support_reads")?;
    let junction_spanning_i = req("junction_spanning_reads")?;
    let rediscovery_i = req("receptor_rediscovery_reads")?;
    let max_i = [
        cell_i,
        id_i,
        naive_i,
        rearr_i,
        receptor_i,
        productivity_i,
        support_i,
        support_umis_i,
        junction_support_i,
        junction_spanning_i,
        rediscovery_i,
    ]
    .into_iter()
    .max()
    .unwrap();

    // Parsing stays serial so malformed input retains precise line-number errors.  The
    // expensive naive-vs-observed alignment is deliberately deferred until after
    // parsing: every candidate is independent and can be evaluated by Rayon.
    let mut parsed = Vec::new();
    for (line_no, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= max_i {
            bail!(
                "{}:{} has {} fields, need at least {}",
                path.display(),
                line_no + 2,
                fields.len(),
                max_i + 1
            );
        }
        let cell = fields[cell_i].trim();
        let id = fields[id_i].trim();
        let naive = fields[naive_i].trim();
        let observed_rearrangement = fields[rearr_i].trim();
        let observed_receptor = fields[receptor_i].trim();
        if cell.is_empty()
            || id.is_empty()
            || naive.is_empty()
            || observed_rearrangement.is_empty()
            || observed_receptor.is_empty()
        {
            continue;
        }
        parsed.push((
            (cell.to_string(), id.to_string()),
            ReconstructionCandidate {
                naive: naive.to_string(),
                observed_rearrangement: observed_rearrangement.to_string(),
                observed_receptor: observed_receptor.to_string(),
                productive: productive_value(fields[productivity_i]),
                mutation_events: usize::MAX,
                informative_pairs: 0,
                support_features: parse_count(fields[support_i]),
                support_umis: parse_count(fields[support_umis_i]),
                junction_support_reads: parse_count(fields[junction_support_i]),
                junction_spanning_reads: parse_count(fields[junction_spanning_i]),
                receptor_rediscovery_reads: parse_count(fields[rediscovery_i]),
            },
        ));
    }

    parsed.par_iter_mut().for_each(|(_, candidate)| {
        if let Some(measurement) =
            align_fragment(&candidate.naive, &candidate.observed_rearrangement)
        {
            candidate.mutation_events = measurement.mutation_events();
            candidate.informative_pairs = measurement.informative_pairs;
        }
    });

    let mut candidates: HashMap<(String, String), Vec<ReconstructionCandidate>> = HashMap::new();
    for (key, candidate) in parsed {
        candidates.entry(key).or_default().push(candidate);
    }

    let mut out = HashMap::new();
    let mut ambiguous_keys = 0usize;
    let mut discarded_rows = 0usize;
    for (key, mut rows) in candidates {
        let alternatives = rows.len();
        if alternatives > 1 {
            ambiguous_keys += 1;
            discarded_rows += alternatives - 1;
        }
        rows.sort_by(|a, b| reconstruction_rank(a).cmp(&reconstruction_rank(b)));
        let best = rows.remove(0);
        let evidence_reads = best
            .receptor_rediscovery_reads
            .max(best.junction_spanning_reads)
            .max(best.junction_support_reads);
        out.insert(
            key,
            SelectedReconstruction {
                naive: best.naive,
                observed_receptor: best.observed_receptor,
                alternatives,
                evidence_reads,
                support_umis: best.support_umis,
            },
        );
    }
    eprintln!(
        "ClonoMap receptor normalization: {ambiguous_keys} cell/recombination calls had alternative reconstructions; selected one and discarded {discarded_rows} alternatives"
    );
    Ok(out)
}

fn reconstruction_rank(
    c: &ReconstructionCandidate,
) -> (
    bool,
    usize,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    String,
    String,
) {
    (
        !c.productive,
        c.mutation_events,
        std::cmp::Reverse(c.junction_support_reads),
        std::cmp::Reverse(c.junction_spanning_reads),
        std::cmp::Reverse(c.support_features),
        std::cmp::Reverse(c.receptor_rediscovery_reads),
        std::cmp::Reverse(c.informative_pairs),
        c.naive.clone(),
        c.observed_receptor.clone(),
    )
}

fn parse_count(s: &str) -> usize {
    s.trim().parse().unwrap_or(0)
}

fn productive_value(s: &str) -> bool {
    matches!(
        s.trim().to_ascii_lowercase().as_str(),
        "productive" | "true" | "t" | "yes" | "1"
    )
}

fn fmt_usize_option(value: Option<usize>) -> String {
    value.map_or_else(|| "-".to_string(), |v| v.to_string())
}

fn fmt_float_option(value: Option<f64>) -> String {
    value.map_or_else(|| "-".into(), |x| format!("{x:.2}"))
}

fn fmt_max3(values: &[usize]) -> String {
    if values.is_empty() {
        return "[]".into();
    }
    format!(
        "[{}]",
        values
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn write_reference_candidate_report(path: &Path, models: &ReferenceModels) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "candidate_id\tchain\toriginal_v\tobservations_this_run\tevidence_reads_this_run\tresolved\tsequence"
    )?;
    for entry in models.entries() {
        let candidate = models
            .curator()
            .candidate(&entry.id)
            .expect("session candidate must exist in curator");
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            entry.id,
            entry.chain,
            entry.original_v,
            entry.observations,
            entry.evidence_reads,
            candidate.resolved.is_some(),
            String::from_utf8_lossy(&candidate.sequence)
        )?;
    }
    Ok(())
}

fn hc_c_counts(family: &Family) -> String {
    let mut counts = BTreeMap::<String, usize>::new();
    for member in &family.members {
        let c = if member.cell.hc.c.is_empty() {
            "unknown"
        } else {
            member.cell.hc.c.as_str()
        };
        *counts.entry(c.to_string()).or_default() += 1;
    }
    counts
        .into_iter()
        .map(|(c, n)| format!("{c}={n}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn write_family_report(path: &Path, families: &[Family]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "family\tcells\thc_c_counts\tmax_hc_mutations\tlc_clones\tlc_n\tmean_lc_mutations\tsd_lc_mutations\tlc_max_1\tlc_max_2\tlc_max_3\tpearson_abundance_hc_n\tpearson_abundance_hc_r\tpearson_abundance_hc_p\tpearson_abundance_lc_n\tpearson_abundance_lc_r\tpearson_abundance_lc_p\tpearson_abundance_paired_n\tpearson_abundance_paired_r\tpearson_abundance_paired_p\tpearson_hc_lc_n\tpearson_hc_lc_r\tpearson_hc_lc_p\tmst_rooted_isotype_svg_hex\tmst_rooted_light_chain_svg_hex\tmst_rooted_hc_depth_svg_hex\tmst_rooted_lc_depth_svg_hex\tmst_rooted_paired_depth_svg_hex"
    )?;
    let mut rows: Vec<_> = families.iter().collect();
    rows.sort_by_key(|f| std::cmp::Reverse(f.members.len()));
    for f in rows {
        let r = f.mutation_report();
        let a = family_plot_analysis(f);
        let mut max3 = r.lc_mutations.max3.iter().rev();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.family,
            r.cells,
            hc_c_counts(f),
            r.max_hc_mutations.map_or("NA".into(), |x| x.to_string()),
            r.lc_clones,
            r.lc_mutations.n,
            r.lc_mutations
                .mean
                .map_or("NA".into(), |x| format!("{x:.4}")),
            r.lc_mutations.sd.map_or("NA".into(), |x| format!("{x:.4}")),
            max3.next().map_or("NA".into(), |x| x.to_string()),
            max3.next().map_or("NA".into(), |x| x.to_string()),
            max3.next().map_or("NA".into(), |x| x.to_string()),
            a.abundance_hc.n,
            fmt_stat(a.abundance_hc.r),
            fmt_stat(a.abundance_hc.p),
            a.abundance_lc.n,
            fmt_stat(a.abundance_lc.r),
            fmt_stat(a.abundance_lc.p),
            a.abundance_paired.n,
            fmt_stat(a.abundance_paired.r),
            fmt_stat(a.abundance_paired.p),
            a.hc_lc.n,
            fmt_stat(a.hc_lc.r),
            fmt_stat(a.hc_lc.p),
            a.isotype_hex,
            a.light_chain_hex,
            a.hc_depth_hex,
            a.lc_depth_hex,
            a.paired_depth_hex,
        )?;
    }
    Ok(())
}

fn color_for(encoded: &str, cell: &str) -> String {
    encoded
        .split(';')
        .find_map(|x| {
            let (id, color) = x.split_once('=')?;
            (id == cell).then(|| color.to_string())
        })
        .unwrap_or_else(|| "NA".into())
}

fn write_cell_report(
    path: &Path,
    families: &[Family],
    light_only: &[LightOnlyAssignment],
    cell_id_detectors: &[PrimerDetector],
) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "source\tcell\tbd_cell_id\tfamily\tlc_clone\tclone_membership\thc_recombination_id\tlc_recombination_id\thc_c_call\thc_mutation_count\tlc_mutation_count\ttotal_mutation_count\thc_distance\tlc_distance\ttotal_distance\tfamily_pearson_abundance_hc_n\tfamily_pearson_abundance_hc_r\tfamily_pearson_abundance_hc_p\tfamily_pearson_abundance_lc_n\tfamily_pearson_abundance_lc_r\tfamily_pearson_abundance_lc_p\tfamily_pearson_abundance_paired_n\tfamily_pearson_abundance_paired_r\tfamily_pearson_abundance_paired_p\tfamily_pearson_hc_lc_n\tfamily_pearson_hc_lc_r\tfamily_pearson_hc_lc_p\thclc_pearson_abundance_hc_n\thclc_pearson_abundance_hc_r\thclc_pearson_abundance_hc_p\thclc_pearson_abundance_lc_n\thclc_pearson_abundance_lc_r\thclc_pearson_abundance_lc_p\thclc_pearson_abundance_paired_n\thclc_pearson_abundance_paired_r\thclc_pearson_abundance_paired_p\thclc_pearson_hc_lc_n\thclc_pearson_hc_lc_r\thclc_pearson_hc_lc_p\tmst_rooted_isotype_svg_hex\tmst_rooted_light_chain_svg_hex\tmst_rooted_hc_depth_svg_hex\tmst_rooted_lc_depth_svg_hex\tmst_rooted_paired_depth_svg_hex"
    )?;
    for family in families {
        let a = family_plot_analysis(family);
        let hc_by_cell: HashMap<&str, &AlignedCell> = family
            .members
            .iter()
            .map(|m| (m.cell.cell_id.as_str(), m))
            .collect();
        let mut emitted = BTreeSet::new();
        for lc in &family.light_clones {
            let hs = a.hclc.get(&lc.name).copied().unwrap_or_default();
            for lm in &lc.members {
                let Some(hm) = hc_by_cell.get(lm.cell.as_str()) else {
                    continue;
                };
                emitted.insert(lm.cell.clone());
                let hc_count = hm.hc_mutations.mutation_events();
                let lc_count = lm.mutations.mutation_events();
                let hc_dist = mutation_distance(&hm.hc_mutations);
                let lc_dist = mutation_distance(&lm.mutations);
                let (source, cell) = split_scoped_cell_id(&lm.cell);
                let bd_cell_id = bd_cell_id_for_seq(cell, cell_id_detectors);
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    source,
                    cell,
                    bd_cell_id,
                    family.name,
                    lc.name,
                    "paired",
                    hm.cell.hc.id,
                    lm.receptor.id,
                    hm.cell.hc.c,
                    hc_count,
                    lc_count,
                    hc_count + lc_count,
                    hc_dist,
                    lc_dist,
                    hc_dist + lc_dist,
                    a.abundance_hc.n,
                    fmt_stat(a.abundance_hc.r),
                    fmt_stat(a.abundance_hc.p),
                    a.abundance_lc.n,
                    fmt_stat(a.abundance_lc.r),
                    fmt_stat(a.abundance_lc.p),
                    a.abundance_paired.n,
                    fmt_stat(a.abundance_paired.r),
                    fmt_stat(a.abundance_paired.p),
                    a.hc_lc.n,
                    fmt_stat(a.hc_lc.r),
                    fmt_stat(a.hc_lc.p),
                    hs.abundance_hc.n,
                    fmt_stat(hs.abundance_hc.r),
                    fmt_stat(hs.abundance_hc.p),
                    hs.abundance_lc.n,
                    fmt_stat(hs.abundance_lc.r),
                    fmt_stat(hs.abundance_lc.p),
                    hs.abundance_paired.n,
                    fmt_stat(hs.abundance_paired.r),
                    fmt_stat(hs.abundance_paired.p),
                    hs.hc_lc.n,
                    fmt_stat(hs.hc_lc.r),
                    fmt_stat(hs.hc_lc.p),
                    color_for(&a.isotype_hex, &lm.cell),
                    color_for(&a.light_chain_hex, &lm.cell),
                    color_for(&a.hc_depth_hex, &lm.cell),
                    color_for(&a.lc_depth_hex, &lm.cell),
                    color_for(&a.paired_depth_hex, &lm.cell)
                )?;
            }
        }
        for hm in &family.members {
            if emitted.contains(&hm.cell.cell_id) {
                continue;
            }
            let hc_count = hm.hc_mutations.mutation_events();
            let hc_dist = mutation_distance(&hm.hc_mutations);
            let (source, cell) = split_scoped_cell_id(&hm.cell.cell_id);
            let bd_cell_id = bd_cell_id_for_seq(cell, cell_id_detectors);
            let (lc_clone, membership) = match family.light_clones.as_slice() {
                [lc] => (lc.name.as_str(), "hc_only"),
                [] => ("unresolved", "hc_only_unresolved"),
                _ => ("ambiguous", "hc_only_ambiguous"),
            };
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA\t{}\t{}\tNA\t{}\t{}\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t{}\t{}\t{}\t{}\t{}",
                source,
                cell,
                bd_cell_id,
                family.name,
                lc_clone,
                membership,
                hm.cell.hc.id,
                hm.cell.hc.c,
                hc_count,
                hc_count,
                hc_dist,
                hc_dist,
                a.abundance_hc.n,
                fmt_stat(a.abundance_hc.r),
                fmt_stat(a.abundance_hc.p),
                a.abundance_lc.n,
                fmt_stat(a.abundance_lc.r),
                fmt_stat(a.abundance_lc.p),
                a.abundance_paired.n,
                fmt_stat(a.abundance_paired.r),
                fmt_stat(a.abundance_paired.p),
                a.hc_lc.n,
                fmt_stat(a.hc_lc.r),
                fmt_stat(a.hc_lc.p),
                color_for(&a.isotype_hex, &hm.cell.cell_id),
                color_for(&a.light_chain_hex, &hm.cell.cell_id),
                color_for(&a.hc_depth_hex, &hm.cell.cell_id),
                color_for(&a.lc_depth_hex, &hm.cell.cell_id),
                color_for(&a.paired_depth_hex, &hm.cell.cell_id)
            )?;
        }
    }
    for assignment in light_only {
        for observation in &assignment.cells {
            let (source, cell) = split_scoped_cell_id(&observation.cell);
            let bd_cell_id = bd_cell_id_for_seq(cell, cell_id_detectors);
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\tlc_only\tNA\t{}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA",
                source,
                cell,
                bd_cell_id,
                assignment.family,
                assignment.clone_id,
                assignment.receptor.id,
            )?;
        }
    }
    Ok(())
}

fn call_functional(row: &CallRow) -> bool {
    row.productive && receptor_complete(&row.receptor)
}

fn call_rank(
    row: &CallRow,
) -> (
    bool,
    bool,
    std::cmp::Reverse<usize>,
    std::cmp::Reverse<usize>,
    String,
) {
    (
        !call_functional(row),
        !row.productive,
        std::cmp::Reverse(row.receptor.support_umis),
        std::cmp::Reverse(row.receptor.evidence_reads),
        row.receptor.id.clone(),
    )
}

fn best_call<'a>(rows: &[&'a CallRow]) -> Option<&'a CallRow> {
    rows.iter().copied().min_by_key(|r| call_rank(r))
}

fn receptor_fields(row: Option<&CallRow>) -> Vec<String> {
    match row {
        Some(r) => vec![
            r.receptor.id.clone(),
            r.receptor.chain.clone(),
            r.receptor.v.clone(),
            r.receptor.d.clone(),
            r.receptor.j.clone(),
            r.receptor.c.clone(),
            r.receptor.cdr3_nt.clone(),
            r.productive.to_string(),
            call_functional(r).to_string(),
            r.receptor.evidence_reads.to_string(),
            r.receptor.support_umis.to_string(),
        ],
        None => vec!["NA".into(); 11],
    }
}

/// Normalize the receptor-granular family report into the public cell contract.
/// `cells.tsv` is exactly one row per scoped CellId. Competing receptor calls are
/// retained in `potential_duplicates.tsv`; nothing is silently thrown away.
fn write_unique_cells_and_duplicates(
    receptor_rows_path: &Path,
    cells_path: &Path,
    duplicates_path: &Path,
    calls: &[CallRow],
    cell_id_detectors: &[PrimerDetector],
) -> Result<()> {
    let mut calls_by_cell: BTreeMap<String, Vec<&CallRow>> = BTreeMap::new();
    for call in calls {
        calls_by_cell
            .entry(call.cell.clone())
            .or_default()
            .push(call);
    }

    let f = File::open(receptor_rows_path)?;
    let mut lines = BufReader::new(f).lines();
    let old_header = lines.next().context("empty temporary cell report")??;
    let names: Vec<&str> = old_header.split('\t').collect();
    let idx: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let source_i = *idx.get("source").context("cell report missing source")?;
    let cell_i = *idx.get("cell").context("cell report missing cell")?;
    let hc_i = *idx
        .get("hc_recombination_id")
        .context("cell report missing hc_recombination_id")?;
    let lc_i = *idx
        .get("lc_recombination_id")
        .context("cell report missing lc_recombination_id")?;
    let membership_i = *idx
        .get("clone_membership")
        .context("cell report missing clone_membership")?;

    let mut report_by_cell: BTreeMap<String, Vec<Vec<String>>> = BTreeMap::new();
    for line in lines {
        let fields: Vec<String> = line?.split('\t').map(str::to_string).collect();
        if fields.len() != names.len() {
            continue;
        }
        let scoped = scoped_cell_id(&fields[source_i], &fields[cell_i]);
        report_by_cell.entry(scoped).or_default().push(fields);
    }

    let extra_header = [
        "hc_chain",
        "hc_v_call",
        "hc_d_call",
        "hc_j_call",
        "hc_c_call_selected",
        "hc_cdr3_nt",
        "hc_productive",
        "hc_functional",
        "hc_evidence_reads",
        "hc_support_umis",
        "lc_chain",
        "lc_v_call",
        "lc_d_call",
        "lc_j_call",
        "lc_c_call_selected",
        "lc_cdr3_nt",
        "lc_productive",
        "lc_functional",
        "lc_evidence_reads",
        "lc_support_umis",
        "hc_receptors",
        "hc_functional_receptors",
        "lc_receptors",
        "lc_functional_receptors",
        "has_potential_duplicates",
    ];
    let mut out = BufWriter::new(File::create(cells_path)?);
    writeln!(out, "{}\t{}", old_header, extra_header.join("\t"))?;

    let mut dup = BufWriter::new(File::create(duplicates_path)?);
    writeln!(
        dup,
        "source\tcell\tbd_cell_id\tchain\trecombination_id\tv_call\td_call\tj_call\tc_call\tcdr3_nt\tproductive\tfunctional\tevidence_reads\tsupport_umis\tselected_in_cells_tsv"
    )?;

    for (scoped, cell_calls) in &calls_by_cell {
        let hc: Vec<&CallRow> = cell_calls
            .iter()
            .copied()
            .filter(|r| r.receptor.chain == "IGH")
            .collect();
        let lc: Vec<&CallRow> = cell_calls
            .iter()
            .copied()
            .filter(|r| matches!(r.receptor.chain.as_str(), "IGK" | "IGL"))
            .collect();
        let best_hc = best_call(&hc);
        let best_lc = best_call(&lc);
        let hc_ids: BTreeSet<&str> = hc.iter().map(|r| r.receptor.id.as_str()).collect();
        let lc_ids: BTreeSet<&str> = lc.iter().map(|r| r.receptor.id.as_str()).collect();
        let hc_func = hc
            .iter()
            .filter(|r| call_functional(r))
            .map(|r| r.receptor.id.as_str())
            .collect::<BTreeSet<_>>()
            .len();
        let lc_func = lc
            .iter()
            .filter(|r| call_functional(r))
            .map(|r| r.receptor.id.as_str())
            .collect::<BTreeSet<_>>()
            .len();
        let has_dups = hc.len() > 1 || lc.len() > 1;

        let chosen_report = report_by_cell.get(scoped).and_then(|rows| {
            let hc_id = best_hc.map(|r| r.receptor.id.as_str());
            let lc_id = best_lc
                .filter(|r| call_functional(r))
                .map(|r| r.receptor.id.as_str());
            rows.iter()
                .find(|r| {
                    hc_id.map_or(r[hc_i] == "NA", |id| r[hc_i] == id)
                        && lc_id.map_or(true, |id| r[lc_i] == id)
                })
                .or_else(|| rows.iter().find(|r| hc_id.map_or(true, |id| r[hc_i] == id)))
                .or_else(|| rows.first())
        });

        let mut base = if let Some(row) = chosen_report {
            row.clone()
        } else {
            let (source, cell) = split_scoped_cell_id(scoped);
            let bd = bd_cell_id_for_seq(cell, cell_id_detectors);
            let mut row = vec!["NA".to_string(); names.len()];
            row[source_i] = source.to_string();
            row[cell_i] = cell.to_string();
            if let Some(i) = idx.get("bd_cell_id") {
                row[*i] = bd;
            }
            row[membership_i] = "unresolved".into();
            row
        };
        // The selected receptor columns are authoritative even when family inference
        // intentionally left a cell unresolved.
        if let Some(h) = best_hc {
            base[hc_i] = h.receptor.id.clone();
        }
        if let Some(l) = best_lc {
            base[lc_i] = l.receptor.id.clone();
        }

        let mut extras = receptor_fields(best_hc);
        // receptor_fields includes id; the legacy table already owns recombination_id.
        extras.remove(0);
        let mut lfields = receptor_fields(best_lc);
        lfields.remove(0);
        extras.extend(lfields);
        extras.extend([
            hc_ids.len().to_string(),
            hc_func.to_string(),
            lc_ids.len().to_string(),
            lc_func.to_string(),
            has_dups.to_string(),
        ]);
        writeln!(out, "{}\t{}", base.join("\t"), extras.join("\t"))?;

        if has_dups {
            let (source, cell) = split_scoped_cell_id(scoped);
            let bd = bd_cell_id_for_seq(cell, cell_id_detectors);
            for r in cell_calls {
                let selected = best_hc.map_or(false, |x| std::ptr::eq(x, *r))
                    || best_lc.map_or(false, |x| std::ptr::eq(x, *r));
                writeln!(
                    dup,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    source,
                    cell,
                    bd,
                    r.receptor.chain,
                    r.receptor.id,
                    r.receptor.v,
                    r.receptor.d,
                    r.receptor.j,
                    r.receptor.c,
                    r.receptor.cdr3_nt,
                    r.productive,
                    call_functional(r),
                    r.receptor.evidence_reads,
                    r.receptor.support_umis,
                    selected
                )?;
            }
        }
    }
    Ok(())
}

fn mutation_signature(m: &MutationMeasurement) -> String {
    let substitutions = m
        .substitution_positions
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let indels = m
        .indels
        .iter()
        .map(|x| {
            format!(
                "{}{}:{}",
                if x.inserted { "+" } else { "-" },
                x.naive_pos,
                x.len
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("S[{substitutions}]I[{indels}]")
}

fn source_counts<'a>(cells: impl Iterator<Item = &'a str>) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for cell in cells {
        let (source, _) = split_scoped_cell_id(cell);
        *counts.entry(source.to_string()).or_default() += 1;
    }
    counts
}

fn counts_text(counts: &BTreeMap<String, usize>) -> String {
    counts
        .iter()
        .map(|(source, n)| format!("{source}={n}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn write_overlap_reports(out: &Path, families: &[Family]) -> Result<()> {
    let mut events: Vec<(String, String, String, String, BTreeMap<String, usize>)> = Vec::new();
    for family in families {
        let hc_counts = source_counts(family.members.iter().map(|m| m.cell.cell_id.as_str()));
        if hc_counts.len() > 1 {
            events.push((
                "HC".into(),
                family.name.clone(),
                "NA".into(),
                "NA".into(),
                hc_counts,
            ));
        }
        let hc_by_cell: HashMap<&str, &AlignedCell> = family
            .members
            .iter()
            .map(|m| (m.cell.cell_id.as_str(), m))
            .collect();
        for lc in &family.light_clones {
            let lc_counts = source_counts(lc.members.iter().map(|m| m.cell.as_str()));
            if lc_counts.len() > 1 {
                events.push((
                    "HC_LC".into(),
                    family.name.clone(),
                    lc.name.clone(),
                    "NA".into(),
                    lc_counts,
                ));
            }
            let mut signatures: BTreeMap<String, Vec<&str>> = BTreeMap::new();
            for lm in &lc.members {
                let Some(hm) = hc_by_cell.get(lm.cell.as_str()) else {
                    continue;
                };
                let signature = format!(
                    "HC:{}|LC:{}",
                    mutation_signature(&hm.hc_mutations),
                    mutation_signature(&lm.mutations)
                );
                signatures
                    .entry(signature)
                    .or_default()
                    .push(lm.cell.as_str());
            }
            for (signature, cells) in signatures {
                let counts = source_counts(cells.into_iter());
                if counts.len() > 1 {
                    events.push((
                        "HC_LC_MUTATION_SET".into(),
                        family.name.clone(),
                        lc.name.clone(),
                        signature,
                        counts,
                    ));
                }
            }
        }
    }
    if events.is_empty() {
        return Ok(());
    }

    let mut tsv = BufWriter::new(File::create(out.join("overlap_events.tsv"))?);
    writeln!(
        tsv,
        "event_type\tfamily\tlc_clone\tmutation_set\tsources\tsource_counts\ttotal_cells"
    )?;
    for (kind, family, lc, mutations, counts) in &events {
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            kind,
            family,
            lc,
            mutations,
            counts.keys().cloned().collect::<Vec<_>>().join(";"),
            counts_text(counts),
            counts.values().sum::<usize>()
        )?;
    }

    let mut md = BufWriter::new(File::create(out.join("overlap_events.md"))?);
    writeln!(md, "# Cross-source ClonoMap overlaps\n")?;
    writeln!(
        md,
        "Joint family construction found **{}** cross-source overlap events. These are derived summaries; `cells.tsv` remains the cell-level authoritative output.\n",
        events.len()
    )?;
    for (kind, family, lc, mutations, counts) in &events {
        writeln!(
            md,
            "- **{}** — `{}`{}: {} cells ({}){}",
            kind,
            family,
            if lc == "NA" {
                String::new()
            } else {
                format!(" / `{lc}`")
            },
            counts.values().sum::<usize>(),
            counts_text(counts),
            if mutations == "NA" {
                String::new()
            } else {
                format!("; mutation set `{mutations}`")
            }
        )?;
    }
    Ok(())
}

fn bd_cell_id_for_seq(cell: &str, detectors: &[PrimerDetector]) -> String {
    detectors
        .iter()
        .find_map(|detector| detector.cell_id_for_seq(cell.as_bytes()))
        .map_or_else(|| "NA".to_string(), |id| id.to_string())
}

fn fmt_stat(x: Option<f64>) -> String {
    x.map_or("NA".into(), |v| format!("{v:.6}"))
}

fn write_unassigned(path: &Path, rows: &[UnassignedCell]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "source\tcell\treceptor_id\tfrom_family")?;
    for r in rows {
        let (source, cell) = split_scoped_cell_id(&r.cell.cell_id);
        writeln!(
            w,
            "{}\t{}\t{}\t{}",
            source, cell, r.cell.hc.id, r.from_family
        )?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct CloneSequenceAudit {
    family_id: String,
    lc_clone_id: String,
    #[allow(dead_code)]
    chain: String,
    representative_recombination_id: String,
    naive_aa: String,
    mature_aa: String,
    aa_changes: usize,
}

#[derive(Debug, Clone)]
struct TrajectoryStat {
    family_id: String,
    lc_clone_id: String,
    chain: String,
    model_id: String,
    model_name: String,
    model_kind: &'static str,
    n_states: usize,
    n_cells: usize,
    rho: Option<f64>,
    nominal_p: Option<f64>,
    delta_min: f64,
    delta_max: f64,
    delta_range: f64,
}

fn average_ranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && values[order[j]] == values[order[i]] {
            j += 1;
        }
        let rank = ((i + 1 + j) as f64) / 2.0;
        for k in i..j {
            ranks[order[k]] = rank;
        }
        i = j;
    }
    ranks
}

fn spearman(x: &[f64], y: &[f64]) -> PearsonStat {
    let n = x.len().min(y.len());
    if n < 3 {
        return PearsonStat {
            n,
            ..Default::default()
        };
    }
    pearson(&average_ranks(&x[..n]), &average_ranks(&y[..n]))
}

fn summarize_trajectory_points(points: &[CloneModelPoint]) -> Vec<TrajectoryStat> {
    let mut groups: BTreeMap<
        (String, String, String, String, String, &'static str),
        Vec<&CloneModelPoint>,
    > = BTreeMap::new();
    for p in points {
        groups
            .entry((
                p.family_id.clone(),
                p.lc_clone_id.clone(),
                p.chain.clone(),
                p.model_id.clone(),
                p.model_name.clone(),
                p.model_kind,
            ))
            .or_default()
            .push(p);
    }
    let mut out = Vec::new();
    for ((family_id, lc_clone_id, chain, model_id, model_name, model_kind), rows) in groups {
        if rows.len() < 4 {
            continue;
        }
        let x: Vec<f64> = rows.iter().map(|p| p.aa_changes as f64).collect();
        let y: Vec<f64> = rows.iter().map(|p| p.delta).collect();
        let stat = spearman(&x, &y);
        let delta_min = y.iter().copied().fold(f64::INFINITY, f64::min);
        let delta_max = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        out.push(TrajectoryStat {
            family_id,
            lc_clone_id,
            chain,
            model_id,
            model_name,
            model_kind,
            n_states: rows.len(),
            n_cells: rows.iter().map(|p| p.clone_cells).max().unwrap_or(0),
            rho: stat.r,
            nominal_p: stat.p,
            delta_min,
            delta_max,
            delta_range: delta_max - delta_min,
        });
    }
    out
}

#[derive(Debug, Clone)]
struct CloneModelPoint {
    family_id: String,
    lc_clone_id: String,
    clone_cells: usize,
    clone_states: usize,
    chain: String,
    representative_recombination_id: String,
    mutation_events: usize,
    aa_changes: usize,
    model_id: String,
    model_name: String,
    model_kind: &'static str,
    naive_score: f64,
    mature_score: f64,
    delta: f64,
    precision: Option<f64>,
    recall: Option<f64>,
    f1: Option<f64>,
}

/// Number of distinct paired HC/LC mutation states already present in one
/// finalized light-chain clone.  This is deliberately a clone-level selector:
/// cell abundance does not increase the number of vault queries.
fn light_clone_state_depth(family: &Family, clone: &LightClone) -> usize {
    let hc_by_cell: HashMap<&str, usize> = family
        .members
        .iter()
        .map(|m| (m.cell.cell_id.as_str(), m.hc_mutations.mutation_events()))
        .collect();
    clone
        .members
        .iter()
        .filter_map(|m| {
            hc_by_cell
                .get(m.cell.as_str())
                .map(|hc| (*hc, m.mutations.mutation_events()))
        })
        .collect::<BTreeSet<_>>()
        .len()
}

fn score_clone_chain(
    vault: &AaModelVault,
    family_id: &str,
    lc_clone_id: &str,
    clone_cells: usize,
    clone_states: usize,
    receptor: &Receptor,
    mutations: &MutationMeasurement,
    naive_nt: &str,
) -> Result<(Vec<CloneModelPoint>, CloneSequenceAudit)> {
    // Build the inferred mature receptor on the full effective germline/naive
    // coordinate system.  Only observed substitutions replace germline bases;
    // every unobserved position therefore remains the germline expectation.
    // HC and LC are reconstructed and scored independently by separate calls.
    let mature_nt = reconstruct_substitutions(naive_nt, &receptor.observed)
        .with_context(|| format!("germline-filling mature receptor {}", receptor.id))?;
    let naive_aa = translate_imgt_alignment(naive_nt);
    let mature_aa = translate_imgt_alignment(&mature_nt);
    if naive_aa.is_empty() || mature_aa.is_empty() || naive_aa.len() != mature_aa.len() {
        bail!(
            "matched AA reconstruction failed for receptor {}: naive={} mature={}",
            receptor.id,
            naive_aa.len(),
            mature_aa.len()
        );
    }
    let aa_changes = naive_aa
        .bytes()
        .zip(mature_aa.bytes())
        .filter(|(a, b)| a != b)
        .count();
    // Strong invariant for an AA model: no AA changes means byte-identical HMM
    // input.  Catch projection/coverage bugs before they become model results.
    if aa_changes == 0 && naive_aa != mature_aa {
        bail!(
            "AA projection invariant failed for receptor {}: aa_changes=0 but naive/mature AA differ",
            receptor.id
        );
    }
    let audit = CloneSequenceAudit {
        family_id: family_id.to_owned(),
        lc_clone_id: lc_clone_id.to_owned(),
        chain: receptor.chain.clone(),
        representative_recombination_id: receptor.id.clone(),
        naive_aa: naive_aa.clone(),
        mature_aa: mature_aa.clone(),
        aa_changes,
    };
    let (naive, _) = vault.predict_sequence(naive_aa.as_bytes())?;
    let (mature, _) = vault.predict_sequence(mature_aa.as_bytes())?;
    let mature_by_id: HashMap<&str, _> = mature.iter().map(|x| (x.model_id.as_str(), x)).collect();
    let mut out = Vec::with_capacity(naive.len());
    for n in naive {
        let m = mature_by_id.get(n.model_id.as_str()).with_context(|| {
            format!(
                "vault model {} disappeared between matched naive/mature scoring",
                n.model_id
            )
        })?;
        out.push(CloneModelPoint {
            family_id: family_id.to_owned(),
            lc_clone_id: lc_clone_id.to_owned(),
            clone_cells,
            clone_states,
            chain: receptor.chain.clone(),
            representative_recombination_id: receptor.id.clone(),
            mutation_events: mutations.mutation_events(),
            aa_changes,
            model_id: n.model_id,
            model_name: n.model_name,
            model_kind: n.model_kind,
            naive_score: n.raw_score,
            mature_score: m.raw_score,
            delta: m.raw_score - n.raw_score,
            precision: n.precision,
            recall: n.recall,
            f1: n.f1,
        });
    }
    Ok((out, audit))
}

fn write_model_space(
    out: &Path,
    families: &[Family],
    references: &ReferenceModels,
    vault: &AaModelVault,
) -> Result<()> {
    create_dir_all(out).with_context(|| format!("create {}", out.display()))?;

    // Vault unit of work: ONE finalized LC clone with a long SHM history.
    // Never score cells and never score graph nodes.  A million cells supporting
    // the same clone still produce exactly one HC and one LC vault query.
    const MIN_CLONE_STATES: usize = 4;
    let family_points: Result<Vec<(Vec<CloneModelPoint>, Vec<CloneSequenceAudit>)>> = families
        .par_iter()
        .map(|family| {
            let hc_by_cell: HashMap<&str, &AlignedCell> = family
                .members
                .iter()
                .map(|m| (m.cell.cell_id.as_str(), m))
                .collect();
            let mut points = Vec::new();
            let mut audits = Vec::new();

            for lc in &family.light_clones {
                let clone_states = light_clone_state_depth(family, lc);
                if clone_states < MIN_CLONE_STATES {
                    continue;
                }
                let clone_cells = lc.members.len();
                let clone_cell_ids: BTreeSet<&str> =
                    lc.members.iter().map(|m| m.cell.as_str()).collect();

                // Pick one already-accepted, best-observed HC representative from
                // this clone.  This is selection only: no new alignment or family
                // membership decision is made here.
                if let Some(hc) = clone_cell_ids
                    .iter()
                    .filter_map(|cell| hc_by_cell.get(*cell).copied())
                    .max_by_key(|m| {
                        (
                            m.hc_mutations.informative_pairs,
                            m.hc_mutations.matching_pairs,
                            m.cell.hc.observed.len(),
                            m.hc_mutations.mutation_events(),
                        )
                    })
                {
                    let naive = references.effective_naive(&hc.cell.hc);
                    let (chain_points, audit) = score_clone_chain(
                        vault,
                        &family.root.id,
                        &lc.root.id,
                        clone_cells,
                        clone_states,
                        &hc.cell.hc,
                        &hc.hc_mutations,
                        &naive,
                    )?;
                    points.extend(chain_points);
                    audits.push(audit);
                }

                // And exactly one best-observed LC representative for the same
                // clone.  Missing LC simply means there is no LC vault query.
                if let Some(light) = lc.members.iter().max_by_key(|m| {
                    (
                        m.mutations.informative_pairs,
                        m.mutations.matching_pairs,
                        m.receptor.observed.len(),
                        m.mutations.mutation_events(),
                    )
                }) {
                    let naive = references.effective_naive(&light.receptor);
                    let (chain_points, audit) = score_clone_chain(
                        vault,
                        &family.root.id,
                        &lc.root.id,
                        clone_cells,
                        clone_states,
                        &light.receptor,
                        &light.mutations,
                        &naive,
                    )?;
                    points.extend(chain_points);
                    audits.push(audit);
                }
            }
            Ok((points, audits))
        })
        .collect();
    let family_results = family_points?;
    let mut points = Vec::new();
    let mut audits = Vec::new();
    for (family_points, family_audits) in family_results {
        points.extend(family_points);
        audits.extend(family_audits);
    }

    // Trajectory statistics are deliberately computed from distinct mutation-depth
    // states, never from cells.  The compact clone_hmm.tsv below remains the endpoint
    // audit table; these transient state scores exist only to answer "look, I found a
    // goody" questions about monotonic model-space movement.
    let trajectory_results: Result<Vec<Vec<CloneModelPoint>>> = families
        .par_iter()
        .map(|family| {
            let hc_by_cell: HashMap<&str, &AlignedCell> = family
                .members
                .iter()
                .map(|m| (m.cell.cell_id.as_str(), m))
                .collect();
            let mut out = Vec::new();
            let paired_cells: BTreeSet<&str> = family
                .light_clones
                .iter()
                .flat_map(|lc| lc.members.iter().map(|m| m.cell.as_str()))
                .collect();

            // HC-only is an independent family background, not a graph/clone member.
            let mut hc_only: BTreeMap<usize, &AlignedCell> = BTreeMap::new();
            for hc in family
                .members
                .iter()
                .filter(|m| !paired_cells.contains(m.cell.cell_id.as_str()))
            {
                let depth = hc.hc_mutations.mutation_events();
                let replace = hc_only.get(&depth).map_or(true, |old| {
                    (
                        hc.hc_mutations.informative_pairs,
                        hc.hc_mutations.matching_pairs,
                        hc.cell.hc.observed.len(),
                    ) > (
                        old.hc_mutations.informative_pairs,
                        old.hc_mutations.matching_pairs,
                        old.cell.hc.observed.len(),
                    )
                });
                if replace {
                    hc_only.insert(depth, hc);
                }
            }
            if hc_only.len() >= MIN_CLONE_STATES {
                for hc in hc_only.values() {
                    let naive = references.effective_naive(&hc.cell.hc);
                    let (mut rows, _) = score_clone_chain(
                        vault,
                        &family.root.id,
                        "HC_ONLY",
                        hc_only.len(),
                        hc_only.len(),
                        &hc.cell.hc,
                        &hc.hc_mutations,
                        &naive,
                    )?;
                    out.append(&mut rows);
                }
            }

            for lc in &family.light_clones {
                let clone_states = light_clone_state_depth(family, lc);
                if clone_states < MIN_CLONE_STATES {
                    continue;
                }
                let mut states: BTreeMap<(usize, usize), &LightMember> = BTreeMap::new();
                for light in &lc.members {
                    let Some(hc) = hc_by_cell.get(light.cell.as_str()).copied() else {
                        continue;
                    };
                    let key = (
                        hc.hc_mutations.mutation_events(),
                        light.mutations.mutation_events(),
                    );
                    let replace = states.get(&key).map_or(true, |old| {
                        (
                            light.mutations.informative_pairs,
                            light.mutations.matching_pairs,
                            light.receptor.observed.len(),
                        ) > (
                            old.mutations.informative_pairs,
                            old.mutations.matching_pairs,
                            old.receptor.observed.len(),
                        )
                    });
                    if replace {
                        states.insert(key, light);
                    }
                }
                for light in states.values() {
                    let Some(hc) = hc_by_cell.get(light.cell.as_str()).copied() else {
                        continue;
                    };
                    let naive_hc = references.effective_naive(&hc.cell.hc);
                    let (mut hc_rows, _) = score_clone_chain(
                        vault,
                        &family.root.id,
                        &lc.root.id,
                        lc.members.len(),
                        states.len(),
                        &hc.cell.hc,
                        &hc.hc_mutations,
                        &naive_hc,
                    )?;
                    out.append(&mut hc_rows);
                    let naive_lc = references.effective_naive(&light.receptor);
                    let (mut lc_rows, _) = score_clone_chain(
                        vault,
                        &family.root.id,
                        &lc.root.id,
                        lc.members.len(),
                        states.len(),
                        &light.receptor,
                        &light.mutations,
                        &naive_lc,
                    )?;
                    out.append(&mut lc_rows);
                }
            }
            Ok(out)
        })
        .collect();
    let trajectory_points: Vec<CloneModelPoint> =
        trajectory_results?.into_iter().flatten().collect();
    let mut trajectories = summarize_trajectory_points(&trajectory_points);
    trajectories.sort_by(|a, b| {
        b.rho
            .unwrap_or(0.0)
            .abs()
            .total_cmp(&a.rho.unwrap_or(0.0).abs())
            .then_with(|| b.delta_range.total_cmp(&a.delta_range))
    });

    let trajectory_path = out.join("clone_trajectories.tsv");
    let mut tw = BufWriter::new(File::create(&trajectory_path)?);
    writeln!(
        tw,
        "family_hc_recombination_id\tlc_clone_recombination_id\tchain\tmodel_id\tmodel_name\tmodel_kind\tn_states\tn_cells\tspearman_rho\tnominal_p\tdelta_min\tdelta_max\tdelta_range"
    )?;
    for t in &trajectories {
        writeln!(
            tw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.8}\t{:.8}\t{:.8}",
            t.family_id,
            t.lc_clone_id,
            t.chain,
            t.model_id,
            t.model_name,
            t.model_kind,
            t.n_states,
            t.n_cells,
            fmt_stat(t.rho),
            fmt_stat(t.nominal_p),
            t.delta_min,
            t.delta_max,
            t.delta_range
        )?;
    }

    // Pair the HC and LC responses from the *same* HC+LC clone trajectory.
    // Both correlations are therefore evaluated over the same clone states;
    // n_states is a property of that paired trajectory, not of either chain.
    let paired_path = out.join("paired_clone_trajectories.tsv");
    let mut pw = BufWriter::new(File::create(&paired_path)?);
    let mut paired_groups: BTreeMap<
        (String, String, String, String, String),
        Vec<&TrajectoryStat>,
    > = BTreeMap::new();
    for t in &trajectories {
        if t.lc_clone_id != "HC_ONLY" {
            paired_groups
                .entry((
                    t.family_id.clone(),
                    t.lc_clone_id.clone(),
                    t.model_id.clone(),
                    t.model_name.clone(),
                    t.model_kind.to_string(),
                ))
                .or_default()
                .push(t);
        }
    }
    // First materialize the paired rows.  The permutation test below treats the
    // model vault as an exploratory panel: within each biological HC+LC clone,
    // shuffle which LC model response is paired with which HC model response and
    // record the *largest* |HC rho * LC rho| obtainable anywhere in the panel.
    // The resulting maxT empirical P therefore already pays for asking hundreds
    // of models the same question.  It is a sanity metric, not a claim that the
    // winning InterPro label is an antibody domain.
    #[derive(Clone)]
    struct PairedRow {
        family: String,
        lc_clone: String,
        lc_chain: String,
        model_id: String,
        model_name: String,
        model_kind: String,
        n_states: usize,
        n_cells: usize,
        hc_rho: f64,
        hc_p: Option<f64>,
        lc_rho: f64,
        lc_p: Option<f64>,
        concordance: f64,
        hc_range: f64,
        lc_range: f64,
    }
    let mut paired_rows = Vec::<PairedRow>::new();
    for ((family, lc_clone, model_id, model_name, model_kind), rows) in paired_groups {
        let Some(hc) = rows.iter().copied().find(|t| t.chain == "IGH") else {
            continue;
        };
        let Some(lc) = rows
            .iter()
            .copied()
            .find(|t| t.chain == "IGK" || t.chain == "IGL")
        else {
            continue;
        };
        // HC and LC were scored from the same distinct paired clone states. If
        // this invariant is ever broken, do not manufacture a paired statistic.
        if hc.n_states != lc.n_states {
            continue;
        }
        let (Some(hc_rho), Some(lc_rho)) = (hc.rho, lc.rho) else {
            continue;
        };
        let concordance = hc_rho * lc_rho;
        paired_rows.push(PairedRow {
            family,
            lc_clone,
            lc_chain: lc.chain.clone(),
            model_id,
            model_name,
            model_kind,
            n_states: hc.n_states,
            n_cells: hc.n_cells.max(lc.n_cells),
            hc_rho,
            hc_p: hc.nominal_p,
            lc_rho,
            lc_p: lc.nominal_p,
            concordance,
            hc_range: hc.delta_range,
            lc_range: lc.delta_range,
        });
    }

    const N_PERMUTATIONS: usize = 1000;
    let mut null_max_by_clone: BTreeMap<(String, String), Vec<f64>> = BTreeMap::new();
    let mut paired_by_clone: BTreeMap<(String, String), Vec<usize>> = BTreeMap::new();
    for (i, r) in paired_rows.iter().enumerate() {
        paired_by_clone
            .entry((r.family.clone(), r.lc_clone.clone()))
            .or_default()
            .push(i);
    }
    // Tiny deterministic PRNG: no dependency and identical results across runs.
    let mut rng: u64 = 0x434c_4f4e_4f4d_4150;
    let mut next_u64 = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    for (key, idx) in &paired_by_clone {
        if idx.len() < 2 {
            continue;
        }
        let hc: Vec<f64> = idx.iter().map(|&i| paired_rows[i].hc_rho).collect();
        let lc: Vec<f64> = idx.iter().map(|&i| paired_rows[i].lc_rho).collect();
        let mut null = Vec::with_capacity(N_PERMUTATIONS);
        for _ in 0..N_PERMUTATIONS {
            let mut perm: Vec<usize> = (0..lc.len()).collect();
            for i in (1..perm.len()).rev() {
                let j = (next_u64() as usize) % (i + 1);
                perm.swap(i, j);
            }
            let max_abs = hc
                .iter()
                .enumerate()
                .map(|(i, &h)| (h * lc[perm[i]]).abs())
                .fold(0.0_f64, f64::max);
            null.push(max_abs);
        }
        null_max_by_clone.insert(key.clone(), null);
    }

    writeln!(
        pw,
        "family_hc_recombination_id\tlc_clone_recombination_id\tlc_chain\tmodel_id\tmodel_name\tmodel_kind\tn_states\tn_cells\thc_rho\thc_nominal_p\tlc_rho\tlc_nominal_p\tconcordance\tpermutation_maxT_p\thc_delta_range\tlc_delta_range"
    )?;
    let mut paired_goodies: Vec<(f64, String)> = Vec::new();
    let mut perm_sig = 0usize;
    for r in &paired_rows {
        let key = (r.family.clone(), r.lc_clone.clone());
        let perm_p = null_max_by_clone.get(&key).map(|null| {
            let obs = r.concordance.abs();
            (1 + null.iter().filter(|&&x| x >= obs).count()) as f64 / (null.len() + 1) as f64
        });
        if perm_p.is_some_and(|p| p <= 0.05) {
            perm_sig += 1;
        }
        writeln!(
            pw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{:.6}\t{}\t{:.8}\t{:.8}",
            r.family,
            r.lc_clone,
            r.lc_chain,
            r.model_id,
            r.model_name,
            r.model_kind,
            r.n_states,
            r.n_cells,
            r.hc_rho,
            fmt_stat(r.hc_p),
            r.lc_rho,
            fmt_stat(r.lc_p),
            r.concordance,
            fmt_stat(perm_p),
            r.hc_range,
            r.lc_range
        )?;
        let movement = r.hc_range.max(r.lc_range);
        paired_goodies.push((r.concordance.abs() * movement,
            format!("{} {} {} {}: n={} HC rho={:+.3} LC rho={:+.3} concordance={:+.3} perm-maxT p={} max delta range={:.3}",
                r.family, r.lc_clone, r.lc_chain, r.model_name, r.n_states, r.hc_rho, r.lc_rho,
                r.concordance, fmt_stat(perm_p), movement)));
    }
    paired_goodies.sort_by(|a, b| b.0.total_cmp(&a.0));

    // Family/model comparison: HC-only background versus LC-conditioned trajectories.
    let mut by_family_model: BTreeMap<
        (String, String, String, String, String),
        Vec<&TrajectoryStat>,
    > = BTreeMap::new();
    for t in &trajectories {
        by_family_model
            .entry((
                t.family_id.clone(),
                t.chain.clone(),
                t.model_id.clone(),
                t.model_name.clone(),
                t.model_kind.to_string(),
            ))
            .or_default()
            .push(t);
    }
    let summary_path = out.join("family_model_summary.tsv");
    let mut sw = BufWriter::new(File::create(&summary_path)?);
    writeln!(
        sw,
        "family_hc_recombination_id\tchain\tmodel_id\tmodel_name\tmodel_kind\thc_only_rho\tbest_positive_lc_clone\tbest_positive_rho\tbest_negative_lc_clone\tbest_negative_rho\tlc_rho_spread\tmax_abs_lc_minus_hc_only"
    )?;
    let mut goodies: Vec<(f64, String)> = Vec::new();
    for ((family, chain, model_id, model_name, model_kind), rows) in by_family_model {
        let hc_rho = rows
            .iter()
            .find(|t| t.lc_clone_id == "HC_ONLY")
            .and_then(|t| t.rho);
        let lc: Vec<_> = rows
            .iter()
            .copied()
            .filter(|t| t.lc_clone_id != "HC_ONLY" && t.rho.is_some())
            .collect();
        if lc.is_empty() {
            continue;
        }
        let pos = lc
            .iter()
            .copied()
            .max_by(|a, b| a.rho.unwrap().total_cmp(&b.rho.unwrap()));
        let neg = lc
            .iter()
            .copied()
            .min_by(|a, b| a.rho.unwrap().total_cmp(&b.rho.unwrap()));
        let spread = match (pos, neg) {
            (Some(a), Some(b)) => Some(a.rho.unwrap() - b.rho.unwrap()),
            _ => None,
        };
        let departure = hc_rho.map(|h| {
            lc.iter()
                .map(|t| (t.rho.unwrap() - h).abs())
                .fold(0.0, f64::max)
        });
        writeln!(
            sw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            family,
            chain,
            model_id,
            model_name,
            model_kind,
            fmt_stat(hc_rho),
            pos.map_or("NA", |t| t.lc_clone_id.as_str()),
            fmt_stat(pos.and_then(|t| t.rho)),
            neg.map_or("NA", |t| t.lc_clone_id.as_str()),
            fmt_stat(neg.and_then(|t| t.rho)),
            fmt_stat(spread),
            fmt_stat(departure)
        )?;
        if let Some(s) = spread {
            let range = lc.iter().map(|t| t.delta_range).fold(0.0, f64::max);
            goodies.push((
                s.abs() * range,
                format!(
                    "{} {} {}: LC rho spread={:.3}, max delta range={:.3}",
                    family, chain, model_name, s, range
                ),
            ));
        }
    }
    goodies.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!("\nInteresting model-space trajectories (>= {MIN_CLONE_STATES} states)");
    for (_, line) in goodies.iter().take(10) {
        println!("  {line}");
    }
    println!("  trajectory statistics: {}", trajectory_path.display());
    println!("  family/model ranking:  {}", summary_path.display());
    println!("  paired HC/LC trajectories: {}", paired_path.display());
    println!("\nStrongest paired HC/LC model trajectories");
    for (_, line) in paired_goodies.iter().take(10) {
        println!("  {line}");
    }
    println!(
        "  permutation sanity test: {} model rows with maxT empirical p <= 0.05 ({} permutations/clone)",
        perm_sig, N_PERMUTATIONS
    );

    // Exact protein strings handed to the vault.  This is intentionally one
    // naive/mature pair per qualifying clone chain, not one entry per model.
    let faa_path = out.join("clone_sequences.faa");
    let mut faa = BufWriter::new(File::create(&faa_path)?);
    for a in &audits {
        writeln!(
            faa,
            ">{}|{}|{}|naive|aa_changes={}",
            a.family_id, a.lc_clone_id, a.representative_recombination_id, a.aa_changes
        )?;
        writeln!(faa, "{}", a.naive_aa)?;
        writeln!(
            faa,
            ">{}|{}|{}|mature_germline_filled|aa_changes={}",
            a.family_id, a.lc_clone_id, a.representative_recombination_id, a.aa_changes
        )?;
        writeln!(faa, "{}", a.mature_aa)?;
    }

    // Deliberately only one compact persistent table.  The previous
    // responses.tsv/family_trajectories.tsv represented cell/model Cartesian
    // products and are not part of the clone-level contract.
    let path = out.join("clone_hmm.tsv");
    let mut w = BufWriter::new(File::create(&path)?);
    writeln!(
        w,
        "family_hc_recombination_id\tlc_clone_recombination_id\tclone_cells\tclone_states\tchain\trepresentative_recombination_id\tmutation_events\taa_changes\tmodel_id\tmodel_name\tmodel_kind\tnaive_raw_log_likelihood\tmature_raw_log_likelihood\tdelta_raw_log_likelihood\theldout_precision\theldout_recall\theldout_f1\tprojection"
    )?;
    for p in &points {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.8}\t{:.8}\t{:.8}\t{}\t{}\t{}\tgermline_fill_observed_substitutions",
            p.family_id,
            p.lc_clone_id,
            p.clone_cells,
            p.clone_states,
            p.chain,
            p.representative_recombination_id,
            p.mutation_events,
            p.aa_changes,
            p.model_id,
            p.model_name,
            p.model_kind,
            p.naive_score,
            p.mature_score,
            p.delta,
            fmt_stat(p.precision),
            fmt_stat(p.recall),
            fmt_stat(p.f1),
        )?;
    }
    let scored_clones = points
        .iter()
        .map(|p| (p.family_id.as_str(), p.lc_clone_id.as_str()))
        .collect::<BTreeSet<_>>()
        .len();
    println!(
        "Clone-level model-space analysis: {} qualifying LC clones; {} chain/model responses written to {}",
        scored_clones,
        points.len(),
        path.display()
    );
    Ok(())
}
