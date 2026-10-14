mod gene_model;
mod experiment;
mod regulatory;
mod module_plots;
mod cluster_profiles;
mod wt_match;
use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use hmm::{CategoricalEmission, Hmm, StateId};
use ommverse::Ommverse;
use rayon::prelude::*;
use sc_analysis::{AnalysisConfig, SingleCellData, analyze_exon_matrix};
use sc_analysis::cli::{AnalysisCli, run_analysis};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
mod ai;


/// Every analysis mode is explicit. No flag can silently select another mode.
#[derive(Debug, Parser)]
#[command(name = "traconmech", about = "Trace biological context to candidate mechanisms")]
struct App {
    #[command(subcommand)]
    mode: Mode,
}

#[derive(Debug, Subcommand)]
enum Mode {
    /// Run the existing sc-analyse CLI (all its options are forwarded).
    Analysis(AnalysisCli),
    /// TP53 -> STRING -> MSigDB regulatory analysis with per-cluster WT fallback.
    Investigate(Cli),
    /// Full investigation with portable gene-model training or fixed-model projection.
    GeneModel(GeneModelCli),
    /// Compare a case Louvain population to an adequate WT reference.
    MatchWt(Cli),
    /// Export population expression summaries and correlation tables.
    ProfileClusters(Cli),
    /// Genomic light tables and HMM analysis.
    Genomic(Cli),
    /// Evidence-based Ollama investigation.
    Ai(Cli),
    /// Inspect H5AD annotations and metadata.
    Inspect(Cli),
    /// Inspect module and gene effects in a saved investigation.
    InspectResults(InspectResultsCli),
    /// Print evidence contract registry.
    Contracts(Cli),
    /// Print the current experimental scope.
    Scope(Cli),
}

#[derive(Debug, Args)]
struct GeneModelCli {
    #[command(subcommand)]
    action: GeneModelAction,
}
#[derive(Debug, Subcommand)]
enum GeneModelAction {
    /// Investigate training data and save a reusable model.
    Create(Cli),
    /// Investigate external data and project a frozen model (no refitting).
    Apply(Cli),
}

#[derive(Debug, Args)]
struct InspectResultsCli {
    /// Existing TraConMech investigation folder.
    results: PathBuf,
    /// Only print genes from this module (default: all modules).
    #[arg(long)]
    module: Option<usize>,
    /// Maximum genes printed per module; 0 prints all.
    #[arg(long, default_value_t = 20)]
    top: usize,
}

#[derive(Debug, Args)]
struct Cli {
    /// Compare case cells with WT in the same Louvain cluster; otherwise nearest eligible cluster.
    #[arg(long)]
    match_wt_clusters: bool,
    /// Louvain cluster containing the DBA population of interest.
    #[arg(long)]
    target_cluster: Option<String>,
    /// Minimum WT cells required for a reference population.
    #[arg(long, default_value_t = 100)]
    min_wt_cells: usize,

    /// Export expression summaries and correlations for all cell populations.
    #[arg(long)]
    profile_clusters: bool,
    /// Reuse gene module assignments from an earlier effect_clusters.tsv.
    #[arg(long)]
    reference_modules: Option<PathBuf>,
    /// Seed gene for deterministic regulatory-module analysis (e.g. TP53).
    #[arg(long = "gene", visible_alias = "investigate-gene")]
    investigate_gene: Option<String>,
    /// STRING protein info table (protein ID to preferred gene symbol).
    #[arg(long)]
    string_info: Option<PathBuf>,
    /// STRING protein links table (protein1 protein2 combined_score).
    #[arg(long)]
    string_links: Option<PathBuf>,
    /// Population annotation for regulatory projection.
    #[arg(long = "cluster")]
    cluster: Option<String>,
    /// Condition annotation (e.g. sname).
    #[arg(long, default_value = "sname")]
    condition: String,
    /// Comma-separated case sample values (exact matching).
    #[arg(long)]
    case: Option<String>,
    /// Comma-separated control sample values (exact matching).
    #[arg(long)]
    control: Option<String>,
    /// Do not ask interactive questions; fail if experiment design is ambiguous.
    #[arg(long)]
    non_interactive: bool,
    /// Cluster gene targets by their shared MSigDB set membership.
    #[arg(long)]
    effect_cluster: bool,
    /// X-axis label rotation in degrees (0-85); auto for long names.
    #[arg(long)]
    label_angle: Option<f64>,
    /// Cell-count placement: horizontal above labels, or rotated with labels.
    #[arg(long, default_value = "horizontal", value_parser = ["horizontal", "rotated"])]
    count_label_layout: String,
    /// Biological cluster IDs to compare the WT/DBA effect profiles across.
    #[arg(long, num_args = 1..)]
    co_variant_clusters: Vec<String>,
    /// Minimum STRING combined score (0-1000).
    #[arg(long, default_value_t = 700)]
    string_min_score: u32,

    /// Print the versioned template, callback and evidence registry as JSON.
    #[arg(long)]
    list_contracts: bool,

    /// Ask Ollama a biological question using real expression evidence.
    /// Without --ask, --ai starts an interactive question/answer session.
    #[arg(long)]
    ask: Option<String>,

    /// Start an interactive Ollama investigation (or explicitly request it for investigate).
    #[arg(long)]
    ai: bool,

    /// Skip automatic developmental AI annotation during effect-cluster investigation.
    #[arg(long)]
    no_developmental_ai: bool,

    /// Investigation directory (conversation.md, evidence.jsonl, error.log).
    #[arg(long, default_value = "investigation")]
    out: PathBuf,

    /// Ollama model used for evidence-driven investigation.
    #[arg(long, default_value = "gpt-oss:20b")]
    model: String,

    /// Experimental context for developmental cell-type annotation.
    #[arg(long)]
    ai_primer: Option<String>,

    /// AI backend: ollama (local) or openai (cloud API).
    #[arg(long, value_parser = ["ollama", "openai"], default_value = "ollama")]
    ai_provider: String,

    /// Local Ollama API endpoint.
    #[arg(long, default_value = "http://localhost:11434")]
    ollama_url: String,

    /// Inspect an H5AD expression analysis and exit (no AI inference).
    #[arg(long)]
    inspect_h5ad: Option<PathBuf>,

    /// Grouping annotation for H5AD inspection.
    #[arg(long)]
    main_grouping: Option<String>,

    /// Print the current experimental scope and exit.
    #[arg(long)]
    scope: bool,

    /// Ommverse reference to use for the genomic projection.
    #[arg(long)]
    ommverse: Option<PathBuf>,

    /// H5AD expression file or sc-analysis directory (expression.mtx.gz, genes.tsv.gz, cells.tsv.gz).
    #[arg(long)]
    data: Option<PathBuf>,

    /// Raw Nelrune/10x exonic MEX folder. When supplied, sc-analysis performs
    /// log-normalization, PCA, over-clustering, correlation merging, and a
    /// within-cluster top-UMI pseudo-deep holdout before TraConMech runs.
    #[arg(long, conflicts_with = "data")]
    exon_matrix: Option<PathBuf>,

    /// Cell annotation used to define expression populations for pre-analysed --data.
    #[arg(long, default_value = "celltype")]
    group: String,

    /// Variable genes used for crude PCA of --exon-matrix.
    #[arg(long, default_value_t = 2000)]
    analysis_variable_genes: usize,

    /// PCA dimensions used before neighbourhood embedding.
    #[arg(long, default_value_t = 30)]
    analysis_pcs: usize,

    /// Target final cluster count (0 = automatic spatial clustering).
    #[arg(long, default_value_t = 0)]
    analysis_overclusters: usize,

    /// Merge over-clusters iteratively when mean-expression Pearson correlation reaches this value.
    #[arg(long, default_value_t = 0.95)]
    analysis_merge_pearson: f64,

    /// Output directory for normalized matrix, annotations, statistics and plots.
    #[arg(long, default_value = "sc_analysis")]
    analysis_output: PathBuf,

    /// Highest-UMI fraction held out within each merged cluster as pseudo-deep validation cells.
    #[arg(long, default_value_t = 0.10)]
    analysis_holdout_fraction: f64,

    /// Restrict cells by annotation value, e.g. --filter tech=smartseq2.
    /// May be supplied more than once; multiple filters are combined with AND.
    #[arg(long = "filter")]
    filters: Vec<String>,

    /// Detection fraction at or above which a gene is called ON in a population.
    #[arg(long, default_value_t = 0.10)]
    on_fraction: f64,

    /// Genomically ordered gene light table. Defaults to ./lights.tsv.
    #[arg(long)]
    output: Option<PathBuf>,

    /// Genomic nearest-neighbour table. Defaults to ./genomic_neighbors.tsv.
    #[arg(long)]
    neighbor_output: Option<PathBuf>,

    /// Number of nearest genes to connect to each gene on the same chromosome.
    #[arg(long, default_value_t = 5)]
    genomic_neighbors: usize,

    /// HMM state persistence used identically for real and shuffled genomic order.
    #[arg(long, default_value_t = 0.90)]
    hmm_stay_probability: f64,

    /// Baum-Welch iterations used to learn OPEN/CLOSED detection and transition probabilities.
    #[arg(long, default_value_t = 10)]
    hmm_train_iterations: usize,

    /// Per-bin HMM output. Defaults to ./hmm_bins.tsv.
    #[arg(long)]
    hmm_output: Option<PathBuf>,

    /// Hold one technology out of HMM fitting and use it as an independent
    /// dropout-recovery validation set, e.g. --hmm-holdout-tech smartseq2.
    #[arg(long)]
    hmm_holdout_tech: Option<String>,

    /// Held-out HMM validation table. Defaults to ./hmm_holdout_validation.tsv.
    #[arg(long)]
    hmm_validation_output: Option<PathBuf>,
}

#[derive(Debug, Clone)]
struct GeneLocus {
    chromosome: String,
    chromosome_id: usize,
    start: u32,
    end: u32,
}

fn main() -> Result<()> {
    let mut cli = match App::parse().mode {
        Mode::Analysis(args) => {
            return run_analysis(args);
        }
        Mode::InspectResults(args) => {
            return inspect_results(&args);
        }
        Mode::GeneModel(args) => {
            let (mut cli, create) = match args.action {
                GeneModelAction::Create(cli) => (cli, true),
                GeneModelAction::Apply(cli) => (cli, false),
            };
            anyhow::ensure!(cli.investigate_gene.is_some(), "gene-model requires --gene");
            anyhow::ensure!(cli.data.is_some() && cli.ommverse.is_some(), "gene-model requires --data and --ommverse");
            anyhow::ensure!(cli.string_info.is_some() && cli.string_links.is_some(), "gene-model requires STRING info and links");
            // In gene-model mode, --model names the portable file, not the AI backend.
            let model_path = PathBuf::from(&cli.model);
            anyhow::ensure!(cli.model != "gpt-oss:20b", "gene-model requires --model PATH");
            let input = cli.data.as_ref().unwrap();
            let data = if input.is_dir() { SingleCellData::from_mtx_dir(input)? }
                else { SingleCellData::from_h5ad(input)? };
            if !create { gene_model::validate(&model_path, cli.investigate_gene.as_deref().unwrap(), &data)?; }
            // The existing investigate path still produces the full reference plots.
            // Apply does not re-fit effect modules: frozen model projection is separate.
            if !create { cli.effect_cluster = false; }
            cli.model = "gpt-oss:20b".to_owned();
            cli.match_wt_clusters = false;
            cli.profile_clusters = false;
            let design = experiment::resolve(&data, &cli.out,
                cli.cluster.as_deref().unwrap_or(&cli.group), &cli.condition,
                cli.case.as_deref(), cli.control.as_deref(), cli.min_wt_cells,
                cli.non_interactive)?;
            regulatory::run(regulatory::Config {
                gene: cli.investigate_gene.as_deref().unwrap(),
                index: cli.ommverse.as_deref().unwrap(), data: &data, out: &cli.out,
                string_info: cli.string_info.as_deref().unwrap(),
                string_links: cli.string_links.as_deref().unwrap(),
                cluster: cli.cluster.as_deref().unwrap_or(&cli.group),
                condition: &cli.condition,
                cases: &design.cases.join(","), controls: &design.controls.join(","),
                effect_cluster: cli.effect_cluster, label_angle: cli.label_angle,
                count_label_layout: &cli.count_label_layout,
                co_variant: &cli.co_variant_clusters, min_score: cli.string_min_score,
                min_wt_cells: cli.min_wt_cells,
                ai_development: !cli.no_developmental_ai,
                ai_model: &cli.model,
                ai_url: if cli.ai_provider == "openai" { "openai" } else { &cli.ollama_url },
                ai_primer: cli.ai_primer.as_deref(),
            })?;
            experiment::record_references(&cli.out)?;
            if create {
                gene_model::create(&data, cli.cluster.as_deref().unwrap_or(&cli.group),
                    cli.investigate_gene.as_deref().unwrap(), &model_path, &cli.out)?;
            } else {
                gene_model::apply(&data, &model_path, &cli.condition, &cli.out)?;
                gene_model::plot_frozen_modules(&data, &model_path,
                    cli.cluster.as_deref().unwrap_or(&cli.group), &cli.condition,
                    &design.cases.join(","), &design.controls.join(","), &cli.out,
                    cli.label_angle, &cli.count_label_layout)?;
            }
            return Ok(());
        }
        Mode::Investigate(mut cli) => {
            // --gene is the preferred spelling; --investigate-gene remains an alias.
            anyhow::ensure!(cli.investigate_gene.is_some(),
                "investigate requires --gene TP53");
            anyhow::ensure!(cli.data.is_some(), "investigate requires --data");
            anyhow::ensure!(cli.ommverse.is_some(), "investigate requires --ommverse");
            anyhow::ensure!(cli.string_info.is_some() && cli.string_links.is_some(),
                "investigate requires --string-info and --string-links");
            cli.match_wt_clusters = false;
            cli.profile_clusters = false;
            cli
        }
        Mode::MatchWt(mut cli) => {
            anyhow::ensure!(cli.data.is_some(), "match-wt requires --data");
            anyhow::ensure!(cli.target_cluster.is_some(), "match-wt requires --target-cluster");
            anyhow::ensure!(cli.case.is_some() && cli.control.is_some(),
                "match-wt requires --case and --control");
            cli.match_wt_clusters = true;
            cli.profile_clusters = false;
            cli.investigate_gene = None;
            cli
        }
        Mode::ProfileClusters(mut cli) => {
            anyhow::ensure!(cli.data.is_some(), "profile-clusters requires --data");
            cli.profile_clusters = true;
            cli.match_wt_clusters = false;
            cli.investigate_gene = None;
            cli
        }
        Mode::Genomic(mut cli) => {
            anyhow::ensure!(cli.data.is_some() || cli.exon_matrix.is_some(),
                "genomic requires --data or --exon-matrix");
            cli.investigate_gene = None;
            cli.match_wt_clusters = false;
            cli.profile_clusters = false;
            cli.ai = false;
            cli.ask = None;
            cli
        }
        Mode::Ai(mut cli) => {
            anyhow::ensure!(cli.data.is_some(), "ai requires --data");
            cli.investigate_gene = None;
            cli.match_wt_clusters = false;
            cli.profile_clusters = false;
            cli.ai = true;
            cli
        }
        Mode::Inspect(mut cli) => {
            anyhow::ensure!(cli.inspect_h5ad.is_some(), "inspect requires --inspect-h5ad");
            cli
        }
        Mode::Contracts(mut cli) => { cli.list_contracts = true; cli }
        Mode::Scope(mut cli) => { cli.scope = true; cli }
    };
    if cli.list_contracts {
        println!("{}", serde_json::to_string_pretty(&traconmech::contracts::registry_json()?)?);
        return Ok(());
    }
    println!("{} — {}", traconmech::NAME, traconmech::TAGLINE);

    if let Some(path) = &cli.inspect_h5ad {
        let (overview, warnings) = traconmech::exprs::inspect_overview(path)?;
        println!("Dataset: {} cells x {} features", overview.cells, overview.features);
        for (label, columns) in [("obs (cell annotations)", &overview.obs),
                                 ("var (feature annotations)", &overview.var)] {
            println!("\n{label}:");
            for col in columns {
                let distinct = col.distinct_values.map(|n| format!(", {n} distinct"))
                    .unwrap_or_default();
                println!("  {}: {}{}", col.name, col.data_type, distinct);
                if let Some(counts) = &col.counts {
                    for (value, count) in counts { println!("    {value}\t{count}"); }
                }
            }
        }
        if let Some(grouping) = &cli.main_grouping {
            let summary = traconmech::exprs::inspect(path, grouping)?;
            println!("\nSelected grouping: {grouping}");
            for (group, count) in summary.groups { println!("{group}\t{count}"); }
        }
        for warning in warnings { eprintln!("H5AD warning: {warning}"); }
        return Ok(());
    }
    if cli.scope {
        println!(
            "Current experiment: test whether genomically neighbouring genes share expression states across populations."
        );
        println!(
            "The neighbour graph uses genomic position only and never connects genes across chromosomes."
        );
    }

    if cli.match_wt_clusters {
        let path = cli.data.as_ref().context("--match-wt-clusters requires --data")?;
        let data = if path.is_dir() { SingleCellData::from_mtx_dir(path)? }
            else { SingleCellData::from_h5ad(path)? };
        return wt_match::run(&data, &cli.out,
            cli.cluster.as_deref().unwrap_or(&cli.group), &cli.condition,
            cli.case.as_deref().context("--match-wt-clusters requires --case")?,
            cli.control.as_deref().context("--match-wt-clusters requires --control")?,
            cli.target_cluster.as_deref().context("--match-wt-clusters requires --target-cluster")?,
            cli.min_wt_cells);
    }
    if cli.profile_clusters {
        let path = cli.data.as_ref().context("--profile-clusters requires --data")?;
        let data = if path.is_dir() { SingleCellData::from_mtx_dir(path)? }
            else if path.extension().is_some_and(|e| e == "h5ad") { SingleCellData::from_h5ad(path)? }
            else { anyhow::bail!("--data must be H5AD or sc-analysis directory"); };
        return cluster_profiles::run(&data, &cli.out,
            cli.cluster.as_deref().unwrap_or(&cli.group),
            if cli.case.is_some() || cli.control.is_some() { Some(&cli.condition) } else { None },
            cli.case.as_deref().unwrap_or(""),
            cli.control.as_deref().unwrap_or(""),
            cli.reference_modules.as_deref());
    }
    let Some(ommverse_path) = cli.ommverse.as_ref() else {
        if cli.scope && cli.data.is_none() && cli.exon_matrix.is_none() {
            return Ok(());
        }
        anyhow::bail!("--ommverse is required with --data or --exon-matrix");
    };
    if let Some(gene) = &cli.investigate_gene {
        let path = cli.data.as_ref().context("--investigate-gene requires --data")?;
        let data = if path.is_dir() { SingleCellData::from_mtx_dir(path)? }
            else { SingleCellData::from_h5ad(path)? };
        let design = experiment::resolve(&data, &cli.out,
            cli.cluster.as_deref().unwrap_or(&cli.group), &cli.condition,
            cli.case.as_deref(), cli.control.as_deref(), cli.min_wt_cells,
            cli.non_interactive)?;
        let result = regulatory::run(regulatory::Config {
            gene, index: ommverse_path, data: &data, out: &cli.out,
            string_info: cli.string_info.as_deref().context("--string-info required")?,
            string_links: cli.string_links.as_deref().context("--string-links required")?,
            cluster: cli.cluster.as_deref().unwrap_or(&cli.group),
            condition: &cli.condition,
            cases: &design.cases.join(","),
            controls: &design.controls.join(","),
            effect_cluster: cli.effect_cluster,
            label_angle: cli.label_angle,
            count_label_layout: &cli.count_label_layout,
            co_variant: &cli.co_variant_clusters,
            min_score: cli.string_min_score,
            min_wt_cells: cli.min_wt_cells,
            ai_development: !cli.no_developmental_ai, ai_model: &cli.model, ai_url: if cli.ai_provider == "openai" { "openai" } else { &cli.ollama_url }, ai_primer: cli.ai_primer.as_deref(),
        });
        result?;
        experiment::record_references(&cli.out)?;
        return Ok(());
    }
    if let Some(data_path) = cli.data.as_ref() {
        if cli.ai || cli.ask.is_some() {
            let data = if data_path.is_dir() {
                SingleCellData::from_mtx_dir(data_path)
            } else if data_path.extension().is_some_and(|ext| ext == "h5ad") {
                SingleCellData::from_h5ad(data_path)
            } else {
                anyhow::bail!("--data must be an H5AD file or sc-analysis MEX directory");
            }.with_context(|| format!("loading {}", data_path.display()))?;
            // Validate the reference path here; expression tools do not yet use its contents.
            anyhow::ensure!(ommverse_path.exists(), "Ommverse reference not found: {}", ommverse_path.display());
            return ai::investigate(&data, &cli.group, &cli.out, &cli.model, if cli.ai_provider == "openai" { "openai" } else { &cli.ollama_url }, cli.ask.as_deref());
        }

        let data = if data_path.is_dir() {
            SingleCellData::from_mtx_dir(data_path)
        } else if data_path.extension().is_some_and(|ext| ext == "h5ad") {
            SingleCellData::from_h5ad(data_path)
        } else {
            anyhow::bail!(
                "--data must be an H5AD file or sc-analysis MEX directory: {}",
                data_path.display()
            );
        }
            .with_context(|| format!("loading single-cell data {}", data_path.display()))?;

        run_lights(
            &cli,
            ommverse_path,
            data,
            &cli.group,
            cli.hmm_holdout_tech.as_deref(),
        )
    } else if let Some(exon_path) = cli.exon_matrix.as_ref() {
        anyhow::ensure!(!cli.ai && cli.ask.is_none(), "AI investigation currently requires --data");
        if cli.hmm_holdout_tech.is_some() {
            anyhow::bail!("--hmm-holdout-tech is automatic with --exon-matrix (pseudo_deep)");
        }
        if !cli.filters.is_empty() {
            anyhow::bail!("--filter is not supported with raw --exon-matrix analysis");
        }
        let config = AnalysisConfig {
            variable_genes: cli.analysis_variable_genes,
            pca_components: cli.analysis_pcs,
            umap_neighbors: 15,
            exclude_vdj_before_normalization: false,
            min_umi_count: Some(400),
            umap_epochs: 200,
            target_clusters: (cli.analysis_overclusters > 0).then_some(cli.analysis_overclusters),
            merge_pearson: cli.analysis_merge_pearson,
            holdout_fraction: cli.analysis_holdout_fraction,
            pseudo_samples: 5,
            output_dir: cli.analysis_output.clone(),
            ..AnalysisConfig::default()
        };
        println!(
            "sc-analysis: crude raw-matrix analysis of {}",
            exon_path.display()
        );
        println!(
            "  QC: mitochondrial/ribosomal counts retained as covariates and excluded before normalization"
        );
        println!(
            "  normalization: log1p(10,000 * count / surviving cell total); PCA={} dims; spatial clustering target={}",
            config.pca_components, config.target_clusters.map_or("auto".to_string(), |n| n.to_string())
        );
        println!(
            "  iterative mean-expression Pearson merge >= {:.3}; top {:.1}% surviving UMI/cluster held out",
            config.merge_pearson,
            100.0 * config.holdout_fraction
        );
        println!("  analysis output: {}", config.output_dir.display());
        let (data, summary) = analyze_exon_matrix(exon_path, &config)?;
        println!(
            "  matrix: {} genes x {} cells; variable genes: {}; over-clusters: {}; merged clusters: {}; pseudo-deep holdout: {} cells",
            summary.input_features,
            summary.input_cells,
            summary.variable_genes,
            summary.initial_clusters,
            summary.merged_clusters,
            summary.holdout_cells
        );
        run_lights(&cli, ommverse_path, data, "cluster", Some("pseudo_deep"))
    } else {
        anyhow::bail!("supply either --data or --exon-matrix");
    }
}

fn run_lights(
    cli: &Cli,
    ommverse_path: &Path,
    data: SingleCellData,
    group_column: &str,
    holdout_tech: Option<&str>,
) -> Result<()> {
    let model = Ommverse::load(ommverse_path)
        .with_context(|| format!("loading Ommverse {}", ommverse_path.display()))?;
    let mut groups = data.annotations.groups(group_column)?;

    let mut selected = vec![true; data.n_cells()];
    for filter in &cli.filters {
        let (column, wanted) = filter
            .split_once('=')
            .with_context(|| format!("invalid --filter {filter:?}; expected COLUMN=VALUE"))?;
        if column.is_empty() || wanted.is_empty() {
            anyhow::bail!("invalid --filter {filter:?}; expected non-empty COLUMN=VALUE");
        }
        let values = data
            .annotations
            .column(column)
            .with_context(|| format!("cell annotation column {column:?} does not exist"))?;
        for (cell_idx, value) in values.iter().enumerate() {
            selected[cell_idx] &= value == wanted;
        }
    }
    if !cli.filters.is_empty() {
        for cells in groups.values_mut() {
            cells.retain(|&cell_idx| selected[cell_idx]);
        }
        groups.retain(|_, cells| !cells.is_empty());
        let retained = selected.iter().filter(|&&keep| keep).count();
        println!("cell filter: {}", cli.filters.join(" AND "));
        println!("cells retained: {retained} / {}", data.n_cells());
        if retained == 0 {
            anyhow::bail!("cell filter retained no cells");
        }
    }

    println!("assembly: {}", model.assembly);
    println!("features: {}", data.n_features());
    println!("cells: {}", data.n_cells());
    println!("groups from {:?}: {}", group_column, groups.len());
    for (name, cells) in &groups {
        println!("  {name}: {}", cells.len());
    }

    let mut mapped = Vec::<(usize, GeneLocus)>::new();
    let mut matched_without_locus = Vec::<String>::new();
    let mut unmapped = Vec::<String>::new();

    for (feature_idx, feature) in data.features.iter().enumerate() {
        match model.gene(feature) {
            Some(gene) => match gene_locus(&model, gene.id) {
                Some(locus) => mapped.push((feature_idx, locus)),
                None => matched_without_locus.push(feature.clone()),
            },
            None => unmapped.push(feature.clone()),
        }
    }
    mapped.sort_by_key(|(_, locus)| (locus.chromosome_id, locus.start, locus.end));

    let splice_aliases = model
        .splice
        .genes
        .iter()
        .map(|g| g.names.len())
        .sum::<usize>();
    println!("Ommverse splice genes: {}", model.splice.genes.len());
    println!("Ommverse splice gene names/aliases: {splice_aliases}");
    println!("Ommverse protein records: {}", model.proteins.len());
    println!(
        "gene identifiers resolved with loci: {}{}",
        mapped.len(),
        mapped_example_suffix(&mapped, &data.features, 20)
    );
    println!(
        "matched but without genomic locus: {}{}",
        matched_without_locus.len(),
        example_suffix(&matched_without_locus, 12)
    );
    println!(
        "unmapped expression features: {}{}",
        unmapped.len(),
        example_suffix(&unmapped, 20)
    );
    println!(
        "genes mapped to Ommverse: {} / {} ({:.2}%)",
        mapped.len(),
        data.n_features(),
        100.0 * mapped.len() as f64 / data.n_features() as f64
    );

    let group_names = groups.keys().cloned().collect::<Vec<_>>();
    let mut cell_group = vec![usize::MAX; data.n_cells()];
    let mut group_sizes = vec![0usize; group_names.len()];
    for (group_idx, name) in group_names.iter().enumerate() {
        for &cell_idx in &groups[name] {
            cell_group[cell_idx] = group_idx;
        }
        group_sizes[group_idx] = groups[name].len();
    }

    let output = cli
        .output
        .clone()
        .unwrap_or_else(|| PathBuf::from("lights.tsv"));
    let mut writer = BufWriter::new(File::create(&output)?);
    write!(writer, "chromosome\tstart\tend\tgene")?;
    for name in &group_names {
        write!(writer, "\t{name}")?;
    }
    writeln!(writer)?;

    let mut on_counts = vec![0usize; group_names.len()];
    let mut gene_fractions = Vec::<Vec<f64>>::with_capacity(mapped.len());
    for (feature_idx, locus) in &mapped {
        let mut detected = vec![0usize; group_names.len()];
        if let Some(row) = data.matrix.outer_view(*feature_idx) {
            for (cell_idx, value) in row.iter() {
                let group_idx = cell_group[cell_idx];
                if *value != 0.0 && group_idx != usize::MAX {
                    detected[group_idx] += 1;
                }
            }
        }
        write!(
            writer,
            "{}\t{}\t{}\t{}",
            locus.chromosome, locus.start, locus.end, data.features[*feature_idx]
        )?;
        let mut fractions = Vec::with_capacity(group_names.len());
        for group_idx in 0..group_names.len() {
            let fraction = if group_sizes[group_idx] == 0 {
                0.0
            } else {
                detected[group_idx] as f64 / group_sizes[group_idx] as f64
            };
            if fraction >= cli.on_fraction {
                on_counts[group_idx] += 1;
            }
            fractions.push(fraction);
            write!(writer, "\t{fraction:.6}")?;
        }
        gene_fractions.push(fractions);
        writeln!(writer)?;
    }
    writer.flush()?;

    println!("ON threshold: detection fraction >= {:.3}", cli.on_fraction);
    for (idx, name) in group_names.iter().enumerate() {
        println!("  {name}: {} mapped genes ON", on_counts[idx]);
    }
    println!(
        "wrote genomically ordered light table: {}",
        output.display()
    );

    if cli.genomic_neighbors == 0 {
        anyhow::bail!("--genomic-neighbors must be >= 1");
    }

    let neighbor_output = cli
        .neighbor_output
        .clone()
        .unwrap_or_else(|| PathBuf::from("genomic_neighbors.tsv"));
    let mut neighbor_writer = BufWriter::new(File::create(&neighbor_output)?);
    writeln!(
        neighbor_writer,
        "chromosome\tgene\tneighbor\trank\tdistance_bp\tlight_agreement\tdetection_pearson"
    )?;

    let mut edges = 0usize;
    let distance_bins = vec![
        ("<10kb", 0u64, 10_000u64),
        ("10-25kb", 10_000, 25_000),
        ("25-50kb", 25_000, 50_000),
        ("50-100kb", 50_000, 100_000),
        ("100-250kb", 100_000, 250_000),
        ("250kb-1Mb", 250_000, 1_000_000),
        (">=1Mb", 1_000_000, u64::MAX),
    ];
    let mut bin_stats = vec![(0usize, 0.0f64, 0.0f64, 0usize); distance_bins.len()];
    let mut rank_stats = vec![(0usize, 0.0f64, 0.0f64, 0usize); cli.genomic_neighbors];

    let mut chr_begin = 0usize;
    while chr_begin < mapped.len() {
        let chr_id = mapped[chr_begin].1.chromosome_id;
        let mut chr_end = chr_begin + 1;
        while chr_end < mapped.len() && mapped[chr_end].1.chromosome_id == chr_id {
            chr_end += 1;
        }

        for i in chr_begin..chr_end {
            let center = gene_midpoint(&mapped[i].1);
            let mut candidates = Vec::<(u64, usize)>::new();
            for j in chr_begin..chr_end {
                if i == j {
                    continue;
                }
                let distance = center.abs_diff(gene_midpoint(&mapped[j].1));
                candidates.push((distance, j));
            }
            candidates.sort_unstable_by_key(|&(distance, j)| (distance, j));

            for (rank0, &(distance, j)) in candidates.iter().take(cli.genomic_neighbors).enumerate()
            {
                let rank = rank0 + 1;
                let agreement =
                    light_agreement(&gene_fractions[i], &gene_fractions[j], cli.on_fraction);
                let pearson = pearson(&gene_fractions[i], &gene_fractions[j]);
                let feature_i = mapped[i].0;
                let feature_j = mapped[j].0;
                writeln!(
                    neighbor_writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}",
                    mapped[i].1.chromosome,
                    data.features[feature_i],
                    data.features[feature_j],
                    rank,
                    distance,
                    agreement,
                    pearson
                        .map(|x| format!("{x:.6}"))
                        .unwrap_or_else(|| "NA".to_string())
                )?;
                edges += 1;

                let r = &mut rank_stats[rank0];
                r.0 += 1;
                r.1 += agreement;
                if let Some(value) = pearson {
                    r.2 += value;
                    r.3 += 1;
                }

                if let Some(bin_idx) = distance_bins
                    .iter()
                    .position(|&(_, lo, hi)| distance >= lo && distance < hi)
                {
                    let b = &mut bin_stats[bin_idx];
                    b.0 += 1;
                    b.1 += agreement;
                    if let Some(value) = pearson {
                        b.2 += value;
                        b.3 += 1;
                    }
                }
            }
        }
        chr_begin = chr_end;
    }
    neighbor_writer.flush()?;

    println!(
        "genomic neighbour graph: {} directed edges (up to {} nearest genes per gene; same chromosome only)",
        edges, cli.genomic_neighbors
    );
    println!("by neighbour rank:");
    for (rank0, &(n, agreement_sum, pearson_sum, pearson_n)) in rank_stats.iter().enumerate() {
        if n == 0 {
            continue;
        }
        let mean_pearson = if pearson_n == 0 {
            "NA".to_string()
        } else {
            format!("{:.4}", pearson_sum / pearson_n as f64)
        };
        println!(
            "  rank {}: n={}, mean light agreement={:.4}, mean detection Pearson={}",
            rank0 + 1,
            n,
            agreement_sum / n as f64,
            mean_pearson
        );
    }
    println!("by genomic distance:");
    for (bin_idx, &(label, _, _)) in distance_bins.iter().enumerate() {
        let (n, agreement_sum, pearson_sum, pearson_n) = bin_stats[bin_idx];
        if n == 0 {
            println!("  {label}: n=0");
            continue;
        }
        let mean_pearson = if pearson_n == 0 {
            "NA".to_string()
        } else {
            format!("{:.4}", pearson_sum / pearson_n as f64)
        };
        println!(
            "  {label}: n={n}, mean light agreement={:.4}, mean detection Pearson={mean_pearson}",
            agreement_sum / n as f64
        );
    }
    println!(
        "wrote genomic nearest-neighbour table: {}",
        neighbor_output.display()
    );

    run_splice_bin_hmm(
        cli,
        &model,
        &data,
        &mapped,
        &gene_fractions,
        &group_names,
        &groups,
        holdout_tech,
    )?;

    Ok(())
}

fn gene_midpoint(locus: &GeneLocus) -> u64 {
    locus.start as u64 + (locus.end as u64 - locus.start as u64) / 2
}

fn light_agreement(a: &[f64], b: &[f64], on_fraction: f64) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    let same = a
        .iter()
        .zip(b)
        .filter(|(x, y)| (**x >= on_fraction) == (**y >= on_fraction))
        .count();
    same as f64 / a.len() as f64
}

fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.len() < 2 {
        return None;
    }
    let n = a.len() as f64;
    let mean_a = a.iter().sum::<f64>() / n;
    let mean_b = b.iter().sum::<f64>() / n;
    let mut numerator = 0.0;
    let mut sum_sq_a = 0.0;
    let mut sum_sq_b = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        let da = x - mean_a;
        let db = y - mean_b;
        numerator += da * db;
        sum_sq_a += da * da;
        sum_sq_b += db * db;
    }
    let denominator = (sum_sq_a * sum_sq_b).sqrt();
    (denominator > 0.0).then_some(numerator / denominator)
}

fn example_suffix(values: &[String], limit: usize) -> String {
    if values.is_empty() {
        return String::new();
    }
    let examples = values
        .iter()
        .take(limit)
        .cloned()
        .collect::<Vec<_>>()
        .join(", ");
    format!(" [e.g. {examples}]")
}

fn mapped_example_suffix(
    mapped: &[(usize, GeneLocus)],
    features: &[String],
    limit: usize,
) -> String {
    if mapped.is_empty() {
        return String::new();
    }
    let mut examples = BTreeSet::new();
    for (feature_idx, _) in mapped {
        examples.insert(features[*feature_idx].clone());
        if examples.len() >= limit {
            break;
        }
    }
    format!(
        " [e.g. {}]",
        examples.into_iter().collect::<Vec<_>>().join(", ")
    )
}

fn gene_locus(model: &Ommverse, gene_id: usize) -> Option<GeneLocus> {
    let gene = model.splice.genes.get(gene_id)?;
    let first_tx = *gene.transcript_ids().first()?;
    let chromosome_id = model.splice.transcripts.get(first_tx)?.chr_id;
    let chromosome = model.splice.chr_names.get(chromosome_id)?.clone();
    let mut start = u32::MAX;
    let mut end = 0u32;
    for &tx_id in gene.transcript_ids() {
        let tx = model.splice.transcripts.get(tx_id)?;
        if tx.chr_id != chromosome_id {
            continue;
        }
        for exon in tx.exons() {
            start = start.min(exon.start);
            end = end.max(exon.end);
        }
    }
    (start < end).then_some(GeneLocus {
        chromosome,
        chromosome_id,
        start,
        end,
    })
}

#[derive(Debug, Clone)]
struct HmmSegment {
    chromosome: String,
    start: u32,
    end: u32,
    name: String,
    feature_idx: Option<usize>,
    detected: Vec<usize>,
}

#[derive(Debug)]
struct HmmGroupResult {
    group_idx: usize,
    real_ll: f64,
    shuffled_ll: f64,
    observations: usize,
    detected: usize,
    real_transitions: usize,
    shuffled_transitions: usize,
    open_detection_probability: f64,
    closed_detection_probability: f64,
    open_stay_probability: f64,
    closed_stay_probability: f64,
    reference_gap_bp: f64,
    rows: Vec<(String, u32, u32, usize, usize, usize, f64)>,
}

fn run_splice_bin_hmm(
    cli: &Cli,
    model: &Ommverse,
    data: &SingleCellData,
    mapped: &[(usize, GeneLocus)],
    _gene_fractions: &[Vec<f64>],
    group_names: &[String],
    groups: &BTreeMap<String, Vec<usize>>,
    holdout_tech: Option<&str>,
) -> Result<()> {
    if !(0.5..1.0).contains(&cli.hmm_stay_probability) {
        anyhow::bail!("--hmm-stay-probability must be >= 0.5 and < 1.0");
    }
    if cli.hmm_train_iterations == 0 {
        anyhow::bail!("--hmm-train-iterations must be >= 1");
    }

    // Ommverse defines the HMM gene universe. The expression matrix only supplies
    // evidence onto that fixed genomic backbone. A reference gene that has no
    // matching panc8 feature therefore contributes observation 0: panc8 contains
    // no evidence that it was expressed. This is intentionally stricter than the
    // light/neighbour analyses above, which remain limited to mapped matrix rows.
    let mut feature_for_gene = HashMap::<usize, usize>::new();
    for (feature_idx, _) in mapped {
        if let Some(gene) = model.gene(&data.features[*feature_idx]) {
            feature_for_gene.entry(gene.id).or_insert(*feature_idx);
        }
    }

    let holdout = if let Some(holdout_tech) = holdout_tech {
        if cli.filters.iter().any(|filter| {
            filter
                .split_once('=')
                .is_some_and(|(column, _)| column == "tech")
        }) {
            anyhow::bail!(
                "--hmm-holdout-tech cannot be combined with --filter tech=...; the holdout must remain available for validation"
            );
        }
        let tech = data
            .annotations
            .column("tech")
            .context("--hmm-holdout-tech requires a cell annotation column named \"tech\"")?;
        let mut training_cells = vec![Vec::<usize>::new(); group_names.len()];
        let mut holdout_cells = vec![Vec::<usize>::new(); group_names.len()];
        for (group_idx, name) in group_names.iter().enumerate() {
            for &cell_idx in &groups[name] {
                if tech[cell_idx] == holdout_tech {
                    holdout_cells[group_idx].push(cell_idx);
                } else {
                    training_cells[group_idx].push(cell_idx);
                }
            }
        }
        let mut training_group_for_cell = vec![usize::MAX; data.n_cells()];
        let mut holdout_group_for_cell = vec![usize::MAX; data.n_cells()];
        for group_idx in 0..group_names.len() {
            for &cell_idx in &training_cells[group_idx] {
                training_group_for_cell[cell_idx] = group_idx;
            }
            for &cell_idx in &holdout_cells[group_idx] {
                holdout_group_for_cell[cell_idx] = group_idx;
            }
        }
        Some(HmmHoldout {
            tech: holdout_tech.to_string(),
            training_cells,
            holdout_cells,
            training_group_for_cell,
            holdout_group_for_cell,
        })
    } else {
        None
    };

    let mut all_group_for_cell = vec![usize::MAX; data.n_cells()];
    for (group_idx, name) in group_names.iter().enumerate() {
        for &cell_idx in &groups[name] {
            all_group_for_cell[cell_idx] = group_idx;
        }
    }

    let mut segments = Vec::<HmmSegment>::new();
    let mut reference_without_locus = 0usize;
    let mut reference_with_matrix_feature = 0usize;
    for gene in &model.splice.genes {
        let Some(locus) = gene_locus(model, gene.id) else {
            reference_without_locus += 1;
            continue;
        };
        let feature_idx = feature_for_gene.get(&gene.id).copied();
        if feature_idx.is_some() {
            reference_with_matrix_feature += 1;
        }
        let detected = if let Some(feature_idx) = feature_idx {
            if let Some(split) = &holdout {
                detection_by_membership(
                    data,
                    feature_idx,
                    group_names.len(),
                    &split.training_group_for_cell,
                )
            } else {
                detection_by_membership(data, feature_idx, group_names.len(), &all_group_for_cell)
            }
        } else {
            vec![0usize; group_names.len()]
        };
        let name = gene
            .names
            .first()
            .cloned()
            .unwrap_or_else(|| format!("gene:{}", gene.id));
        segments.push(HmmSegment {
            chromosome: locus.chromosome,
            start: locus.start,
            end: locus.end,
            name,
            feature_idx,
            detected,
        });
    }
    segments.sort_unstable_by(|a, b| {
        let ca = model
            .splice
            .chr_names
            .iter()
            .position(|x| x == &a.chromosome)
            .unwrap_or(usize::MAX);
        let cb = model
            .splice
            .chr_names
            .iter()
            .position(|x| x == &b.chromosome)
            .unwrap_or(usize::MAX);
        (ca, a.start, a.end).cmp(&(cb, b.start, b.end))
    });
    println!(
        "HMM reference backbone: {} Ommverse genes with loci; {} have a mapped expression feature; {} reference genes without loci skipped",
        segments.len(),
        reference_with_matrix_feature,
        reference_without_locus
    );
    println!(
        "gene-order OPEN/CLOSED HMM: {} mapped genes",
        segments.len()
    );
    if let Some(split) = &holdout {
        println!(
            "HMM training excludes tech={}; held-out cells are used only after state inference",
            split.tech
        );
        for (group_idx, name) in group_names.iter().enumerate() {
            println!(
                "  {name}: {} training cells, {} held-out {} cells",
                split.training_cells[group_idx].len(),
                split.holdout_cells[group_idx].len(),
                split.tech
            );
        }
    }
    println!("HMM observation: 1 if the gene was detected in >=1 training cell; 0 otherwise");
    println!(
        "Baum-Welch training iterations: {}",
        cli.hmm_train_iterations
    );
    let stay = cli.hmm_stay_probability;
    let results = (0..group_names.len())
        .into_par_iter()
        .map(|group_idx| infer_hmm_group(group_idx, &segments, stay, cli.hmm_train_iterations))
        .collect::<Result<Vec<_>>>()?;

    if let Some(split) = &holdout {
        write_holdout_validation(cli, data, &segments, group_names, &results, split)?;
    }

    let hmm_output = cli
        .hmm_output
        .clone()
        .unwrap_or_else(|| PathBuf::from("hmm_bins.tsv"));
    let mut writer = BufWriter::new(File::create(&hmm_output)?);
    writeln!(
        writer,
        "group\tchromosome\tstart\tend\tgene\tdetected\tstate\tposterior_open"
    )?;
    for result in &results {
        for (chr, start, end, gene, detected, state, posterior) in &result.rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
                group_names[result.group_idx],
                chr,
                start,
                end,
                segments[*gene].name,
                detected,
                if *state == 1 { "OPEN" } else { "CLOSED" },
                posterior
            )?;
        }
    }
    writer.flush()?;

    println!(
        "distance-aware binary HMM, real genomic order versus deterministic within-chromosome shuffle:"
    );
    println!(
        "  detected %       genes observed in >=1 retained cell; a zero means not observed, not proven inactive"
    );
    println!(
        "  P(det|OPEN)       learned probability of observing a gene while the latent state is OPEN"
    );
    println!(
        "  P(det|CLOSED)     learned detection probability in the lower-detection latent state; this is not proof of closed chromatin"
    );
    println!("  P(O->O @ ref)     OPEN-state persistence across the reference genomic gap");
    println!("  P(C->C @ ref)     CLOSED-state persistence across the reference genomic gap");
    println!("  real/shuffle LL   log likelihood per gene; less negative is better");
    println!(
        "  delta LL          real LL/bin - shuffled LL/bin; positive values mean genomic order explains the observations better"
    );
    println!(
        "  transitions       Viterbi state changes within chromosomes; chromosome boundaries are not counted"
    );
    println!(
        "  ref gap           median positive distance between consecutive mapped genes; distance-aware inference scales transitions from this gap"
    );
    println!();

    print!("  {:<18}", "measurement");
    for result in &results {
        print!(" {:>8}", hmm_group_label(&group_names[result.group_idx]));
    }
    println!();
    print!("  {:<18}", "detected %");
    for result in &results {
        print!(
            " {:>8.1}",
            100.0 * result.detected as f64 / result.observations.max(1) as f64
        );
    }
    println!();
    print!("  {:<18}", "P(det|OPEN)");
    for result in &results {
        print!(" {:>8.3}", result.open_detection_probability);
    }
    println!();
    print!("  {:<18}", "P(det|CLOSED)");
    for result in &results {
        print!(" {:>8.3}", result.closed_detection_probability);
    }
    println!();
    print!("  {:<18}", "P(O->O @ ref)");
    for result in &results {
        print!(" {:>8.3}", result.open_stay_probability);
    }
    println!();
    print!("  {:<18}", "P(C->C @ ref)");
    for result in &results {
        print!(" {:>8.3}", result.closed_stay_probability);
    }
    println!();
    print!("  {:<18}", "real LL/bin");
    for result in &results {
        print!(
            " {:>8.4}",
            result.real_ll / result.observations.max(1) as f64
        );
    }
    println!();
    print!("  {:<18}", "shuffle LL/bin");
    for result in &results {
        print!(
            " {:>8.4}",
            result.shuffled_ll / result.observations.max(1) as f64
        );
    }
    println!();
    print!("  {:<18}", "delta LL/bin");
    for result in &results {
        let n = result.observations.max(1) as f64;
        print!(" {:>+8.4}", (result.real_ll - result.shuffled_ll) / n);
    }
    println!();
    print!("  {:<18}", "real transitions");
    for result in &results {
        print!(" {:>8}", result.real_transitions);
    }
    println!();
    print!("  {:<18}", "shuffle transitions");
    for result in &results {
        print!(" {:>8}", result.shuffled_transitions);
    }
    println!();
    print!("  {:<18}", "ref gap bp");
    for result in &results {
        print!(" {:>8.0}", result.reference_gap_bp);
    }
    println!();
    println!(
        "  labels: act_stel=activated_stellate, endo=endothelial, macro=macrophage, qui_stel=quiescent_stellate, schwann=schwann"
    );
    println!(
        "wrote gene-order OPEN/CLOSED HMM table: {}",
        hmm_output.display()
    );
    Ok(())
}

#[derive(Debug)]
struct HmmHoldout {
    tech: String,
    training_cells: Vec<Vec<usize>>,
    holdout_cells: Vec<Vec<usize>>,
    training_group_for_cell: Vec<usize>,
    holdout_group_for_cell: Vec<usize>,
}

fn detection_by_membership(
    data: &SingleCellData,
    feature_idx: usize,
    n_groups: usize,
    membership: &[usize],
) -> Vec<usize> {
    let mut detected = vec![0usize; n_groups];
    if let Some(row) = data.matrix.outer_view(feature_idx) {
        for (cell_idx, value) in row.iter() {
            let group_idx = membership[cell_idx];
            if *value != 0.0 && group_idx != usize::MAX {
                detected[group_idx] = 1;
            }
        }
    }
    detected
}

fn write_holdout_validation(
    cli: &Cli,
    data: &SingleCellData,
    segments: &[HmmSegment],
    group_names: &[String],
    results: &[HmmGroupResult],
    split: &HmmHoldout,
) -> Result<()> {
    let holdout_detection = segments
        .iter()
        .map(|seg| {
            seg.feature_idx
                .map(|feature_idx| {
                    detection_by_membership(
                        data,
                        feature_idx,
                        group_names.len(),
                        &split.holdout_group_for_cell,
                    )
                })
                .unwrap_or_else(|| vec![0usize; group_names.len()])
        })
        .collect::<Vec<_>>();
    let output = cli
        .hmm_validation_output
        .clone()
        .unwrap_or_else(|| PathBuf::from("hmm_holdout_validation.tsv"));
    let mut writer = BufWriter::new(File::create(&output)?);
    writeln!(
        writer,
        "group\ttraining_cells\tholdout_cells\ttrain_zero_open\tholdout_detected_open\trate_open\ttrain_zero_closed\tholdout_detected_closed\trate_closed\todds_ratio\tfisher_p_one_sided"
    )?;

    println!();
    println!("held-out {} dropout-recovery test", split.tech);
    println!(
        "  The HMM was fitted without held-out cells. Validation considers only genes with zero detections in training."
    );
    println!(
        "  H/L asks whether a training-zero gene inferred OPEN (H) is recovered more often in held-out cells than one inferred CLOSED (L)."
    );
    println!(
        "  OR > 1 and a higher H rate support dropout recovery; Fisher p is one-sided for enrichment in H."
    );
    println!();
    println!(
        "  {:<18} {:>7} {:>7} {:>8} {:>7} {:>8} {:>8} {:>7} {:>8} {:>8} {:>10}",
        "celltype",
        "train",
        "holdout",
        "H zero",
        "H +",
        "H rate",
        "L zero",
        "L +",
        "L rate",
        "OR",
        "Fisher p"
    );

    for result in results {
        let group_idx = result.group_idx;
        let mut h_zero = 0usize;
        let mut h_pos = 0usize;
        let mut l_zero = 0usize;
        let mut l_pos = 0usize;
        for (i, seg) in segments.iter().enumerate() {
            if seg.detected[group_idx] != 0 {
                continue;
            }
            let is_open = result.rows[i].5 == 1;
            let holdout_pos = holdout_detection[i][group_idx] != 0;
            if is_open {
                h_zero += 1;
                h_pos += usize::from(holdout_pos);
            } else {
                l_zero += 1;
                l_pos += usize::from(holdout_pos);
            }
        }
        let h_neg = h_zero.saturating_sub(h_pos);
        let l_neg = l_zero.saturating_sub(l_pos);
        let h_rate = fraction(h_pos, h_zero);
        let l_rate = fraction(l_pos, l_zero);
        let odds_ratio = ((h_pos as f64 + 0.5) * (l_neg as f64 + 0.5))
            / ((h_neg as f64 + 0.5) * (l_pos as f64 + 0.5));
        let fisher = fisher_exact_greater(h_pos, h_neg, l_pos, l_neg);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}\t{:.6e}",
            group_names[group_idx],
            split.training_cells[group_idx].len(),
            split.holdout_cells[group_idx].len(),
            h_zero,
            h_pos,
            h_rate,
            l_zero,
            l_pos,
            l_rate,
            odds_ratio,
            fisher
        )?;
        println!(
            "  {:<18} {:>7} {:>7} {:>8} {:>7} {:>7.1}% {:>8} {:>7} {:>7.1}% {:>8.2} {:>10.3e}",
            hmm_group_label(&group_names[group_idx]),
            split.training_cells[group_idx].len(),
            split.holdout_cells[group_idx].len(),
            h_zero,
            h_pos,
            100.0 * h_rate,
            l_zero,
            l_pos,
            100.0 * l_rate,
            odds_ratio,
            fisher
        );
    }
    writer.flush()?;
    println!("wrote held-out HMM validation table: {}", output.display());
    Ok(())
}

fn fraction(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn ln_choose(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    (1..=k)
        .map(|i| ((n - k + i) as f64).ln() - (i as f64).ln())
        .sum()
}

fn fisher_exact_greater(a: usize, b: usize, c: usize, d: usize) -> f64 {
    let row1 = a + b;
    let row2 = c + d;
    let col1 = a + c;
    let total = row1 + row2;
    if total == 0 {
        return 1.0;
    }
    let max_a = row1.min(col1);
    let denom = ln_choose(total, col1);
    let mut p = 0.0;
    for x in a..=max_a {
        if col1 < x || row2 < col1 - x {
            continue;
        }
        let lp = ln_choose(row1, x) + ln_choose(row2, col1 - x) - denom;
        p += lp.exp();
    }
    p.min(1.0)
}

fn hmm_group_label(name: &str) -> &str {
    match name {
        "activated_stellate" => "act_stel",
        "endothelial" => "endo",
        "macrophage" => "macro",
        "quiescent_stellate" => "qui_stel",
        other => other,
    }
}

fn infer_hmm_group(
    group_idx: usize,
    segments: &[HmmSegment],
    stay: f64,
    iterations: usize,
) -> Result<HmmGroupResult> {
    let values = segments
        .iter()
        .map(|s| s.detected[group_idx])
        .collect::<Vec<_>>();
    if values.len() < 4 {
        anyhow::bail!("too few gene observations for HMM");
    }

    // First estimate the two binary emission regimes and a reference transition
    // matrix with ordinary Baum-Welch. Distance-aware inference below then treats
    // that transition matrix as applying at the median positive inter-gene gap.
    // This preserves the existing generic HMM trainer while removing the false
    // assumption that every neighbouring gene is one equally sized genomic step.
    let mut model = Hmm::new(
        vec![0.5, 0.5],
        vec![stay, 1.0 - stay, 1.0 - stay, stay],
        vec![
            CategoricalEmission::new(vec![0.99, 0.01])?,
            CategoricalEmission::new(vec![0.20, 0.80])?,
        ],
    )?;
    model.baum_welch(&values, iterations, 1e-9)?;

    let emission0 = model.emissions()[0].probabilities();
    let emission1 = model.emissions()[1].probabilities();
    let open_state = usize::from(emission1[1] >= emission0[1]);
    let closed_state = 1 - open_state;
    let transitions = model.transition_probabilities();
    let initial = model.initial_probabilities();
    let emissions = [[emission0[0], emission0[1]], [emission1[0], emission1[1]]];

    let reference_gap_bp = median_positive_gap(segments).max(1.0);
    let real = infer_distance_aware(
        segments,
        &values,
        &initial,
        &transitions,
        &emissions,
        reference_gap_bp,
    )?;
    let shuffled_values = deterministic_chr_shuffle(segments, group_idx);
    let shuffled = infer_distance_aware(
        segments,
        &shuffled_values,
        &initial,
        &transitions,
        &emissions,
        reference_gap_bp,
    )?;
    let real_transitions = transition_count_chr(&real.viterbi, segments);
    let shuffled_transitions = transition_count_chr(&shuffled.viterbi, segments);
    let rows = segments
        .iter()
        .enumerate()
        .map(|(i, seg)| {
            let raw_state = real.viterbi[i].0;
            (
                seg.chromosome.clone(),
                seg.start,
                seg.end,
                i,
                seg.detected[group_idx],
                usize::from(raw_state == open_state),
                real.posterior[i * 2 + open_state],
            )
        })
        .collect();
    Ok(HmmGroupResult {
        group_idx,
        real_ll: real.log_likelihood,
        shuffled_ll: shuffled.log_likelihood,
        observations: values.len(),
        detected: values.iter().sum(),
        real_transitions,
        shuffled_transitions,
        open_detection_probability: model.emissions()[open_state].probabilities()[1],
        closed_detection_probability: model.emissions()[closed_state].probabilities()[1],
        open_stay_probability: transitions[open_state * 2 + open_state],
        closed_stay_probability: transitions[closed_state * 2 + closed_state],
        reference_gap_bp,
        rows,
    })
}

#[derive(Debug)]
struct DistanceHmmResult {
    log_likelihood: f64,
    viterbi: Vec<StateId>,
    posterior: Vec<f64>,
}

fn genomic_gap(left: &HmmSegment, right: &HmmSegment) -> u64 {
    if left.chromosome != right.chromosome {
        return 0;
    }
    if right.start >= left.end {
        (right.start - left.end) as u64
    } else if left.start >= right.end {
        (left.start - right.end) as u64
    } else {
        0
    }
}

fn median_positive_gap(segments: &[HmmSegment]) -> f64 {
    let mut gaps = segments
        .windows(2)
        .filter(|w| w[0].chromosome == w[1].chromosome)
        .map(|w| genomic_gap(&w[0], &w[1]))
        .filter(|&d| d > 0)
        .collect::<Vec<_>>();
    if gaps.is_empty() {
        return 1.0;
    }
    gaps.sort_unstable();
    gaps[gaps.len() / 2] as f64
}

// Convert the learned two-state transition matrix at one reference genomic gap
// into a transition matrix for an arbitrary distance. For a two-state Markov
// chain the non-trivial eigenvalue rho controls decay toward stationarity:
// P(d) = Pi + rho^(d / d_ref) * (I - Pi).
fn distance_transition(base: &[f64], distance_bp: f64, reference_gap_bp: f64) -> [f64; 4] {
    let p01 = base[1].clamp(1e-12, 1.0 - 1e-12);
    let p10 = base[2].clamp(1e-12, 1.0 - 1e-12);
    let sum = p01 + p10;
    let pi0 = p10 / sum;
    let pi1 = p01 / sum;
    let rho = (1.0 - sum).clamp(1e-12, 1.0 - 1e-12);
    let scale = (distance_bp / reference_gap_bp.max(1.0)).max(0.0);
    let r = rho.powf(scale);
    [
        pi0 + r * (1.0 - pi0),
        pi1 * (1.0 - r),
        pi0 * (1.0 - r),
        pi1 + r * (1.0 - pi1),
    ]
}

fn logsum2(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if !m.is_finite() {
        m
    } else {
        m + ((a - m).exp() + (b - m).exp()).ln()
    }
}

fn infer_distance_aware(
    segments: &[HmmSegment],
    values: &[usize],
    initial: &[f64],
    base_transition: &[f64],
    emissions: &[[f64; 2]; 2],
    reference_gap_bp: f64,
) -> Result<DistanceHmmResult> {
    if segments.len() != values.len() || segments.is_empty() {
        anyhow::bail!("invalid distance-HMM input");
    }
    let n = segments.len();
    let mut forward = vec![f64::NEG_INFINITY; n * 2];
    let mut backward = vec![0.0; n * 2];
    let mut vscore = vec![f64::NEG_INFINITY; n * 2];
    let mut backptr = vec![0usize; n * 2];
    let mut viterbi = vec![StateId(0); n];
    let mut total_ll = 0.0;

    let mut begin = 0usize;
    while begin < n {
        let chr = &segments[begin].chromosome;
        let mut end = begin + 1;
        while end < n && segments[end].chromosome == *chr {
            end += 1;
        }

        for state in 0..2 {
            let e = emissions[state][values[begin]].max(1e-300).ln();
            forward[begin * 2 + state] = initial[state].max(1e-300).ln() + e;
            vscore[begin * 2 + state] = forward[begin * 2 + state];
        }
        for i in begin + 1..end {
            let tr = distance_transition(
                base_transition,
                genomic_gap(&segments[i - 1], &segments[i]) as f64,
                reference_gap_bp,
            );
            for dst in 0..2 {
                let a = forward[(i - 1) * 2] + tr[dst].max(1e-300).ln();
                let b = forward[(i - 1) * 2 + 1] + tr[2 + dst].max(1e-300).ln();
                forward[i * 2 + dst] = logsum2(a, b) + emissions[dst][values[i]].max(1e-300).ln();
                let va = vscore[(i - 1) * 2] + tr[dst].max(1e-300).ln();
                let vb = vscore[(i - 1) * 2 + 1] + tr[2 + dst].max(1e-300).ln();
                let src = usize::from(vb > va);
                vscore[i * 2 + dst] =
                    if src == 0 { va } else { vb } + emissions[dst][values[i]].max(1e-300).ln();
                backptr[i * 2 + dst] = src;
            }
        }
        let chr_ll = logsum2(forward[(end - 1) * 2], forward[(end - 1) * 2 + 1]);
        total_ll += chr_ll;

        backward[(end - 1) * 2] = 0.0;
        backward[(end - 1) * 2 + 1] = 0.0;
        for i in (begin..end - 1).rev() {
            let tr = distance_transition(
                base_transition,
                genomic_gap(&segments[i], &segments[i + 1]) as f64,
                reference_gap_bp,
            );
            for src in 0..2 {
                let a = tr[src * 2].max(1e-300).ln()
                    + emissions[0][values[i + 1]].max(1e-300).ln()
                    + backward[(i + 1) * 2];
                let b = tr[src * 2 + 1].max(1e-300).ln()
                    + emissions[1][values[i + 1]].max(1e-300).ln()
                    + backward[(i + 1) * 2 + 1];
                backward[i * 2 + src] = logsum2(a, b);
            }
        }

        let mut state = usize::from(vscore[(end - 1) * 2 + 1] > vscore[(end - 1) * 2]);
        viterbi[end - 1] = StateId(state);
        for i in (begin + 1..end).rev() {
            state = backptr[i * 2 + state];
            viterbi[i - 1] = StateId(state);
        }

        // Normalize posterior chromosome-by-chromosome using that chromosome's LL.
        for i in begin..end {
            let p0 = (forward[i * 2] + backward[i * 2] - chr_ll).exp();
            let p1 = (forward[i * 2 + 1] + backward[i * 2 + 1] - chr_ll).exp();
            let z = (p0 + p1).max(1e-300);
            backward[i * 2] = p0 / z;
            backward[i * 2 + 1] = p1 / z;
        }
        begin = end;
    }

    Ok(DistanceHmmResult {
        log_likelihood: total_ll,
        viterbi,
        posterior: backward,
    })
}

fn deterministic_chr_shuffle(segments: &[HmmSegment], group_idx: usize) -> Vec<usize> {
    // Shuffle detection labels within each chromosome while leaving genomic
    // coordinates untouched. Real and null therefore see exactly the same gap
    // distribution; only the association between detection and position changes.
    let mut out = vec![0usize; segments.len()];
    let mut begin = 0usize;
    while begin < segments.len() {
        let chr = &segments[begin].chromosome;
        let mut end = begin + 1;
        while end < segments.len() && segments[end].chromosome == *chr {
            end += 1;
        }
        let mut idx = (begin..end).collect::<Vec<_>>();
        idx.sort_unstable_by_key(|&i| splitmix64((i as u64) ^ ((group_idx as u64 + 1) << 32)));
        for (dst, src) in (begin..end).zip(idx.into_iter()) {
            out[dst] = segments[src].detected[group_idx];
        }
        begin = end;
    }
    out
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

fn transition_count_chr(states: &[StateId], segments: &[HmmSegment]) -> usize {
    states
        .windows(2)
        .zip(segments.windows(2))
        .filter(|(state, seg)| seg[0].chromosome == seg[1].chromosome && state[0] != state[1])
        .count()
}

/// Inspect saved evidence; never recalculate module assignments or expression effects.
fn inspect_results(args: &InspectResultsCli) -> Result<()> {
    use std::io::{BufRead, BufReader};
    fn table(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines.next().context("empty TSV")??
            .split('\t').map(str::to_owned).collect::<Vec<_>>();
        let mut rows = Vec::new();
        for (i, line) in lines.enumerate() {
            let row = line?.split('\t').map(str::to_owned).collect::<Vec<_>>();
            anyhow::ensure!(row.len() == header.len(),
                "{}: line {} has {} columns, expected {}", path.display(), i + 2, row.len(), header.len());
            rows.push(row);
        }
        Ok((header, rows))
    }
    fn col(header: &[String], name: &str) -> Result<usize> {
        header.iter().position(|h| h == name)
            .with_context(|| format!("missing required column {name}"))
    }
    let (mh, modules) = table(&args.results.join("effect_modules.tsv"))?;
    let mi = col(&mh, "module")?;
    let gi = col(&mh, "genes")?;
    println!("Module effects (mean gene log2(case/control) by population):");
    println!("{}", mh.join("\t"));
    let mut found = false;
    for row in &modules {
        if args.module.is_some_and(|id| row[mi] != id.to_string()) { continue; }
        found = true;
        println!("{}", row.join("\t"));
    }
    anyhow::ensure!(found, "module {:?} not found", args.module);
    let (gh, genes) = table(&args.results.join("effect_clusters.tsv"))?;
    let gene_i = col(&gh, "gene")?;
    let module_i = col(&gh, "effect_cluster")?;
    let cluster_i = col(&gh, "cluster")?;
    let effect_i = col(&gh, "log2_ratio")?;
    for name in ["case_n", "control_n", "case_mean", "control_mean", "case_detected", "control_detected"] {
        col(&gh, name)?;
    }
    let mut grouped: BTreeMap<(usize, String), Vec<&Vec<String>>> = BTreeMap::new();
    for row in &genes {
        let id: usize = row[module_i].parse().with_context(|| format!("invalid module ID for {}", row[gene_i]))?;
        if args.module.is_some_and(|selected| selected != id) { continue; }
        grouped.entry((id, row[cluster_i].clone())).or_default().push(row);
    }
    for ((id, cluster), mut rows) in grouped {
        rows.sort_by(|a,b| {
            let ae = a[effect_i].parse::<f64>().unwrap_or(0.0).abs();
            let be = b[effect_i].parse::<f64>().unwrap_or(0.0).abs();
            be.total_cmp(&ae).then_with(|| a[gene_i].cmp(&b[gene_i]))
        });
        let count = rows.len();
        println!("\nModule {id}, population {cluster}: {count} genes (ranked by absolute log2 effect)");
        println!("{}", gh.join("\t"));
        for row in rows.into_iter().take(if args.top == 0 { count } else { args.top }) {
            println!("{}", row.join("\t"));
        }
    }
    println!("\nNote: these are existing gene-level means and effects; per-gene SD was not saved in this investigation.");
    println!("Module numbers are local to each investigation and are not matched across datasets.");
    let _ = gi;
    Ok(())
}
