use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use sc_analysis::{AnalysisConfig, analyze_exon_matrix};

#[derive(Parser, Debug)]
#[command(
    name = "sc-analyse",
    version,
    about = "Run Lumrik's basic single-cell analysis pipeline on a Norn result or exonic Matrix Market folder"
)]
struct Cli {
    /// Norn sample/result folder, Nelrune output folder, or exonic Matrix Market folder.
    #[arg(value_name = "INPUT")]
    matrix: PathBuf,

    /// Output directory for normalized matrix, annotations, statistics, and plots.
    #[arg(short, long, default_value = "sc_analysis")]
    out: PathBuf,

    /// Retain only cells with at least this many total UMIs before analysis.
    /// By default sc-analyse keeps every cell present in the input matrix.
    #[arg(long)]
    min_umi_count: Option<usize>,

    /// Remove IG/TR V, D, and J segment genes before normalization and downstream GEX analysis.
    /// Raw UMI filtering still counts these molecules; constant-region genes remain.
    #[arg(long)]
    exclude_vdj_before_normalization: bool,

    /// Number of highest-variance surviving genes used for PCA.
    #[arg(long, default_value_t = 2000)]
    variable_genes: usize,

    /// Number of principal components to calculate.
    #[arg(long, default_value_t = 30)]
    pca_components: usize,

    /// Dimensions of the deterministic local-neighbour embedding used for clustering.
    #[arg(long, default_value_t = 10)]
    umap_components: usize,

    /// Number of neighbours used by the clustering embedding and 2-D UMAP visualisation.
    #[arg(long, default_value_t = 15)]
    umap_neighbors: usize,

    /// Number of optimization epochs used by the 2-D UMAP visualisation.
    #[arg(long, default_value_t = 200)]
    umap_epochs: usize,

    /// Number of initial k-means overclusters.
    #[arg(long, default_value_t = 30)]
    overclusters: usize,

    /// Maximum k-means iterations.
    #[arg(long, default_value_t = 50)]
    kmeans_iterations: usize,

    /// Merge clusters while their recomputed mean-expression Pearson correlation reaches this value.
    #[arg(long, default_value_t = 0.95)]
    merge_pearson: f64,

    /// Fraction of the deepest cells within each final cluster reserved as pseudo-deep holdout.
    #[arg(long, default_value_t = 0.10)]
    holdout_fraction: f64,

    /// Target library size for normalization of surviving genes (non-mito/ribo and, when requested, non-VDJ).
    #[arg(long, default_value_t = 10_000.0)]
    normalization_scale: f64,

    /// Number of deterministic pseudo-samples per final cluster for statistics.
    #[arg(long, default_value_t = 5)]
    pseudo_samples: usize,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let config = AnalysisConfig {
        min_umi_count: cli.min_umi_count,
        exclude_vdj_before_normalization: cli.exclude_vdj_before_normalization,
        variable_genes: cli.variable_genes,
        pca_components: cli.pca_components,
        umap_components: cli.umap_components,
        umap_neighbors: cli.umap_neighbors,
        umap_epochs: cli.umap_epochs,
        overclusters: cli.overclusters,
        kmeans_iterations: cli.kmeans_iterations,
        merge_pearson: cli.merge_pearson,
        holdout_fraction: cli.holdout_fraction,
        normalization_scale: cli.normalization_scale,
        pseudo_samples: cli.pseudo_samples,
        output_dir: cli.out.clone(),
    };

    eprintln!("sc-analyse");
    eprintln!("  input: {}", cli.matrix.display());
    eprintln!("  output: {}", cli.out.display());
    match config.min_umi_count {
        Some(n) => eprintln!("  minimum UMI count: {n}"),
        None => eprintln!("  minimum UMI count: disabled"),
    }
    eprintln!(
        "  V/D/J before normalization: {}",
        if config.exclude_vdj_before_normalization { "excluded" } else { "retained" }
    );
    eprintln!("  variable genes: {}", config.variable_genes);
    eprintln!("  PCA components: {}", config.pca_components);
    eprintln!(
        "  spatial patching: recursive PCA boxes, max(50 cells, 1% of retained cells)"
    );
    eprintln!(
        "  UMAP neighbours / epochs: {} / {}",
        config.umap_neighbors, config.umap_epochs
    );
    eprintln!("  merge Pearson threshold: {:.4}", config.merge_pearson);
    eprintln!(
        "  pseudo-deep holdout fraction: {:.3}",
        config.holdout_fraction
    );
    eprintln!("  normalization scale: {:.0}", config.normalization_scale);
    eprintln!("  pseudo-samples per cluster: {}", config.pseudo_samples);

    let (_, summary) = analyze_exon_matrix(&cli.matrix, &config)?;

    eprintln!("\nAnalysis complete");
    eprintln!(
        "  input: {} cells x {} features",
        summary.input_cells, summary.input_features
    );
    if summary.filtered_cells != summary.input_cells {
        eprintln!(
            "  UMI filter: {} -> {} cells",
            summary.input_cells, summary.filtered_cells
        );
    }
    eprintln!(
        "  mitochondrial features: {}",
        summary.mitochondrial_features
    );
    eprintln!("  ribosomal features: {}", summary.ribosomal_features);
    eprintln!(
        "  surviving expression features: {}",
        summary.surviving_features
    );
    eprintln!("  variable genes used: {}", summary.variable_genes);
    eprintln!(
        "  clusters: {} initial -> {} final",
        summary.initial_clusters, summary.merged_clusters
    );
    eprintln!("  final cluster sizes:");
    for (cluster, &cells) in summary.cluster_sizes.iter().enumerate() {
        eprintln!(
            "    cluster_{cluster:03}: {cells} ({:.1}%)",
            100.0 * cells as f64 / summary.filtered_cells.max(1) as f64
        );
    }
    eprintln!("  pseudo-deep holdout cells: {}", summary.holdout_cells);
    eprintln!("  stage timings:");
    for timing in &summary.timings {
        let seconds = timing.elapsed.as_secs_f64();
        if seconds >= 60.0 {
            eprintln!(
                "    {:<28} {}m {:.1}s",
                timing.stage,
                (seconds / 60.0).floor() as u64,
                seconds % 60.0
            );
        } else {
            eprintln!("    {:<28} {:.2}s", timing.stage, seconds);
        }
    }
    eprintln!("  results: {}", cli.out.display());
    Ok(())
}
