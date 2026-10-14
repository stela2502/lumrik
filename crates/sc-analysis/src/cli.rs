use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use crate::{AnalysisConfig, analyze_exon_matrix};

#[derive(Parser, Debug)]
#[command(
    name = "sc-analyse",
    version,
    about = "Run Lumrik's basic single-cell analysis pipeline on a Norn result, expression Matrix Market folder, or H5AD file"
)]
pub struct AnalysisCli {
    /// Norn sample/result folder, Nelrune output folder, expression Matrix Market folder, or H5AD file.
    #[arg(value_name = "INPUT")]
    matrix: PathBuf,

    /// Output directory for normalized matrix, annotations, statistics, and plots.
    #[arg(short, long, default_value = "sc_analysis")]
    out: PathBuf,

    /// Minimum surviving UMIs per cell, and exact fixed-depth resampling target.
    #[arg(long, default_value_t = 800)]
    min_umi_count: usize,

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

    /// Number of neighbours used by the clustering embedding and 2-D UMAP visualisation.
    #[arg(long, default_value_t = 15)]
    umap_neighbors: usize,

    /// Number of optimization epochs used by the 2-D UMAP visualisation.
    #[arg(long, default_value_t = 200)]
    umap_epochs: usize,

    /// Minimum target size for spatial PCA patches (cells).
    #[arg(long, default_value_t = 50)]
    spatial_min_cells: usize,

    /// Fraction of retained cells setting the spatial patch size floor (0 disables).
    #[arg(long, default_value_t = 0.01)]
    spatial_fraction: f64,

    /// Merge clusters while their recomputed mean-expression Pearson correlation reaches this value.
    #[arg(long, default_value_t = 0.95)]
    merge_pearson: f64,

    /// Maximum fraction of cells in a merged cluster (0 disables the limit).
    #[arg(long, default_value_t = 0.10)]
    max_cluster_fraction: f64,

    /// Merge PCA spatial patches by full-expression Pearson to exactly this many clusters.
    /// Overrides the Pearson merge stopping rule.
    #[arg(long)]
    target_clusters: Option<usize>,

    /// Fraction of the deepest cells within each final cluster reserved as pseudo-deep holdout.
    #[arg(long, default_value_t = 0.10)]
    holdout_fraction: f64,

    /// Number of deterministic pseudo-samples per final cluster for statistics.
    #[arg(long, default_value_t = 5)]
    pseudo_samples: usize,

    /// Minimum fraction of all cells in a reference cluster used as a mixture parent.
    #[arg(long, default_value_t = 0.01)]
    multiplet_min_parent_fraction: f64,

    /// Minimum relative improvement over best single-cluster fit for a multiplet explanation.
    #[arg(long, default_value_t = 0.35)]
    multiplet_min_improvement: f64,
}

pub fn run_analysis(cli: AnalysisCli) -> Result<()> {
    anyhow::ensure!(cli.spatial_min_cells > 0, "--spatial-min-cells must be positive");
    anyhow::ensure!(cli.spatial_fraction.is_finite() && (0.0..=1.0).contains(&cli.spatial_fraction), "--spatial-fraction must be between 0 and 1");
    anyhow::ensure!(cli.max_cluster_fraction.is_finite() && (0.0..=1.0).contains(&cli.max_cluster_fraction), "--max-cluster-fraction must be between 0 and 1");
    anyhow::ensure!(cli.target_clusters.is_none_or(|n| n > 0), "--target-clusters must be positive");
    anyhow::ensure!(cli.min_umi_count > 0, "--min-umi-count must be positive");
    anyhow::ensure!(cli.multiplet_min_parent_fraction.is_finite() && (0.0..=1.0).contains(&cli.multiplet_min_parent_fraction), "--multiplet-min-parent-fraction must be between 0 and 1");
    anyhow::ensure!((0.0..1.0).contains(&cli.multiplet_min_improvement), "--multiplet-min-improvement must be in [0,1)");
    let config = AnalysisConfig {
        min_umi_count: Some(cli.min_umi_count),
        exclude_vdj_before_normalization: cli.exclude_vdj_before_normalization,
        variable_genes: cli.variable_genes,
        pca_components: cli.pca_components,
        umap_neighbors: cli.umap_neighbors,
        umap_epochs: cli.umap_epochs,
        spatial_min_cells: cli.spatial_min_cells,
        spatial_fraction: cli.spatial_fraction,
        merge_pearson: cli.merge_pearson,
        max_cluster_fraction: cli.max_cluster_fraction,
        target_clusters: cli.target_clusters,
        holdout_fraction: cli.holdout_fraction,
        pseudo_samples: cli.pseudo_samples,
        multiplet_min_parent_fraction: cli.multiplet_min_parent_fraction,
        multiplet_min_improvement: cli.multiplet_min_improvement,
        output_dir: cli.out.clone(),
    };

    eprintln!("sc-analyse");
    eprintln!("  input: {}", cli.matrix.display());
    eprintln!("  output: {}", cli.out.display());
    eprintln!("  minimum surviving UMI count / resampling target: {}", cli.min_umi_count);
    eprintln!(
        "  V/D/J before normalization: {}",
        if config.exclude_vdj_before_normalization {
            "excluded"
        } else {
            "retained"
        }
    );
    eprintln!("  variable genes: {}", config.variable_genes);
    eprintln!("  PCA components: {}", config.pca_components);
    eprintln!("  spatial patching: recursive PCA boxes, max({} cells, {:.2}% of retained cells)", config.spatial_min_cells, config.spatial_fraction * 100.0);
    eprintln!(
        "  UMAP neighbours / epochs: {} / {}",
        config.umap_neighbors, config.umap_epochs
    );
    if let Some(k) = config.target_clusters {
        eprintln!("  merge mode: full-expression Pearson, target {k} clusters (rho cutoff ignored)");
    } else {
        eprintln!("  merge Pearson threshold: {:.4}", config.merge_pearson);
    }
    eprintln!(
        "  pseudo-deep holdout fraction: {:.3}",
        config.holdout_fraction
    );
    eprintln!("  user-facing expression: log1p scale 10000; PCA/UMAP/Pearson/multiplets: resampled counts, seed 42");
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
