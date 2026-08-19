
use clap::Args;

use crate::{
    background::BackgroundConfig,
    caller::CallConfig,
    model::FitConfig,
};

#[derive(Debug, Clone, Args)]
pub struct GuideModelCli {
    /// Dirichlet pseudocount used when estimating ambient
    /// guide frequencies from raw-only droplets.
    #[arg(long, default_value_t = 0.5)]
    pub ambient_alpha: f64,

    /// Minimum posterior probability required to call a guide genuine.
    #[arg(long, default_value_t = 0.95)]
    pub posterior_threshold: f64,

    /// Maximum Benjamini-Hochberg FDR for a genuine guide call.
    #[arg(long, default_value_t = 0.05)]
    pub fdr: f64,

    /// Maximum number of mixture-model iterations.
    #[arg(long, default_value_t = 500)]
    pub max_iterations: usize,

    /// Number of consecutive iterations with unchanged guide assignments
    /// required for biological convergence.
    #[arg(long, default_value_t = 3)]
    pub stable_iterations: usize,

    /// Minimum posterior probability required to call a guide genuine.
    #[arg(long, default_value_t = 0.95)]
    pub minimum_posterior: f64,

    /// Relative parameter-change tolerance used for mathematical convergence.
    #[arg(long, default_value_t = 1e-5)]
    pub convergence_tolerance: f64,

    /// Initial prior probability that an observed guide is genuine.
    #[arg(long, default_value_t = 0.05)]
    pub initial_prior_real: f64,

    /// Initial negative-binomial dispersion parameter.
    #[arg(long, default_value_t = 10.0)]
    pub initial_dispersion: f64,

    /// Minimum allowed mean for the genuine-guide component.
    #[arg(long, default_value_t = 0.5)]
    pub minimum_true_mean: f64,

    /// Lower bound for the cell-specific ambient burden.
    #[arg(long, default_value_t = 1e-6)]
    pub minimum_lambda: f64,

    /// Alpha parameter of the Beta prior on genuine-guide frequency.
    #[arg(long, default_value_t = 0.5)]
    pub prior_alpha: f64,

    /// Beta parameter of the Beta prior on genuine-guide frequency.
    #[arg(long, default_value_t = 9.5)]
    pub prior_beta: f64,
}


impl GuideModelCli {
    pub fn background_config(&self) -> BackgroundConfig {
        BackgroundConfig {
            alpha: self.ambient_alpha,
        }
    }

    pub fn fit_config(&self) -> FitConfig {
        FitConfig {
            max_iterations: self.max_iterations,
            tolerance: self.convergence_tolerance,

            initial_prior_real: self.initial_prior_real,
            initial_dispersion: self.initial_dispersion,

            minimum_true_mean: self.minimum_true_mean,
            minimum_lambda: self.minimum_lambda,

            prior_alpha: self.prior_alpha,
            prior_beta: self.prior_beta,

            stable_iterations_required: self.stable_iterations,
            minimum_posterior: self.minimum_posterior,

            // If this still exists in FitConfig but isn't part of
            // biological stopping anymore:
            posterior_tolerance: 1e-3,
        }
    }

    pub fn call_config(&self) -> CallConfig {
        CallConfig {
            minimum_posterior: self.posterior_threshold,
            maximum_fdr: self.fdr,
        }
    }
}