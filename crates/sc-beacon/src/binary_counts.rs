use anyhow::{Result, bail};

use crate::stats::{PosteriorEvidence, expected_ambient_given_true};

#[derive(Debug, Clone)]
pub struct BinaryCountFitConfig {
    pub max_iterations: usize,
    pub tolerance: f64,
    pub initial_signal_prior: f64,
    pub initial_dispersion: f64,
    pub minimum_signal_mean: f64,
    pub minimum_background_mean: f64,
    pub prior_alpha: f64,
    pub prior_beta: f64,
}

impl Default for BinaryCountFitConfig {
    fn default() -> Self {
        Self {
            max_iterations: 500,
            tolerance: 1e-5,
            initial_signal_prior: 0.05,
            initial_dispersion: 10.0,
            minimum_signal_mean: 0.5,
            minimum_background_mean: 1e-6,
            prior_alpha: 0.5,
            prior_beta: 9.5,
        }
    }
}

#[derive(Debug, Clone)]
pub struct BinaryCountFit {
    /// Mean of the Poisson background component.
    pub background_mean: f64,
    /// Prior probability that an observation belongs to the signal component.
    pub signal_prior: f64,
    /// Mean NB signal count excluding the Poisson background contribution.
    pub signal_mean: f64,
    /// NB size/dispersion parameter. Large values approach Poisson.
    pub theta: f64,
    /// Posterior P(signal | count), in the same order as the input counts.
    pub posteriors: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
}

/// Fit a two-component count mixture using Beacon's existing evidence model.
///
/// Background observations are modeled as Poisson(background_mean). Signal
/// observations are modeled as Poisson(background_mean) + NB(signal_mean,
/// theta). The returned posterior for each input value is P(signal | count).
pub fn fit_binary_counts(counts: &[u32]) -> Result<BinaryCountFit> {
    fit_binary_counts_with_config(counts, &BinaryCountFitConfig::default())
}

pub fn fit_binary_counts_with_config(
    counts: &[u32],
    cfg: &BinaryCountFitConfig,
) -> Result<BinaryCountFit> {
    if counts.is_empty() {
        bail!("Cannot fit binary count mixture without observations");
    }

    if cfg.max_iterations == 0 {
        bail!("Binary count mixture requires at least one iteration");
    }

    let n = counts.len() as f64;
    let (mut background_mean, mut signal_mean) = initialize(counts, cfg);
    let mut signal_prior = cfg.initial_signal_prior.clamp(1e-9, 1.0 - 1e-9);
    let mut theta = cfg.initial_dispersion.max(0.05);
    let mut posteriors = vec![0.0; counts.len()];
    let mut iterations = 0usize;
    let mut converged = false;

    for iteration in 0..cfg.max_iterations {
        iterations = iteration + 1;

        let old_background_mean = background_mean;
        let old_signal_prior = signal_prior;
        let old_signal_mean = signal_mean;
        let old_theta = theta;

        let mut sum_signal = 0.0;
        let mut sum_expected_ambient = 0.0;
        let mut sum_expected_signal = 0.0;
        let mut inferred_signal = Vec::with_capacity(counts.len());

        for (&count, posterior) in counts.iter().zip(&mut posteriors) {
            let evidence =
                PosteriorEvidence::new(count, background_mean, signal_prior, signal_mean, theta);

            let ambient_if_signal =
                expected_ambient_given_true(count, background_mean, signal_mean, theta);

            let p = evidence.probability;
            let expected_ambient = (1.0 - p) * count as f64 + p * ambient_if_signal;
            let expected_signal = p * (count as f64 - ambient_if_signal).max(0.0);

            *posterior = p;
            sum_signal += p;
            sum_expected_ambient += expected_ambient;
            sum_expected_signal += expected_signal;
            inferred_signal.push(if p > 1e-12 { expected_signal / p } else { 0.0 });
        }

        background_mean = (sum_expected_ambient / n).max(cfg.minimum_background_mean);

        signal_prior = ((sum_signal + cfg.prior_alpha) / (n + cfg.prior_alpha + cfg.prior_beta))
            .clamp(1e-9, 1.0 - 1e-9);

        signal_mean = if sum_signal > 1e-8 {
            (sum_expected_signal / sum_signal).max(cfg.minimum_signal_mean)
        } else {
            cfg.minimum_signal_mean
        };

        if sum_signal > 1e-8 {
            let weighted_variance = posteriors
                .iter()
                .zip(&inferred_signal)
                .map(|(&p, &value)| p * (value - signal_mean).powi(2))
                .sum::<f64>()
                / sum_signal;

            theta = if weighted_variance > signal_mean + 1e-8 {
                (signal_mean * signal_mean / (weighted_variance - signal_mean)).clamp(0.05, 1e6)
            } else {
                1e6
            };
        }

        let delta = relative_delta(old_background_mean, background_mean)
            .max(relative_delta(old_signal_prior, signal_prior))
            .max(relative_delta(old_signal_mean, signal_mean))
            .max(relative_delta(old_theta, theta));

        if delta < cfg.tolerance {
            converged = true;
            break;
        }
    }

    // Refresh posteriors once with the final parameters, since the loop ends
    // after the M-step.
    for (&count, posterior) in counts.iter().zip(&mut posteriors) {
        *posterior =
            PosteriorEvidence::new(count, background_mean, signal_prior, signal_mean, theta)
                .probability;
    }

    Ok(BinaryCountFit {
        background_mean,
        signal_prior,
        signal_mean,
        theta,
        posteriors,
        iterations,
        converged,
    })
}

fn initialize(counts: &[u32], cfg: &BinaryCountFitConfig) -> (f64, f64) {
    let mut sorted = counts.to_vec();
    sorted.sort_unstable();

    let split = (sorted.len() / 2).max(1);
    let lower = &sorted[..split];
    let upper = &sorted[split.min(sorted.len())..];

    let background_mean = (lower.iter().map(|&x| x as f64).sum::<f64>() / lower.len() as f64)
        .max(cfg.minimum_background_mean);

    let upper_mean = if upper.is_empty() {
        sorted[sorted.len() - 1] as f64
    } else {
        upper.iter().map(|&x| x as f64).sum::<f64>() / upper.len() as f64
    };

    let signal_mean = (upper_mean - background_mean).max(cfg.minimum_signal_mean);

    (background_mean, signal_mean)
}

fn relative_delta(old: f64, new: f64) -> f64 {
    (old - new).abs() / old.abs().max(new.abs()).max(1e-9)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn separates_low_background_from_high_signal() {
        let counts = [
            0, 0, 1, 1, 1, 2, 2, 1, 0, 2, 35, 40, 45, 50, 55, 60, 42, 48, 52, 58,
        ];

        let fit = fit_binary_counts(&counts).unwrap();

        assert!(fit.background_mean < 5.0);
        assert!(fit.signal_mean > 20.0);

        for &p in &fit.posteriors[..10] {
            assert!(p < 0.5, "background posterior was {p}");
        }

        for &p in &fit.posteriors[10..] {
            assert!(p > 0.95, "signal posterior was {p}");
        }
    }

    #[test]
    fn preserves_input_order() {
        let counts = [50, 1, 45, 0, 60, 2];
        let fit = fit_binary_counts(&counts).unwrap();

        assert!(fit.posteriors[0] > fit.posteriors[1]);
        assert!(fit.posteriors[2] > fit.posteriors[3]);
        assert!(fit.posteriors[4] > fit.posteriors[5]);
    }

    #[test]
    fn rejects_empty_input() {
        assert!(fit_binary_counts(&[]).is_err());
    }
}
