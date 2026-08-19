use std::collections::HashMap;

use anyhow::{bail, Result};
use rayon::prelude::*;

use crate::background::AmbientModel;
use crate::dataset::{GuideDataset, GuideObservation};
use crate::stats::{expected_ambient_given_true, PosteriorEvidence};

#[derive(Debug, Clone)]
pub struct FitConfig {
    pub max_iterations: usize,
    pub tolerance: f64,

    // Biological convergence diagnostics.
    pub posterior_tolerance: f64,
    pub stable_iterations_required: usize,
    pub minimum_posterior: f64,

    pub initial_prior_real: f64,
    pub initial_dispersion: f64,
    pub minimum_true_mean: f64,
    pub minimum_lambda: f64,
    pub prior_alpha: f64,
    pub prior_beta: f64,
}

#[derive(Debug, Clone)]
pub struct FitIteration {
    pub iteration: usize,
    pub max_parameter_delta: f64,
    pub max_posterior_delta: f64,
    pub mean_posterior_delta: f64,
    pub changed_calls: usize,
    pub total_observations: usize,
}

#[derive(Debug, Clone)]
pub struct GuideExpressionModel {
    /// P(a cell genuinely contains this guide).
    pub prior_real: f64,
    /// Mean true guide-expression count, excluding ambient molecules.
    pub mean: f64,
    /// NB size/dispersion parameter. Large values approach Poisson.
    pub theta: f64,
}

#[derive(Debug, Clone)]
pub struct ObservationFit {
    pub observation: GuideObservation,
    pub evidence: PosteriorEvidence,
    pub expected_ambient: f64,
    pub expected_true: f64,
}

#[derive(Debug, Clone)]
pub struct FittedModel {
    pub ambient: AmbientModel,
    pub guides: Vec<GuideExpressionModel>,
    pub lambda_by_cell: HashMap<u64, f64>,
    pub observations: Vec<ObservationFit>,
    pub iterations: usize,
    pub mathematical_converged: bool,
    pub biological_converged: bool,
    pub diagnostics: Vec<FitIteration>,
}

pub fn fit_mixture(
    data: &GuideDataset,
    ambient: AmbientModel,
    cfg: &FitConfig,
) -> Result<FittedModel> {
    if data.n_guides() != ambient.guide_probability.len() {
        bail!("Guide count differs between ambient model and filtered data");
    }

    if data.n_cells() == 0 {
        bail!("Cannot fit guide model without filtered cells");
    }

    let mut lambda_by_cell: HashMap<u64, f64> = data
        .cell_ids
        .iter()
        .copied()
        .map(|cell_id| {
            (
                cell_id,
                (data.cell_total(cell_id) as f64).max(cfg.minimum_lambda),
            )
        })
        .collect();

    let mut guides = initialize_guides(data, &ambient, &lambda_by_cell, cfg);

    // Keep guide boundaries once. Because observations are flattened in
    // guide-major order, the M-step can work on slices without allocating a
    // Vec<&ObservationFit> for every guide on every iteration.
    let mut guide_ranges = Vec::with_capacity(data.by_guide.len());
    let mut offset = 0usize;
    for guide in &data.by_guide {
        let end = offset + guide.len();
        guide_ranges.push(offset..end);
        offset = end;
    }

    // Flattening preserves guide-major order and gives Rayon one contiguous
    // mutable slice for the E-step.
    let mut observations: Vec<ObservationFit> = data
        .by_guide
        .iter()
        .flat_map(|guide| guide.iter().copied())
        .map(|observation| ObservationFit {
            observation,
            evidence: PosteriorEvidence::default(),
            expected_ambient: observation.count as f64,
            expected_true: 0.0,
        })
        .collect();

    let mut previous_posteriors: Vec<f64> = Vec::with_capacity(observations.len());
    let mut stable_iterations = 0usize;
    let mut mathematical_converged = false;
    let mut biological_converged = false;
    let mut iterations = 0usize;
    let mut diagnostics = Vec::with_capacity(cfg.max_iterations);

    for iter in 0..cfg.max_iterations {
        iterations = iter + 1;

        // E-step: each observation is independent for fixed model parameters.
        update_observations(
            &mut observations,
            &guides,
            &lambda_by_cell,
            &ambient,
        );

        let (changed_calls, max_posterior_delta, mean_posterior_delta) =
            posterior_change_stats(
                &previous_posteriors,
                &observations,
                cfg.minimum_posterior,
            );

        if !previous_posteriors.is_empty() {
            if changed_calls == 0 {
                stable_iterations += 1;
            } else {
                stable_iterations = 0;
            }

            if stable_iterations >= cfg.stable_iterations_required {
                biological_converged = true;
            }
        }

        // Reuse the allocated buffer rather than clear + extend through a
        // second temporary collection.
        previous_posteriors.clear();
        previous_posteriors.extend(
            observations
                .iter()
                .map(|state| state.evidence.probability),
        );

        let old_lambda = lambda_by_cell.clone();
        let old_guides = guides.clone();

        // M-step: cell-specific expected ambient burden.
        for value in lambda_by_cell.values_mut() {
            *value = cfg.minimum_lambda;
        }

        for state in &observations {
            *lambda_by_cell
                .entry(state.observation.cell_id)
                .or_insert(cfg.minimum_lambda) += state.expected_ambient;
        }

        // M-step: guide-specific true-expression component.
        let n_cells = data.n_cells() as f64;

        for guide_id in 0..guides.len() {
            let states = &observations[guide_ranges[guide_id].clone()];

            let sum_z: f64 = states
                .iter()
                .map(|state| state.evidence.probability)
                .sum();

            let sum_true: f64 = states
                .iter()
                .map(|state| state.expected_true)
                .sum();

            let prior_real = (sum_z + cfg.prior_alpha)
                / (n_cells + cfg.prior_alpha + cfg.prior_beta);

            let true_mean = if sum_z > 1e-8 {
                (sum_true / sum_z).max(cfg.minimum_true_mean)
            } else {
                cfg.minimum_true_mean
            };

            let mut weighted_variance = 0.0;

            if sum_z > 1e-8 {
                for state in states {
                    let posterior = state.evidence.probability;
                    let inferred_true = if posterior > 1e-12 {
                        state.expected_true / posterior
                    } else {
                        0.0
                    };

                    weighted_variance +=
                        posterior * (inferred_true - true_mean).powi(2);
                }

                weighted_variance /= sum_z;
            }

            let theta = if weighted_variance > true_mean + 1e-8 {
                (true_mean * true_mean / (weighted_variance - true_mean))
                    .clamp(0.05, 1e6)
            } else {
                // Very high theta means the NB approaches a Poisson.
                1e6
            };

            guides[guide_id] = GuideExpressionModel {
                prior_real: prior_real.clamp(1e-9, 1.0 - 1e-9),
                mean: true_mean,
                theta,
            };
        }

        let max_parameter_delta = parameter_delta(
            &old_lambda,
            &lambda_by_cell,
            &old_guides,
            &guides,
        );

        if max_parameter_delta < cfg.tolerance {
            mathematical_converged = true;
        }

        diagnostics.push(FitIteration {
            iteration: iter + 1,
            max_parameter_delta,
            max_posterior_delta,
            mean_posterior_delta,
            changed_calls,
            total_observations: observations.len(),
        });

        if mathematical_converged || biological_converged {
            break;
        }
    }

    // The loop ends after an M-step, therefore refresh evidence once using
    // the final parameters. This deliberately reuses the same E-step helper.
    update_observations(
        &mut observations,
        &guides,
        &lambda_by_cell,
        &ambient,
    );

    Ok(FittedModel {
        ambient,
        guides,
        lambda_by_cell,
        observations,
        iterations,
        mathematical_converged,
        biological_converged,
        diagnostics,
    })
}

/// Update posterior evidence and expected latent counts for all observations.
///
/// Keeping this in one place prevents the iterative E-step and final refresh
/// from drifting apart.
fn update_observations(
    observations: &mut [ObservationFit],
    guides: &[GuideExpressionModel],
    lambda_by_cell: &HashMap<u64, f64>,
    ambient: &AmbientModel,
) {
    observations.par_iter_mut().for_each(|state| {
        let obs = state.observation;
        let guide = &guides[obs.guide_id as usize];
        let lambda = lambda_by_cell[&obs.cell_id];
        let ambient_mean = lambda * ambient.p(obs.guide_id);

        let evidence = PosteriorEvidence::new(
            obs.count,
            ambient_mean,
            guide.prior_real,
            guide.mean,
            guide.theta,
        );

        let ambient_if_real = expected_ambient_given_true(
            obs.count,
            ambient_mean,
            guide.mean,
            guide.theta,
        );

        state.evidence = evidence;
        state.expected_ambient =
            (1.0 - evidence.probability) * obs.count as f64
                + evidence.probability * ambient_if_real;
        state.expected_true =
            evidence.probability * (obs.count as f64 - ambient_if_real).max(0.0);
    });
}

/// Compare only posterior probabilities between iterations.
/// Log-odds are retained as evidence for guide ranking, but are not part of
/// the historical convergence definition.
fn posterior_change_stats(
    previous: &[f64],
    observations: &[ObservationFit],
    minimum_posterior: f64,
) -> (usize, f64, f64) {
    if previous.is_empty() {
        return (0, f64::INFINITY, f64::INFINITY);
    }

    assert_eq!(
        previous.len(),
        observations.len(),
        "number/order of observations changed during fitting"
    );

    let (changed, max_delta, sum_delta) = previous
        .par_iter()
        .zip(observations.par_iter())
        .map(|(previous, state)| {
            let current = state.evidence.probability;
            let delta = (current - *previous).abs();
            let changed = usize::from(
                (*previous >= minimum_posterior)
                    != (current >= minimum_posterior),
            );
            (changed, delta, delta)
        })
        .reduce(
            || (0usize, 0.0_f64, 0.0_f64),
            |a, b| (a.0 + b.0, a.1.max(b.1), a.2 + b.2),
        );

    (
        changed,
        max_delta,
        sum_delta / observations.len() as f64,
    )
}

fn parameter_delta(
    old_lambda: &HashMap<u64, f64>,
    lambda_by_cell: &HashMap<u64, f64>,
    old_guides: &[GuideExpressionModel],
    guides: &[GuideExpressionModel],
) -> f64 {
    let lambda_delta = lambda_by_cell
        .iter()
        .map(|(&cell_id, &new_value)| {
            let old_value = old_lambda
                .get(&cell_id)
                .copied()
                .unwrap_or(new_value);
            relative_delta(old_value, new_value)
        })
        .fold(0.0_f64, f64::max);

    let guide_delta = old_guides
        .iter()
        .zip(guides)
        .map(|(old, new)| {
            relative_delta(old.prior_real, new.prior_real)
                .max(relative_delta(old.mean, new.mean))
                .max(relative_delta(old.theta, new.theta))
        })
        .fold(0.0_f64, f64::max);

    lambda_delta.max(guide_delta)
}

fn initialize_guides(
    data: &GuideDataset,
    ambient: &AmbientModel,
    lambda_by_cell: &HashMap<u64, f64>,
    cfg: &FitConfig,
) -> Vec<GuideExpressionModel> {
    data.by_guide
        .iter()
        .enumerate()
        .map(|(guide_id, obs)| {
            let mut excess: Vec<f64> = obs
                .iter()
                .map(|x| {
                    let ambient_mean =
                        lambda_by_cell[&x.cell_id] * ambient.p(guide_id as u32);
                    (x.count as f64 - ambient_mean).max(0.0)
                })
                .filter(|x| *x > 0.0)
                .collect();

            excess.sort_by(f64::total_cmp);

            let mean = if excess.is_empty() {
                cfg.minimum_true_mean
            } else {
                // Upper-half mean is deliberately resistant to the sea of
                // low-count ambient observations at initialization.
                let start = excess.len() / 2;
                (excess[start..].iter().sum::<f64>()
                    / (excess.len() - start) as f64)
                    .max(cfg.minimum_true_mean)
            };

            GuideExpressionModel {
                prior_real: cfg.initial_prior_real,
                mean,
                theta: cfg.initial_dispersion,
            }
        })
        .collect()
}

fn relative_delta(old: f64, new: f64) -> f64 {
    (old - new).abs() / old.abs().max(new.abs()).max(1e-9)
}
