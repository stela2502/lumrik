use statrs::distribution::{DiscreteCDF, Poisson};
use statrs::function::gamma::ln_gamma;

const MIN_PROB: f64 = 1e-300;

pub fn poisson_log_pmf(k: u32, lambda: f64) -> f64 {
    let lambda = lambda.max(1e-12);
    k as f64 * lambda.ln() - lambda - ln_gamma(k as f64 + 1.0)
}

pub fn negbin_log_pmf(k: u32, mean: f64, theta: f64) -> f64 {
    let mean = mean.max(1e-12);
    let theta = theta.max(1e-6);
    let k = k as f64;

    ln_gamma(k + theta)
        - ln_gamma(theta)
        - ln_gamma(k + 1.0)
        + theta * (theta / (theta + mean)).ln()
        + k * (mean / (theta + mean)).ln()
}

pub fn logsumexp(values: &[f64]) -> f64 {
    let m = values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    if !m.is_finite() {
        return m;
    }
    m + values.iter().map(|x| (*x - m).exp()).sum::<f64>().ln()
}

/// log P(X=x | X=A+T), where
/// A ~ Poisson(ambient_mean)
/// T ~ NegativeBinomial(true_mean, theta)
pub fn true_convolution_log_pmf(
    x: u32,
    ambient_mean: f64,
    true_mean: f64,
    theta: f64,
) -> f64 {
    let terms: Vec<f64> = (0..=x)
        .map(|ambient| {
            poisson_log_pmf(ambient, ambient_mean)
                + negbin_log_pmf(x - ambient, true_mean, theta)
        })
        .collect();
    logsumexp(&terms)
}

/// E[A | X=x, true-guide component].
pub fn expected_ambient_given_true(
    x: u32,
    ambient_mean: f64,
    true_mean: f64,
    theta: f64,
) -> f64 {
    let mut logs = Vec::with_capacity(x as usize + 1);
    for ambient in 0..=x {
        logs.push(
            poisson_log_pmf(ambient, ambient_mean)
                + negbin_log_pmf(x - ambient, true_mean, theta),
        );
    }

    let norm = logsumexp(&logs);
    if !norm.is_finite() {
        return 0.0;
    }

    logs.iter()
        .enumerate()
        .map(|(ambient, lp)| ambient as f64 * (*lp - norm).exp())
        .sum()
}

#[derive(Debug, Clone, Copy)]
pub struct PosteriorEvidence {
    pub probability: f64,
    pub log_odds: f64,
}

impl Default for PosteriorEvidence {
    fn default() -> Self {
        Self {
            probability: 0.0,
            log_odds: f64::NEG_INFINITY,
        }
    }
}

impl PosteriorEvidence {
    pub fn new(
        x: u32,
        ambient_mean: f64,
        prior_real: f64,
        true_mean: f64,
        theta: f64,
    ) -> Self {
        let prior =
            prior_real.clamp(1e-9, 1.0 - 1e-9);

        let l0 =
            (1.0 - prior).ln()
            + poisson_log_pmf(
                x,
                ambient_mean,
            );

        let l1 =
            prior.ln()
            + true_convolution_log_pmf(
                x,
                ambient_mean,
                true_mean,
                theta,
            );

        let log_odds =
            l1 - l0;

        let denom =
            logsumexp(&[l0, l1]);

        let probability =
            (l1 - denom)
                .exp()
                .clamp(0.0, 1.0);

        Self {
            probability,
            log_odds,
        }
    }
}

pub fn poisson_upper_tail(x: u32, lambda: f64) -> f64 {
    if x == 0 {
        return 1.0;
    }
    let Ok(p) = Poisson::new(lambda.max(1e-12)) else {
        return 1.0;
    };
    // P(X >= x) = 1 - P(X <= x-1)
    (1.0 - p.cdf((x - 1) as u64)).clamp(MIN_PROB, 1.0)
}

pub fn benjamini_hochberg(pvalues: &[f64]) -> Vec<f64> {
    if pvalues.is_empty() {
        return Vec::new();
    }

    let mut order: Vec<usize> = (0..pvalues.len()).collect();
    order.sort_by(|&a, &b| pvalues[a].total_cmp(&pvalues[b]));

    let m = pvalues.len() as f64;
    let mut q = vec![1.0; pvalues.len()];
    let mut running = 1.0_f64;

    for (rank0, &idx) in order.iter().enumerate().rev() {
        let rank = rank0 + 1;
        let candidate = (pvalues[idx] * m / rank as f64).min(1.0);
        running = running.min(candidate);
        q[idx] = running;
    }

    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bh_is_monotonic_in_rank() {
        let p = vec![0.001, 0.01, 0.2, 0.9];
        let q = benjamini_hochberg(&p);
        assert!(q[0] <= q[1]);
        assert!(q[1] <= q[2]);
        assert!(q[2] <= q[3]);
    }

    #[test]
    fn strong_signal_gets_large_posterior() {
        let z = PosteriorEvidence::new(50, 0.2, 0.05, 40.0, 10.0);
        assert!(z.probability > 0.99);
    }

    #[test]
    fn ambient_like_signal_stays_small() {
        let z = PosteriorEvidence::new(2, 2.0, 0.05, 40.0, 10.0);
        assert!(z.probability < 0.5);
    }
}
