//! Generic hidden Markov models in log space.
//!
//! The crate deliberately knows nothing about proteins, genomes, peaks, BigWig,
//! ChIP, or any other biological data type. Biology enters through three small
//! traits: [`ObservationSource`], [`Emission`], and [`ResultProcessor`].

use std::f64::consts::PI;

const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StateId(pub usize);

#[derive(Debug, Clone, PartialEq)]
pub enum HmmError {
    NoStates,
    EmptyObservations,
    InitialLength { expected: usize, got: usize },
    TransitionLength { expected: usize, got: usize },
    InvalidProbability { what: &'static str, value: f64 },
    InvalidTransitionRow { state: usize, sum: f64 },
    InvalidInitialSum { sum: f64 },
    InvalidEmission(String),
}

impl std::fmt::Display for HmmError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for HmmError {}

/// Supplies an ordered observation sequence to the HMM.
///
/// `O` is returned by value on purpose: adapters can return a tiny packed or
/// computed observation without allocating. A protein adapter may return a
/// small amino-acid chemistry bitset; a BigWig adapter may return one signal
/// bin; a ChIP adapter may return a `(signal, background)` observation.
pub trait ObservationSource<O> {
    fn len(&self) -> usize;
    fn observation(&self, index: usize) -> O;
    fn is_empty(&self) -> bool { self.len() == 0 }
}

impl<O: Clone> ObservationSource<O> for [O] {
    fn len(&self) -> usize { <[O]>::len(self) }
    fn observation(&self, index: usize) -> O { self[index].clone() }
}

impl<O: Clone> ObservationSource<O> for Vec<O> {
    fn len(&self) -> usize { self.as_slice().len() }
    fn observation(&self, index: usize) -> O { self[index].clone() }
}

/// Probability of an observation given one hidden state.
pub trait Emission<O> {
    fn log_probability(&self, observation: &O) -> f64;
}

/// Optional training interface for emissions used by Baum-Welch.
pub trait TrainableEmission<O>: Emission<O> {
    type Accumulator;
    fn accumulator(&self) -> Self::Accumulator;
    fn accumulate(&self, accumulator: &mut Self::Accumulator, observation: &O, weight: f64);
    fn update(&mut self, accumulator: Self::Accumulator) -> Result<(), HmmError>;
}

/// Turns the mathematical result back into domain-specific output.
pub trait ResultProcessor<S> {
    type Output;
    fn process(&self, source: &S, result: &HmmResult) -> Self::Output;
}

#[derive(Debug, Clone)]
pub struct HmmResult {
    state_count: usize,
    observation_count: usize,
    log_likelihood: f64,
    viterbi: Vec<StateId>,
    posterior: Vec<f64>,
}

impl HmmResult {
    pub fn log_likelihood(&self) -> f64 { self.log_likelihood }
    pub fn viterbi(&self) -> &[StateId] { &self.viterbi }
    pub fn state_count(&self) -> usize { self.state_count }
    pub fn observation_count(&self) -> usize { self.observation_count }
    pub fn posterior(&self, position: usize, state: StateId) -> f64 {
        self.posterior[position * self.state_count + state.0]
    }
    pub fn posterior_row(&self, position: usize) -> &[f64] {
        let start = position * self.state_count;
        &self.posterior[start..start + self.state_count]
    }
}

/// Dense HMM with row-major transition matrix.
///
/// Public constructors accept ordinary probabilities. Internally all dynamic
/// programming is performed in natural-log space.
#[derive(Debug, Clone)]
pub struct Hmm<E> {
    initial: Vec<f64>,
    transition: Vec<f64>,
    emissions: Vec<E>,
}

impl<E> Hmm<E> {
    pub fn new(initial: Vec<f64>, transition: Vec<f64>, emissions: Vec<E>) -> Result<Self, HmmError> {
        let n = emissions.len();
        if n == 0 { return Err(HmmError::NoStates); }
        if initial.len() != n {
            return Err(HmmError::InitialLength { expected: n, got: initial.len() });
        }
        if transition.len() != n * n {
            return Err(HmmError::TransitionLength { expected: n * n, got: transition.len() });
        }
        validate_distribution(&initial, "initial")?;
        for state in 0..n {
            let row = &transition[state * n..(state + 1) * n];
            let sum: f64 = row.iter().sum();
            if !approx_one(sum) { return Err(HmmError::InvalidTransitionRow { state, sum }); }
            for &p in row { validate_probability(p, "transition")?; }
        }
        Ok(Self {
            initial: initial.into_iter().map(prob_to_log).collect(),
            transition: transition.into_iter().map(prob_to_log).collect(),
            emissions,
        })
    }

    pub fn state_count(&self) -> usize { self.emissions.len() }
    pub fn emissions(&self) -> &[E] { &self.emissions }
    pub fn emissions_mut(&mut self) -> &mut [E] { &mut self.emissions }

    pub fn initial_probabilities(&self) -> Vec<f64> {
        self.initial.iter().map(|x| x.exp()).collect()
    }

    pub fn transition_probabilities(&self) -> Vec<f64> {
        self.transition.iter().map(|x| x.exp()).collect()
    }

    pub fn infer<S, O>(&self, source: &S) -> Result<HmmResult, HmmError>
    where
        S: ObservationSource<O> + ?Sized,
        E: Emission<O>,
    {
        if source.is_empty() { return Err(HmmError::EmptyObservations); }
        let (forward, log_likelihood) = self.forward(source);
        let backward = self.backward(source);
        let posterior = self.posterior_from(&forward, &backward, log_likelihood, source.len());
        let viterbi = self.viterbi(source);
        Ok(HmmResult {
            state_count: self.state_count(),
            observation_count: source.len(),
            log_likelihood,
            viterbi,
            posterior,
        })
    }

    fn forward<S, O>(&self, source: &S) -> (Vec<f64>, f64)
    where S: ObservationSource<O> + ?Sized, E: Emission<O> {
        let n = self.state_count();
        let t_len = source.len();
        let mut f = vec![NEG_INF; t_len * n];
        let obs0 = source.observation(0);
        for s in 0..n {
            f[s] = self.initial[s] + self.emissions[s].log_probability(&obs0);
        }
        let mut scratch = vec![NEG_INF; n];
        for t in 1..t_len {
            let obs = source.observation(t);
            for dst in 0..n {
                for src in 0..n {
                    scratch[src] = f[(t - 1) * n + src] + self.transition[src * n + dst];
                }
                f[t * n + dst] = log_sum_exp(&scratch) + self.emissions[dst].log_probability(&obs);
            }
        }
        let ll = log_sum_exp(&f[(t_len - 1) * n..t_len * n]);
        (f, ll)
    }

    fn backward<S, O>(&self, source: &S) -> Vec<f64>
    where S: ObservationSource<O> + ?Sized, E: Emission<O> {
        let n = self.state_count();
        let t_len = source.len();
        let mut b = vec![0.0; t_len * n];
        let mut scratch = vec![NEG_INF; n];
        for t in (0..t_len - 1).rev() {
            let next_obs = source.observation(t + 1);
            for src in 0..n {
                for dst in 0..n {
                    scratch[dst] = self.transition[src * n + dst]
                        + self.emissions[dst].log_probability(&next_obs)
                        + b[(t + 1) * n + dst];
                }
                b[t * n + src] = log_sum_exp(&scratch);
            }
        }
        b
    }

    fn posterior_from(&self, f: &[f64], b: &[f64], ll: f64, t_len: usize) -> Vec<f64> {
        let n = self.state_count();
        let mut out = vec![0.0; t_len * n];
        for t in 0..t_len {
            let mut sum = 0.0;
            for s in 0..n {
                let p = (f[t * n + s] + b[t * n + s] - ll).exp();
                out[t * n + s] = p;
                sum += p;
            }
            if sum > 0.0 {
                for s in 0..n { out[t * n + s] /= sum; }
            }
        }
        out
    }

    fn viterbi<S, O>(&self, source: &S) -> Vec<StateId>
    where S: ObservationSource<O> + ?Sized, E: Emission<O> {
        let n = self.state_count();
        let t_len = source.len();
        let mut prev = vec![NEG_INF; n];
        let mut curr = vec![NEG_INF; n];
        let mut back = vec![0usize; t_len * n];
        let obs0 = source.observation(0);
        for s in 0..n { prev[s] = self.initial[s] + self.emissions[s].log_probability(&obs0); }
        for t in 1..t_len {
            let obs = source.observation(t);
            for dst in 0..n {
                let mut best_state = 0;
                let mut best = NEG_INF;
                for src in 0..n {
                    let score = prev[src] + self.transition[src * n + dst];
                    if score > best { best = score; best_state = src; }
                }
                curr[dst] = best + self.emissions[dst].log_probability(&obs);
                back[t * n + dst] = best_state;
            }
            std::mem::swap(&mut prev, &mut curr);
        }
        let mut state = argmax(&prev);
        let mut path = vec![StateId(0); t_len];
        path[t_len - 1] = StateId(state);
        for t in (1..t_len).rev() {
            state = back[t * n + state];
            path[t - 1] = StateId(state);
        }
        path
    }
}

impl<E> Hmm<E> {
    /// One or more Baum-Welch iterations over a single observation sequence.
    ///
    /// This updates initial probabilities, transitions, and any emission that
    /// implements [`TrainableEmission`]. `floor` prevents exact zero transition
    /// probabilities from permanently removing a path.
    pub fn baum_welch<S, O>(&mut self, source: &S, iterations: usize, floor: f64) -> Result<(), HmmError>
    where
        S: ObservationSource<O> + ?Sized,
        E: TrainableEmission<O>,
    {
        if source.is_empty() { return Err(HmmError::EmptyObservations); }
        if !(floor.is_finite() && floor >= 0.0) {
            return Err(HmmError::InvalidProbability { what: "floor", value: floor });
        }
        let n = self.state_count();
        let t_len = source.len();
        for _ in 0..iterations {
            let (f, ll) = self.forward(source);
            let b = self.backward(source);
            let gamma = self.posterior_from(&f, &b, ll, t_len);

            for s in 0..n { self.initial[s] = gamma[s].max(floor).ln(); }
            normalize_log_row(&mut self.initial);

            let mut trans_counts = vec![0.0; n * n];
            for t in 0..t_len.saturating_sub(1) {
                let next_obs = source.observation(t + 1);
                for src in 0..n {
                    for dst in 0..n {
                        let log_xi = f[t * n + src]
                            + self.transition[src * n + dst]
                            + self.emissions[dst].log_probability(&next_obs)
                            + b[(t + 1) * n + dst]
                            - ll;
                        trans_counts[src * n + dst] += log_xi.exp();
                    }
                }
            }
            for src in 0..n {
                let row = &mut trans_counts[src * n..(src + 1) * n];
                for x in row.iter_mut() { *x = (*x).max(floor); }
                let sum: f64 = row.iter().sum();
                if sum > 0.0 {
                    for dst in 0..n { self.transition[src * n + dst] = (row[dst] / sum).ln(); }
                }
            }

            let mut accumulators: Vec<E::Accumulator> = self.emissions.iter().map(E::accumulator).collect();
            for t in 0..t_len {
                let obs = source.observation(t);
                for s in 0..n {
                    self.emissions[s].accumulate(&mut accumulators[s], &obs, gamma[t * n + s]);
                }
            }
            for (emission, accumulator) in self.emissions.iter_mut().zip(accumulators) {
                emission.update(accumulator)?;
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct CategoricalEmission {
    log_probabilities: Vec<f64>,
}

impl CategoricalEmission {
    pub fn new(probabilities: Vec<f64>) -> Result<Self, HmmError> {
        validate_distribution(&probabilities, "categorical emission")?;
        Ok(Self { log_probabilities: probabilities.into_iter().map(prob_to_log).collect() })
    }
    pub fn probabilities(&self) -> Vec<f64> { self.log_probabilities.iter().map(|x| x.exp()).collect() }
}

impl Emission<usize> for CategoricalEmission {
    fn log_probability(&self, observation: &usize) -> f64 {
        self.log_probabilities.get(*observation).copied().unwrap_or(NEG_INF)
    }
}

#[derive(Debug, Clone)]
pub struct CategoricalAccumulator { counts: Vec<f64> }

impl TrainableEmission<usize> for CategoricalEmission {
    type Accumulator = CategoricalAccumulator;
    fn accumulator(&self) -> Self::Accumulator { CategoricalAccumulator { counts: vec![0.0; self.log_probabilities.len()] } }
    fn accumulate(&self, acc: &mut Self::Accumulator, observation: &usize, weight: f64) {
        if let Some(x) = acc.counts.get_mut(*observation) { *x += weight; }
    }
    fn update(&mut self, acc: Self::Accumulator) -> Result<(), HmmError> {
        let sum: f64 = acc.counts.iter().sum();
        if sum <= 0.0 { return Ok(()); }
        self.log_probabilities = acc.counts.into_iter().map(|x| if x == 0.0 { NEG_INF } else { (x / sum).ln() }).collect();
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GaussianEmission { mean: f64, variance: f64 }

impl GaussianEmission {
    pub fn new(mean: f64, variance: f64) -> Result<Self, HmmError> {
        if !mean.is_finite() || !variance.is_finite() || variance <= 0.0 {
            return Err(HmmError::InvalidEmission("Gaussian requires finite mean and variance > 0".into()));
        }
        Ok(Self { mean, variance })
    }
    pub fn mean(&self) -> f64 { self.mean }
    pub fn variance(&self) -> f64 { self.variance }
}

impl Emission<f64> for GaussianEmission {
    fn log_probability(&self, x: &f64) -> f64 {
        -0.5 * ((2.0 * PI * self.variance).ln() + (x - self.mean).powi(2) / self.variance)
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct GaussianAccumulator { weight: f64, sum: f64, sum_sq: f64 }

impl TrainableEmission<f64> for GaussianEmission {
    type Accumulator = GaussianAccumulator;
    fn accumulator(&self) -> Self::Accumulator { GaussianAccumulator::default() }
    fn accumulate(&self, acc: &mut Self::Accumulator, x: &f64, weight: f64) {
        acc.weight += weight; acc.sum += weight * x; acc.sum_sq += weight * x * x;
    }
    fn update(&mut self, acc: Self::Accumulator) -> Result<(), HmmError> {
        if acc.weight <= 0.0 { return Ok(()); }
        let mean = acc.sum / acc.weight;
        let variance = (acc.sum_sq / acc.weight - mean * mean).max(1e-12);
        self.mean = mean; self.variance = variance;
        Ok(())
    }
}

/// Discrete histogram emission for already-binned continuous data.
/// This intentionally mirrors the useful idea in Lumrik's 2008 Perl HMM:
/// arbitrary experimental values can be reduced to an empirical probability
/// function while the HMM remains oblivious to their biological meaning.
#[derive(Debug, Clone)]
pub struct HistogramEmission { inner: CategoricalEmission }

impl HistogramEmission {
    pub fn new(probabilities: Vec<f64>) -> Result<Self, HmmError> { Ok(Self { inner: CategoricalEmission::new(probabilities)? }) }
    pub fn probabilities(&self) -> Vec<f64> { self.inner.probabilities() }
}
impl Emission<usize> for HistogramEmission {
    fn log_probability(&self, observation: &usize) -> f64 { self.inner.log_probability(observation) }
}
impl TrainableEmission<usize> for HistogramEmission {
    type Accumulator = CategoricalAccumulator;
    fn accumulator(&self) -> Self::Accumulator { self.inner.accumulator() }
    fn accumulate(&self, acc: &mut Self::Accumulator, obs: &usize, weight: f64) { self.inner.accumulate(acc, obs, weight) }
    fn update(&mut self, acc: Self::Accumulator) -> Result<(), HmmError> { self.inner.update(acc) }
}

fn validate_distribution(values: &[f64], what: &'static str) -> Result<(), HmmError> {
    for &p in values { validate_probability(p, what)?; }
    let sum: f64 = values.iter().sum();
    if !approx_one(sum) { return Err(HmmError::InvalidInitialSum { sum }); }
    Ok(())
}
fn validate_probability(p: f64, what: &'static str) -> Result<(), HmmError> {
    if !p.is_finite() || !(0.0..=1.0).contains(&p) { return Err(HmmError::InvalidProbability { what, value: p }); }
    Ok(())
}
fn approx_one(x: f64) -> bool { (x - 1.0).abs() <= 1e-9 }
fn prob_to_log(p: f64) -> f64 { if p == 0.0 { NEG_INF } else { p.ln() } }
fn argmax(xs: &[f64]) -> usize {
    let mut best = 0; for i in 1..xs.len() { if xs[i] > xs[best] { best = i; } } best
}
fn log_sum_exp(xs: &[f64]) -> f64 {
    let max = xs.iter().copied().fold(NEG_INF, f64::max);
    if max == NEG_INF { return NEG_INF; }
    max + xs.iter().map(|x| (x - max).exp()).sum::<f64>().ln()
}
fn normalize_log_row(row: &mut [f64]) {
    let z = log_sum_exp(row); if z.is_finite() { for x in row { *x -= z; } }
}
