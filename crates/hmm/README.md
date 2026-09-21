# `hmm`

A small, generic, log-space Hidden Markov Model engine for Lumrik.

This crate is intentionally **not a protein HMM**, **not a ChIP peak caller**, and **not a BigWig model**. It is the reusable mathematical engine underneath all of those things.

The design grew out of an older HMM implementation written for experimental genomics in 2008. The useful idea from that implementation survives here: the HMM should ask only two biological questions — *how likely is a state transition?* and *how likely is this observation in this state?* Everything else belongs outside the mathematical core.

The 2026 Rust version makes that separation explicit and fast:

```text
 biological data
      │
      ▼
ObservationSource<O>
      │
      ▼
 observations O
      │
      ▼
┌────────────────────────────┐
│            Hmm<E>          │
│                            │
│  forward / backward        │
│  posterior probabilities   │
│  Viterbi path              │
│  Baum-Welch training       │
└──────────────┬─────────────┘
               │
               ▼
           HmmResult
               │
               ▼
      ResultProcessor<S>
               │
       biological answer
```

## Why this crate exists

Lumrik increasingly needs to recognize *patterns along ordered biological data* rather than only exact strings or fixed thresholds.

Examples include:

- a protein represented by amino-acid identity plus per-residue biochemical properties;
- a learned protein/domain pattern scanned across another proteome;
- BigWig signal over a genomic interval;
- ChIP/CUT&RUN counts together with an input/background track;
- multiple chromatin tracks such as ATAC, H3K4me3, H3K27ac and H3K27me3;
- methylation or copy-number segmentation;
- motif/state sequences;
- any future ordered observation stream for which hidden states explain local structure.

Those applications should **not** each implement their own forward/backward algorithm. They should implement small adapters that tell this crate what an observation is and how a state emits it.

## Design goals

1. **Biology-independent core.** No dependency on `IntToProt`, BigWig, BED, BAM, or genomic coordinate types.
2. **Generic observations.** An observation may be an integer symbol, a floating-point signal, a packed chemistry bitset, a count pair, or a multi-track struct.
3. **Log-space numerics.** Dynamic programming stays in log space to avoid probability underflow.
4. **Dense hot loops.** States are integer positions and transitions are a row-major `Vec<f64>`. No state-name hash lookups occur in inference loops.
5. **Training is optional.** Inference requires only `Emission<O>`. Baum-Welch additionally requires `TrainableEmission<O>`.
6. **Interpretation is outside the core.** The HMM returns mathematical results; a `ResultProcessor` turns them into domains, peaks, genomic intervals, candidate states, etc.
7. **No forced allocations in adapters.** `ObservationSource` returns an observation value. Small computed/packed observations can therefore be generated on demand.

## Core traits

### `ObservationSource<O>`

```rust
pub trait ObservationSource<O> {
    fn len(&self) -> usize;
    fn observation(&self, index: usize) -> O;
}
```

This is the input boundary. The HMM does not care where observations came from.

A protein adapter might return:

```rust
#[derive(Clone, Copy)]
struct AminoAcidChemistry {
    residue: u8,
    chemistry: u16,
}
```

A ChIP adapter might return:

```rust
#[derive(Clone, Copy)]
struct ChipObservation {
    chip: u32,
    input: u32,
}
```

A chromatin adapter could return several tracks at once:

```rust
#[derive(Clone, Copy)]
struct ChromatinObservation {
    atac: f32,
    h3k4me3: f32,
    h3k27ac: f32,
    h3k27me3: f32,
}
```

None of these types belong in `hmm`; they belong in the biological adapter crate that understands them.

Slices and `Vec<O>` already implement `ObservationSource<O>` for cloneable observations, which keeps tests and simple callers trivial.

### `Emission<O>`

```rust
pub trait Emission<O> {
    fn log_probability(&self, observation: &O) -> f64;
}
```

This is the key abstraction inherited from the old HMM's `prob_for_observed()` idea.

An emission answers exactly one question:

> Given that the hidden state is this state, what is the log probability of observing `O`?

The HMM never needs to know why.

### `TrainableEmission<O>`

```rust
pub trait TrainableEmission<O>: Emission<O> {
    type Accumulator;

    fn accumulator(&self) -> Self::Accumulator;
    fn accumulate(
        &self,
        accumulator: &mut Self::Accumulator,
        observation: &O,
        posterior_weight: f64,
    );
    fn update(&mut self, accumulator: Self::Accumulator) -> Result<(), HmmError>;
}
```

This is the boundary used by Baum-Welch training. The HMM computes posterior state weights; the emission decides how those weighted observations update its own parameters.

That separation lets a Gaussian learn a mean/variance while a categorical model learns symbol frequencies and a future negative-binomial model learns count parameters — without changing the HMM engine.

### `ResultProcessor<S>`

```rust
pub trait ResultProcessor<S> {
    type Output;

    fn process(&self, source: &S, result: &HmmResult) -> Self::Output;
}
```

This is deliberately a small output-side adapter. The mathematical result should not know what a protein domain or genomic peak is.

Possible processors include:

```text
ProteinResultProcessor   -> protein/domain ranges
BigWigResultProcessor    -> genomic intervals
ChipResultProcessor      -> enriched/bound regions
ChromatinResultProcessor -> chromatin-state segmentation
```

## Creating a model

`Hmm::new` accepts ordinary probabilities because that is the least surprising public API. They are validated and immediately converted to natural-log probabilities internally.

For two states:

```rust
use hmm::{GaussianEmission, Hmm};

let emissions = vec![
    GaussianEmission::new(0.0, 0.25)?,
    GaussianEmission::new(5.0, 0.25)?,
];

let model = Hmm::new(
    vec![0.99, 0.01],
    vec![
        0.98, 0.02, // state 0 -> 0, 1
        0.02, 0.98, // state 1 -> 0, 1
    ],
    emissions,
)?;
```

Transitions are row-major:

```text
transition[source_state * N + destination_state]
```

Every transition row must sum to one, as must the initial-state probabilities.

## Inference

```rust
let signal = vec![0.0, 0.1, -0.2, 0.2, 4.8, 5.1, 5.0, 5.2];
let result = model.infer(&signal)?;

println!("log P(observations) = {}", result.log_likelihood());
println!("best path = {:?}", result.viterbi());

for position in 0..result.observation_count() {
    println!("{}: {:?}", position, result.posterior_row(position));
}
```

`HmmResult` contains:

- total sequence log-likelihood;
- Viterbi most-likely state path;
- posterior probability of every state at every position.

Posterior probabilities are often more biologically useful than the Viterbi path. A boundary can legitimately be uncertain, and callers should be able to retain that uncertainty rather than forcing every position into one hard state.

## Training

Trainable emissions can be optimized together with initial and transition probabilities using Baum-Welch:

```rust
let mut model = Hmm::new(initial, transitions, emissions)?;
model.baum_welch(&observations, 10, 1e-9)?;
```

The final argument is a probability floor used for transitions. It prevents an exact zero introduced during training from permanently deleting a possible transition.

The current trainer operates on one observation sequence. Multiple-sequence training should accumulate sufficient statistics across sequences before each M-step rather than repeatedly fitting one sequence after another; that is a planned extension.

## Built-in emissions

### `CategoricalEmission`

For discrete symbols represented as `usize`:

```rust
use hmm::CategoricalEmission;

let mostly_a = CategoricalEmission::new(vec![0.90, 0.05, 0.05])?;
```

This is useful for simple alphabets and also provides the machinery underneath `HistogramEmission`.

### `GaussianEmission`

For scalar continuous observations:

```rust
use hmm::GaussianEmission;

let background = GaussianEmission::new(0.0, 1.0)?;
let enriched   = GaussianEmission::new(4.0, 2.0)?;
```

It is trainable by posterior-weighted mean and variance.

### `HistogramEmission`

For continuous measurements that have already been converted into bins. This deliberately preserves one useful feature of the old Perl implementation: empirical distributions can model awkward experimental signals without requiring a closed-form distribution.

## Protein chemistry example

A future protein adapter can expose both residue identity and biochemical properties without unpacking a whole protein into strings:

```rust
#[derive(Clone, Copy)]
struct AaObservation {
    aa: u8,
    chemistry: u16,
}

struct ProteinChemistryEmission {
    // learned state parameters
}

impl hmm::Emission<AaObservation> for ProteinChemistryEmission {
    fn log_probability(&self, aa: &AaObservation) -> f64 {
        // score exact identity, charge, hydrophobicity, aromaticity,
        // size, polarity, etc. according to this hidden state's model
        todo!()
    }
}
```

Then the same engine can learn a state pattern from known proteins/domains and scan another protein or proteome.

The important point is that `hmm` itself never depends on `IntToProt`.

```text
IntToProt / protein adapter
          │
          ▼
    AaObservation
          │
          ▼
      Hmm<ProteinChemistryEmission>
```

## BigWig / genomic signal example

A BigWig adapter can expose one value or one multi-track observation per genomic bin:

```rust
#[derive(Clone, Copy)]
struct SignalObservation {
    signal: f32,
}
```

or:

```rust
#[derive(Clone, Copy)]
struct MultiTrackObservation {
    atac: f32,
    h3k4me3: f32,
    h3k27ac: f32,
}
```

The only new HMM-side work is an appropriate emission model. Coordinate handling, missing bins, chromosome boundaries and conversion of state runs back to genomic intervals belong in the adapter/result processor.

## ChIP plus background example

Do **not** force ChIP/input into a precomputed ratio merely because the HMM expects a scalar. It does not.

```rust
#[derive(Clone, Copy)]
struct ChipObservation {
    chip: u32,
    input: u32,
}
```

A future count emission can model the joint/count relationship directly, for example with Poisson or negative-binomial components. This keeps the raw evidence available to the statistical model.

## Multi-track / TF-state models

Nothing restricts an observation to one measurement. A model can combine several measurements at every genomic position/bin:

```rust
struct LocusObservation {
    atac: f32,
    tf_a_chip: f32,
    tf_b_chip: f32,
    motif_a: f32,
    motif_b: f32,
}
```

The posterior can then answer questions of the form:

```text
position/bin 1234
  background      0.004
  TF-A-like       0.021
  TF-B-like       0.913
  TF-C-like       0.038
  unexplained     0.024
```

The names are biological interpretation supplied by the caller. Internally they are just `StateId(0)`, `StateId(1)`, ... so state labels never enter the numerical hot loop.

## Unknown / X states

An `X`, `unknown`, or `background` state is not magic open-set recognition. It should be represented as an ordinary state with a deliberately broad emission distribution. Sequence-level log-likelihood and diffuse posterior probabilities are also useful signals that none of the modeled states explain an observation well.

This distinction matters: an HMM should not confidently invent a known biological label merely because all alternatives are bad.

## Numerical implementation

The forward and backward algorithms operate in log space. Sums of probabilities use log-sum-exp:

```text
log(exp(a) + exp(b) + ...)
```

computed after subtracting the largest log probability. This avoids the catastrophic underflow that occurs when multiplying many small probabilities in ordinary floating-point space.

The Viterbi algorithm also operates in log space but replaces summation with maximization.

Transitions are stored densely and row-major. This is intentional. Biological models commonly have a modest number of hidden states, and contiguous arrays avoid hash/string lookups in the inner dynamic-programming loops.

## Complexity

For `T` observations and `N` states:

- forward: `O(T * N^2)`;
- backward: `O(T * N^2)`;
- Viterbi: `O(T * N^2)`;
- full posterior storage: `O(T * N)`;
- dense transition matrix: `O(N^2)`.

`infer()` currently stores forward, backward and posterior matrices. That is appropriate for model development and moderate regions, but not the final answer for whole-genome scans. A later scanning API should use bounded-memory/chunked inference where the consumer only needs hits or state runs.

## What this crate deliberately does not do yet

### Profile HMM silent states

Classic protein profile HMMs use match/insert/delete topology, including silent delete states. A conventional HMM over one emitted observation per step does not represent silent states directly.

The current crate is therefore **not yet a drop-in HMMER replacement**. Protein chemistry/state models already fit the engine, but a full profile-HMM layer should add explicit topology/silent-state support rather than smuggling delete states into ordinary emissions.

### Hidden semi-Markov models

An ordinary HMM implies geometrically distributed state durations. Some genomic states have stronger duration structure. If that matters, an HSMM/duration extension belongs above or beside this core rather than complicating the initial implementation.

### Negative-binomial / Poisson / multivariate emissions

These are natural next emissions for sequencing counts and multi-track data. They should be added as independent emission implementations, not special cases in the HMM algorithms.

### Multi-sequence Baum-Welch

Training on many proteins, chromosomes, peaks or genomic intervals should aggregate expected sufficient statistics across independent sequences before updating the model. The current single-sequence implementation establishes the trait boundary and mathematics first.

### Serialization

Model persistence should be introduced after the model schema settles. When added, it should have an explicit format/version rather than serializing implementation details blindly.

## Biological initialization / weak supervision

The old HMM had a useful `HMM_hypothesis` idea: initialize state distributions from weak biological knowledge such as high/low quantiles or known regions, then let the HMM refine the model.

That idea should return, but **not inside this crate's mathematical core**.

Examples:

```text
ChIP:
  background <- low-signal / control-like regions
  bound      <- known positive or high-enrichment regions

Protein:
  TM-like    <- known transmembrane segments
  signal     <- known signal peptides
  other      <- remaining residues
```

A future initializer can construct `Hmm<E>` values from those hints. The core HMM remains ignorant of what the hints mean.

## Performance direction

The first priority is correctness and a stable abstraction. The layout is nevertheless chosen with performance in mind:

- integer state IDs;
- contiguous transition matrix;
- contiguous dynamic-programming matrices;
- no state-name maps in hot loops;
- observations supplied on demand;
- monomorphized Rust traits rather than a required dynamic-dispatch layer.

Likely later optimizations include:

- reusing scratch buffers across inference calls;
- sparse/topology-aware transitions where appropriate;
- chunked/streaming genome scans;
- parallel inference across independent proteins, chromosomes or regions;
- specialized profile-HMM kernels;
- SIMD only after profiling demonstrates a useful target.

Do not optimize the biological adapter API away: the central goal is to keep one correct inference engine reusable across very different data.

## Testing philosophy

The crate starts with synthetic tests where the correct answer is obvious:

1. a categorical two-regime sequence must switch states at the expected boundary;
2. a continuous Gaussian low/high signal must use the same inference engine and recover the same conceptual segmentation;
3. Baum-Welch must move poorly initialized Gaussian emissions toward the generating regimes.

Those tests prove the abstraction before real proteins, BigWigs or ChIP data are introduced.

Real integration tests should then live with the adapters and assert biological facts against tiny retained fixtures.

## Relationship to the 2008 implementation

The old Perl implementation already contained several sound ideas:

- log-space forward calculations;
- backward probabilities;
- posterior state probabilities;
- transition re-estimation;
- empirical/histogram emissions;
- weak biological hypotheses for initializing state distributions.

Its main limitation was architectural rather than conceptual: states, transitions, emissions, training state, plotting, serialization and experiment-specific code were tightly coupled, and hash/method lookups dominated the inner loops.

This crate keeps the useful statistical model and gives each responsibility a narrow Rust interface.

In short:

```text
2008: biology + HMM + histogram + training + files + plotting in the same object graph

2026: ObservationSource -> Hmm<Emission> -> HmmResult -> ResultProcessor
```

That is the contract future Lumrik biology should build against.
