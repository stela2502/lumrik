# sc_beacon

Ambient-aware CRISPR guide calling from 10x Feature Barcode matrices.

`sc_beacon` distinguishes genuine guide RNA (gRNA) signal from ambient guide contamination using the information already present in the **raw** and **filtered** 10x feature-barcode matrices.

The central idea is simple:

> Use droplets rejected by the 10x cell-calling step to learn the ambient guide distribution, then ask for every guide observed in every real cell whether its count is better explained by ambient contamination or genuine guide expression.

Importantly, `sc_beacon` does **not** use a winner-takes-all assignment. A cell may contain zero, one, two, or multiple genuine guides.

## Input

The program expects the standard 10x output pair:

```text
raw_feature_bc_matrix/
├── matrix.mtx.gz
├── barcodes.tsv.gz
└── features.tsv.gz

filtered_feature_bc_matrix/
├── matrix.mtx.gz
├── barcodes.tsv.gz
└── features.tsv.gz
```

By default, features with the 10x feature type

```text
CRISPR Guide Capture
```

are analyzed.

The feature type can be changed with `--feature-type`.

## How it works

The analysis consists of three stages.

### 1. Learn the ambient guide composition

Barcodes occurring in the raw matrix but **not** in the filtered matrix are treated as background droplets.

These droplets provide an empirical measurement of the relative abundance of every guide in the ambient pool.

For guide `g`:

```text
p_g = fraction of ambient guide molecules belonging to guide g
```

A guide that is very abundant in the ambient pool therefore requires stronger evidence before it is called genuine in a cell than a guide that is rarely observed in the background.

The ambient guide composition is estimated before fitting the model to the filtered cells.

### 2. Model guide counts in real cells

For each filtered cell `c`, `sc_beacon` estimates a cell-specific ambient guide burden:

```text
lambda_c
```

The expected ambient count for guide `g` in cell `c` is then

```text
lambda_c * p_g
```

and the ambient component is modeled as

```text
A_cg ~ Poisson(lambda_c * p_g)
```

Genuine guide expression is modeled independently for each guide using a negative-binomial distribution:

```text
T_cg ~ NegativeBinomial(mu_g, theta_g)
```

where:

- `mu_g` is the guide-specific mean genuine expression;
- `theta_g` describes guide-specific overdispersion.

A genuinely guide-positive cell can still contain ambient molecules of the same guide.

Therefore the genuine-guide model is not a replacement for the ambient model. The observed count is modeled as

```text
X_cg = A_cg + T_cg
```

and the likelihood of the genuine-guide state is calculated from the Poisson/negative-binomial convolution.

For every observed `(cell, guide)` pair the model estimates

```text
P(genuine guide | observed count)
```

### 3. Call guides independently

Each guide in each cell is evaluated independently.

There is deliberately **no winner-takes-all step**.

For example:

```text
cell_1
    guide_A    genuine
    guide_B    ambient

cell_2
    guide_A    genuine
    guide_C    genuine

cell_3
    no genuine guide
```

The second cell is retained as a genuine two-guide cell rather than being forced to choose between `guide_A` and `guide_C`.

This makes the caller suitable for experiments in which multiple perturbations per cell are possible or expected.

## Statistical calling

For every non-zero cell-guide observation, `sc_beacon` reports both a model posterior and an ambient-only significance test.

The ambient-only hypothesis is

```text
X_cg ~ Poisson(lambda_c * p_g)
```

from which an upper-tail probability is calculated:

```text
P(X >= observed count | ambient)
```

These p-values are corrected across observations using the Benjamini-Hochberg procedure.

By default, a guide is called genuine when both:

```text
posterior_real >= 0.95
q_value        <= 0.05
```

are satisfied.

This provides two complementary pieces of evidence:

1. the mixture model strongly favors genuine guide expression;
2. the observed count is inconsistent with the empirically measured ambient model.

## Biological versus numerical convergence

The mixture model iteratively estimates parameters including:

```text
lambda_c
mu_g
theta_g
prior_real_g
```

Numerical parameters can continue to change slightly even after the inferred biological assignments have stopped changing.

For this reason, `sc_beacon` distinguishes between:

### Mathematical convergence

The fitted model parameters no longer change beyond the configured numerical tolerance.

### Assignment convergence

The inferred cell-guide assignments reproduce themselves across consecutive iterations.

The latter is particularly relevant for this application: the purpose of the model is to determine which guides are genuinely present in which cells, rather than to optimize nuisance parameters to arbitrary numerical precision after those assignments have become stable.

Several consecutive stable iterations are required before assignment convergence is accepted.

This also provides a useful diagnostic distinction between:

```text
model parameters are still refining
```

and

```text
the inferred biological result is still changing
```

## Sparse data representation

`sc_beacon` uses [`scdata`](https://github.com/stela2502/scdata) for its cell-major sparse representation.

During MatrixMarket import, a complementary guide-major representation is constructed:

```rust
guide -> Vec<GuideObservation>
```

This provides both views required by the model:

```text
cell -> guide counts
guide -> cell observations
```

without converting the input into an intermediate dense matrix.

Only the requested guide features are retained.

## Installation

Clone the repository and build the release binary:

```bash
git clone https://github.com/stela2502/sc_beacon.git
cd sc_beacon

cargo build --release
```

The executable is then available as:

```bash
target/release/lumrik-guides
```

Alternatively:

```bash
cargo install --path .
```

## Usage

```bash
lumrik-guides \
    --raw /path/to/raw_feature_bc_matrix \
    --filtered /path/to/filtered_feature_bc_matrix \
    --out guide_calls
```

For example:

```bash
lumrik-guides \
    --raw outs/raw_feature_bc_matrix \
    --filtered outs/filtered_feature_bc_matrix \
    --out outs/lumrik_guide_calls
```

Useful options include:

```text
--feature-type
--threads
--ambient-alpha
--posterior
--fdr
--max-iterations
```

Run

```bash
lumrik-guides --help
```

for the current command-line interface.

## Output

The output directory contains several complementary reports.

### `ambient_guides.tsv`

The empirical ambient guide composition.

Typical columns include:

```text
guide_id
guide_name
ambient_umis
p_g
```

This file shows how strongly each guide is represented in the background pool.

### `guide_models.tsv`

The fitted genuine-expression model for each guide.

Typical parameters include:

```text
guide_id
guide_name
prior_real
true_mean
theta
```

### `guide_calls.tsv`

The detailed cell-guide calls.

For every observed cell-guide pair this includes information such as:

```text
barcode
guide_id
guide_name
umi_count
lambda_c
p_g
expected_ambient
posterior_real
ambient_p
q_value
called
```

This file deliberately retains rejected observations as well as positive calls so that the classification can be inspected rather than providing only a final assignment.

## QC statistics

The program reports simple assignment statistics after calling, including:

```text
cells_total
cells_no_guide
cells_single_guide
cells_multi_guide
cells_2_guides
cells_3_guides
cells_4plus_guides
called_guides_total
```

For example:

```text
Error Type              Count
cells_total             4,776
cells_no_guide            872
cells_single_guide      3,899
cells_multi_guide           5
cells_2_guides              5
called_guides_total      3,909
```

These statistics are generated using [`mapping_info`](https://github.com/stela2502/mapping_info).

## Why use the raw matrix?

The raw 10x matrix contains a large population of droplets that did not pass cell calling.

For CRISPR guide data these droplets are useful rather than simply being discarded: they provide an empirical measurement of the guide molecules present outside called cells.

Consequently, `sc_beacon` does not need to infer the ambient guide composition solely from the cells it is trying to classify.

Conceptually:

```text
raw-only droplets
        |
        v
ambient guide composition p_g
        |
        +-------------------+
                            |
filtered cells              |
        |                   |
        v                   v
cell-specific lambda_c + guide-specific background
        |
        v
ambient + genuine-guide mixture
        |
        v
P(genuine | count)
        |
        v
independent multi-guide calls
```

## Design goals

`sc_beacon` is intended to be:

- **ambient-aware** — background guide abundance is measured from the experiment;
- **cell-aware** — ambient burden is allowed to differ between cells;
- **guide-aware** — different guides may have different genuine-expression distributions;
- **multi-guide aware** — cells are never forced into a single-guide assignment;
- **transparent** — posterior probabilities, ambient expectations, p-values and q-values are retained;
- **sparse** — guide matrices are processed without unnecessary dense representations;
- **standalone** — no runtime model compilation or external statistical environment is required.

## Current status

`sc_beacon` is under active development.

The statistical model and calling strategy are functional, but the implementation and model validation are still being developed. In particular, performance optimization and systematic benchmarking against simulated and experimental CRISPR datasets are ongoing.

Results should therefore be validated carefully before using the software as the sole basis for biological conclusions.

## Lumrik

`sc_beacon` is part of the developing **Lumrik** bioinformatics ecosystem, with an emphasis on efficient, transparent and reusable analysis tools implemented in Rust.

Related projects:

- [`scdata`](https://github.com/stela2502/scdata) — sparse single-cell data structures
- [`mapping_info`](https://github.com/stela2502/mapping_info) — analysis and processing statistics

## License

AGPL-3.0-or-later.

See the repository license for details.