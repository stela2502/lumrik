# sc-regnet

**Trainable, biology-constrained regulatory-network inference for single-cell data.**

**Status:** design proposal, October 2026. No training implementation is claimed here.

`sc-regnet` is a proposed independent Rust crate in the Lumrik workspace. Its purpose is to learn which transcription-factor (TF) activities, TF partnerships, and regulatory-element routes are consistent with observed cell-specific expression responses. The model is **not** intended merely to reconstruct an expression matrix: it should produce ranked, testable regulatory hypotheses, preserve competing explanations, and identify observations it cannot explain.

## Motivation

TraConMech already identifies interesting response genes and candidate TF partnerships through its inspection workflow. The next question is not simply *which TFs could be involved?* but:

> Which combinations of TFs acting through which promoters and enhancers could generate the observed responses, and which predictions distinguish competing mechanisms?

Single-cell TF RNA detection is unreliable as a direct readout of TF activity. Therefore **TF expression must not be a required input**. Instead, infer latent TF and regulatory-element activities from the joint response of target genes, constrained by genomic annotations.

The problem is not a small gene-by-TF regression. Hundreds of responding genes can expand to thousands of transcripts and alternative promoters, thousands of candidate enhancers, overlapping TF-binding sites, and many-to-many enhancer–promoter links. Add approximately **300–10,000 cells**, and the task becomes structured learning over a large, heterogeneous regulatory graph. Those cells are repeated, noisy observations—not 10,000 independent perturbation experiments.

## Separation of responsibilities

| Component | Responsibility |
| --- | --- |
| **Ommverse** | Assembly-specific genes, transcripts, genomic coordinates, protein knowledge, and available reference regulatory evidence |
| **TraConMech** | Select response genes, cell states/contrasts, construct biological candidate graphs, interpret and report hypotheses |
| **sc-regnet** | Graph representation, latent-state inference, training, prediction, model comparison, uncertainty and perturbation simulations |
| **sc-beacon** | Cell calling/QC and other existing responsibilities; not the regulatory training engine |

`sc-regnet` must not depend on Ommverse or TraConMech types. They export a versioned, typed training graph and observation dataset. All genomic evidence must retain assembly, coordinates, source, provenance and confidence where available.

## Biological graph

A heterogeneous, directed, typed graph may contain:

- TF proteins and candidate complexes; partner/cofactor and competition relations.
- Motif matches and experimentally observed binding sites (distinguished explicitly).
- Promoters, transcription start sites, distal enhancers and other regulatory elements.
- Enhancer–promoter links with evidence and tissue/context annotations.
- Transcripts and genes, including alternative promoter usage and transcript-to-gene aggregation.

A binding site is **not** evidence of functional regulation by itself. An enhancer–gene link is **not** automatically active in the studied cell type. Activating/repressing signs may be unknown. Missing edges must not be treated as proven negatives.

### Existing Lumrik sources and limitations

The supplied October 15 Lumrik source archive describes Ommverse's FANTOM5 enhancer annotations and an **ENCODE4 merged candidate binding-region union**. It explicitly states that the full ENCODE4 TF rPeak BigBed is kept in the source cache rather than embedded wholesale in the `.ommverse` index. Consequently, **TF-specific binding evidence must not be assumed to be directly queryable from the current compact index**. Implement an evidence extraction/import layer before claiming TF-specific binding-site models. Validate enhancer–promoter associations, genomic assembly consistency, cell-type relevance and provenance separately.

## Inputs and outputs

**Inputs**

1. Sparse cell-by-gene counts or expression measurements, with cell metadata, conditions, batch/donor identifiers and optional developmental state or pseudotime.
2. A selected response-gene set plus the full expressed-gene background, including non-responders for calibration.
3. A regulatory graph constructed from genomic evidence, including optional transcript-level observations when available.
4. Optional perturbation labels, allele/mutation status, protein or chromatin measurements; **none required for the initial baseline**.

**Outputs**

- Predicted expression distributions and held-out prediction metrics, by gene and cell/state.
- Inferred latent TF/program activities with uncertainty and clear identifiability warnings.
- Ranked TF partnerships, regulatory paths, candidate promoters and enhancers, with evidence and contradictions.
- Alternative mechanisms that fit similarly well; residual genes and unexplained states.
- Model-based virtual perturbation predictions, explicitly **not asserted as causal effects** without suitable interventions.
- Reproducible model weights, graph/evidence versions, configuration, seeds, training diagnostics and evaluation splits.

## Modeling strategy: simple baseline first, trainable graph second

### Stage A — falsifiable baseline

Aggregate cells by biologically justified states/conditions (or use a count-aware hierarchical model), and fit sparse, weighted, graph-constrained regression to response measurements. Include uncertainty, correlated TF target sets, matched expressed-gene controls and null models. This establishes whether genomic priors actually improve out-of-sample predictions over expression-only and shuffled-graph baselines.

### Stage B — cell-specific latent regulatory model

Use **one shared regulatory graph** across all cells. Give each cell a compact latent regulatory state. Learn shared parameters connecting TF/program activity to regulatory-element state, promoter/transcript output and gene-level expression. Use an appropriate count likelihood (e.g., negative binomial with library-size offsets) for raw UMI counts, or a justified alternative for transformed/pseudobulk data. Model technical batch and biological state separately when the design permits.

A conceptual generative route is:

```text
Genomic priors + TF/cofactor relationships
                  |
        Shared regulatory graph
                  |
   Cell-specific latent TF/program state
                  |
  Enhancer/promoter activity (mostly latent)
                  |
      Transcript initiation / abundance
                  |
      Observed gene counts per cell
```

Where transcript observations are absent, transcript/promoter paths may be non-identifiable from gene-level counts. Retain uncertainty or collapse indistinguishable paths; do not pretend that 10x gene counts resolve alternative promoters.

### Stage C — heterogeneous graph learning

Only after Stage B beats baselines, consider typed message passing, attention or other graph-neural-network components. Constrain parameter growth using shared edge-type transformations, sparsity, evidence-based priors and regularization. Do not create a freely trainable parameter for every cell × binding site × enhancer × transcript combination.

### Stage D — simulations and active experimental design

Perturb a candidate TF/program or regulatory edge *in the fitted model*, predict expression changes and rank experiments that could distinguish otherwise equivalent explanations. These are **hypotheses**, not experimentally validated causal predictions.

## How to train it in Rust

**TensorFlow is not required.** We need differentiable tensor operations, automatic differentiation, optimizers, sparse/graph operations and checkpointing. Evaluate a Rust-native training backend such as **Burn** first; choose the backend only after verifying current support for the required sparse operations, autodiff, hardware and deployment targets. A small baseline can use ordinary Rust numerical/statistical crates without deep-learning machinery. GPU acceleration is an optimization, not a prerequisite for the first experiment.

Proposed training loop:

1. Freeze and version the evidence graph and observation schema.
2. Split data by **donor, perturbation, experiment or other independent biological unit** when available. Held-out cells from the same donor are not sufficient evidence of generalization.
3. Initialize biologically constrained parameters and per-cell latent states.
4. Predict observed expression or response; calculate a count-aware likelihood plus complexity/graph-prior penalties.
5. Backpropagate through trainable operations, update parameters, checkpoint, and monitor validation performance.
6. Compare against simple baselines and graph-randomization/null controls.
7. Export ranked hypotheses, calibrated uncertainty and held-out diagnostics.

Sketch of the objective:

\[
\mathcal{L} = -\log p(Y\mid G,Z,\theta) + \lambda_1\Omega_{\mathrm{sparse}}(\theta) + \lambda_2\Omega_{\mathrm{prior}}(\theta,G)
\]

Here `Y` is observed expression, `G` the evidence graph, `Z` latent cell states and `θ` shared model parameters. This is a **design equation**, not an implemented loss.

## First-month milestone

| Week | Deliverable | Exit criterion |
| --- | --- | --- |
| 1 | Define typed graph/observation formats; synthetic fixtures; reference-data export audit | Known simulated TF/program activities are recoverable when identifiable |
| 2 | Sparse baseline, count-aware or pseudobulk response model, null controls | Honest held-out metrics and meaningful uncertainty |
| 3 | First trainable cell-specific latent model; CPU training; checkpoints | Training reproduces synthetic responses and beats trivial baselines |
| 4 | TraConMech integration and TP53/erythroid case study | Ranked mechanisms, residuals, evidence trails and explicit limitations |

**Scope discipline:** A one-month milestone is a working, testable inference prototype, **not** a validated causal regulatory foundation model. ENCODE/FANTOM ingestion gaps and data identifiability may be the primary bottlenecks rather than compute.

## Suggested crate layout (proposal)

```text
crates/sc-regnet/
├── Cargo.toml
├── README.md
└── src/
    ├── lib.rs
    ├── graph.rs          # typed regulatory graph, edge evidence
    ├── observations.rs   # sparse expression, metadata, splits
    ├── baseline.rs       # simple interpretable comparison model
    ├── latent.rs         # shared regulatory model, cell latent states
    ├── train.rs          # optimization and checkpoints
    ├── evaluate.rs       # held-out tests, null models, calibration
    ├── simulate.rs       # explicitly hypothetical perturbations
    └── export.rs         # versioned results and provenance
```

The CLI belongs initially in TraConMech, with a narrowly named command such as `traconmech regulatory-model` rather than an overly broad `explain`. Final command names and argument contracts remain open.

## Scientific acceptance criteria

A credible model must (1) beat simple baselines on independent held-out biological units; (2) distinguish evidence-backed edges from speculative links; (3) report when promoter/enhancer or TF alternatives cannot be resolved; (4) avoid leaking cell labels, conditions or target responses through graph construction; (5) produce reproducible predictions; and (6) nominate experimentally discriminating tests rather than overclaim causal mechanisms.

**Long-term aim:** move from observed single-cell responses through genomic regulatory architecture to experimentally testable TF-complex, promoter and enhancer hypotheses—while keeping the mathematics reusable outside TraConMech.
