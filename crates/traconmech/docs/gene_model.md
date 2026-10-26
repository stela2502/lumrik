# Frozen gene models

**Implementation:** `src/main.rs` (`Mode::GeneModel`), `src/gene_model.rs`, `src/experiment.rs`.

## `gene-model create`

The CLI loads an H5AD file or Matrix Market directory, resolves the experimental design, and runs the regulatory investigation with effect-module discovery enabled by `--effect-cluster`. `create` then reads `module_membership.tsv` from the output directory; model creation requires nonempty modules. The portable JSON stores the nominated gene, training grouping name, feature names, frozen gene-to-module assignments, and one expression centroid per training population:

\\[
\mathbf m_k(g)=\frac{1}{|C_k|}\sum_{c\in C_k}x_{gc}.
\\]

These are **means of supplied log1p values**, not re-log-transformed or back-transformed counts. The model format is `traconmech.gene-model.centroids`, version `1`. The saved centroids are reference metadata; the subsequent module-scoring implementation does **not** classify external cells by nearest centroid.

## `gene-model apply`

The CLI validates the model's gene identity and loads the external data. At least **10** shared features and at least **50%** overlap with all model features are required. The model's original module memberships are used unchanged. `apply` writes per-cell background-corrected scores, uncorrected module-gene means, condition summaries and feature coverage. It also generates module UMAPs when the input is an H5AD file, and frozen-module trajectory plots using the original model gene sets.

**Crucial processing-path detail:** `gene-model apply` can additionally run a fresh **regulatory investigation** on the external data when both comparison conditions are present, but sets `effect_cluster=false` for that path. That investigation is not the same as refitting the saved model. When a condition is absent, the differential regulatory investigation is skipped and frozen-model outputs still run.

## Experimental configuration and reproducibility

`experiment::resolve` uses condition annotations and explicit or resolved case/control labels, storing decisions in `experiment.yaml`. Keep that file with the model, input provenance, reference versions and command invocation. **Never delete or silently replace it during reruns.** Some generated outputs are rewritten when commands run; protect run directories and archive input configuration before regeneration.

**Key files:** portable model JSON, `gene_model_training_clusters.tsv`, `gene_model_scores.tsv`, `gene_model_module_means.tsv`, `gene_model_condition_scores.tsv`, `gene_model_module_coverage.tsv`, `gene_model_coverage.json`, `frozen_module_coverage.tsv`, `module_trajectories.tsv`, `module_*_trajectory.svg`.
