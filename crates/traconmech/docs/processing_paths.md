# Processing paths and outputs

This is an implementation map, not a promise that every architectural feature described in the crate README is available.

| CLI path | Source entry point | Implemented processing | Principal outputs |
| --- | --- | --- | --- |
| `gene-model create` | `main.rs` → `experiment::resolve` → `regulatory::run` → `gene_model::create` | Candidate-gene discovery, per-population effects, optional effect modules, frozen centroids | `experiment.yaml`, `effect_clusters.tsv`, `module_membership.tsv`, portable model JSON, training-cluster table |
| `gene-model apply` | `main.rs` → `gene_model::validate` → optional `regulatory::run` → `gene_model::apply` → `model_umap::write` → `gene_model::plot_frozen_modules` | Frozen scoring and plots; optional fresh investigation **without** refitting modules | `gene_model_scores.tsv`, `gene_model_module_means.tsv`, coverage files, module UMAPs, trajectories |
| `investigate` | `main.rs` → `regulatory::run` | STRING/MSigDB hypothesis investigation; effect-module generation only when enabled | Partner and candidate evidence, per-population effects, optional modules/plots |
| `analysis` | `main.rs` → `sc_analysis::cli::run_analysis` | Delegated single-cell analysis | Depends on `sc_analysis` CLI |
| `match-wt` | `main.rs` → `wt_match::run` | Observed within-population WT comparison | WT-matching reports |
| `profile-clusters` | `main.rs` → `cluster_profiles::run` | Cluster expression/detection summaries and correlations | Cluster profile/correlation tables |
| `genomic` | `main.rs` (`run_lights`, HMM functions) | Genomic-bin/light-state and HMM processing | Genomic/HMM result tables |
| `ai`, `inspect`, `inspect-results`, `contracts`, `scope` | `main.rs`, `ai.rs`, `inspect_model.rs`, `contracts.rs` | Evidence/inspection/registry utilities | Mode-dependent |

**Scope of detailed mathematical review:** this handbook audits the regulatory → module-discovery → frozen-model → module-scoring → trajectory/statistics path. The other CLI modes are identified above but their algorithms are **not** presented as fully audited here. In particular, the genomic HMM is a separate statistical procedure and should receive its own methods chapter before publication.

## Output interpretation

- `module_membership.tsv`: training gene → exclusive module assignment (including possible module 0).
- `gene_model_scores.tsv`: per-cell background-corrected module score, plus highest-scoring module.
- `gene_model_module_means.tsv`: per-cell **uncorrected** module mean on the supplied log1p scale.
- `module_trajectories.tsv`: condition-wise summary of those uncorrected means; one Welch comparison per module × category, repeated on the two condition rows.
- `module_*_trajectory.svg`: means, approximate individual-group 95% CIs, cell counts and adjusted p-values.
- `experiment.yaml`: analysis-design decisions; retain with results.

## Future mdBook site

The directory is intentionally laid out as a book source. From the crate root, after installing `mdbook`, the handbook can be built with `mdbook build docs` using the provided `docs/book.toml`. Equations are written in LaTeX notation; **stock mdBook does not automatically typeset every LaTeX block**. Choose and configure a math renderer (such as KaTeX) before publishing a public site, and verify the rendered equations. This is documentation scaffolding, not a claim that a site has been built.
