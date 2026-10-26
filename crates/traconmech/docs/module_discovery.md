# Investigation and module discovery

**Implementation:** `src/regulatory.rs` (`run`, `stats`, `effect_modules`, `stability_weight`); `src/experiment.rs`.

## Candidate genes

The nominated gene is resolved against the supplied STRING protein-info table. Direct incident edges in the supplied STRING links table are retained when `combined_score >= --string-min-score`; a headered full links file or headerless standard three-column file is supported. Protein identifiers are mapped to preferred gene names. Partner symbols are matched to configured regulatory gene-set names by underscore-delimited symbol tokens. The union of the matching sets' genes forms the candidate list. A gene may occur in several sets, but this is **not** a statistical enrichment test.

## Within-population comparisons

For gene \\(g\\), population \\(k\\), and condition \\(h\\), let \\(x_{gc}\\) denote the **supplied** expression-matrix value (expected to be log1p-transformed). Define

\\[
\mu_{gkh}=\frac{1}{n_{kh}}\sum_{c\in C_{kh}}x_{gc},\qquad
r_{gk}=\log_2\left(\frac{\mu_{gk,\mathrm{case}}+0.01}{\mu_{gk,\mathrm{control}}+0.01}\right).
\\]

**Important:** \\(r_{gk}\\) is the log2 ratio of **means of log1p values**; it is *not* a log2 fold change of raw counts. The pseudocount `0.01` is applied to those means. The code also reports detection fraction (nonzero values), and a population-standard-deviation with divisor \\(n\\) for the gene-level summaries.

Only populations with both observed conditions enter these effects. The implementation explicitly does **not** substitute WT cells from a different population. `--co-variant-clusters` can select and order a specific sequence, and fails if a selected population lacks a comparison.

## Discovering modules from effect profiles

For each gene, form the effect vector \\(\mathbf r_g=(r_{g1},\ldots,r_{gK})\\). Center it and normalize its Euclidean length:

\\[
\mathbf z_g=\frac{\mathbf r_g-\bar r_g\mathbf 1}
{\lVert\mathbf r_g-\bar r_g\mathbf 1\rVert_2}.
\\]

Genes with fewer than two entries or centered norm below `1e-9` enter **module 0** (flat/uninformative trajectories). Remaining genes undergo deterministic spherical k-means-like clustering: at most **eight** initial centers, farthest-first initialization, cosine/dot-product assignment, normalized centroid updates, at most **40** iterations, deterministic gene-name iteration and tie resolution. Modules are numbered from 1. **Eight is a cap, not a guarantee of eight modules.** Module IDs are labels, not effect magnitudes.

A separate *stability weight* used in reporting/ranking is

\\[
D=\tfrac12\left(\frac{\sigma_{\rm case}^2}{\mu_{\rm case}+0.1}
+\frac{\sigma_{\rm control}^2}{\mu_{\rm control}+0.1}\right),\qquad
w=\frac{1}{1+D}.
\\]

The code sets \\(w=0\\) when the maximum condition mean **or** detection fraction is below `0.05`. These weights rank/report gene evidence; they do **not** define the module assignments, which use the centered effect profiles above.

**Key files:** `effect_clusters.tsv`, `module_membership.tsv`, `effect_modules.tsv`, `module_measurements.tsv`, `stable_module_genes.tsv`, `wt_reference_clusters.tsv`, plus STRING/partner evidence tables. Module 0 remains auditable but is not plotted as an ordinary module.
