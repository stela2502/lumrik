# Module scoring: two different measurements

**Implementation:** `src/module_score.rs` (`score`); used by `src/gene_model.rs` and `src/module_plots.rs`.

Let \\(x_{gc}\\) be the supplied expression of gene \\(g\\) in cell \\(c\\), **already on the log1p scale**. For a frozen module \\(M_j\\), only genes present in the current dataset are used: \\(M'_j=M_j\cap G_{\rm observed}\\).

## Uncorrected module-gene mean

For each cell, TraConMech computes

\\[
A_{cj}=\frac{1}{|M'_j|}\sum_{g\in M'_j}x_{gc}.
\\]

Zeros contribute to the denominator. The software does **not** apply `ln_1p` a second time. This **uncorrected** per-cell mean is written to `gene_model_module_means.tsv` and is the input to the trajectory statistics. If a module has no observed genes, scoring fails.

## Expression-matched control background

To calculate a separate AddModuleScore-style score, the program computes each gene's mean across **all** cells (including zeros), ranks genes by this mean (breaking ties by gene index), and partitions the rank order into **24** approximately equal-sized bins.

For **each module gene**, it selects up to **100** genes from that gene's bin, sampled without replacement using an explicitly seeded deterministic xorshift-like PRNG. The selections are **unioned and deduplicated across module genes** to form a module-specific control set \\(B_j\\). The implementation does not explicitly exclude module genes from this control set; do not describe it as a disjoint negative control.

The background and corrected score are

\\[
B_{cj}=\frac{1}{|B_j|}\sum_{g\in B_j}x_{gc},\qquad
S_{cj}=A_{cj}-B_{cj}.
\\]

`gene_model_scores.tsv` contains \\(S_{cj}\\); `best_module` is simply the largest score in that cell. `gene_model_condition_scores.tsv` summarizes scores per condition and uses a population (divisor-\\(n\\)) SD. **The plotted trajectory means, confidence intervals and p-values use \\(A_{cj}\\), not \\(S_{cj}\\).**

The method is *Seurat-style*, not a claim of bit-for-bit equivalence with Seurat's `AddModuleScore`: exact binning, sampling PRNG, duplicate handling and gene universes can differ.
