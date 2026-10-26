# Statistical analysis and limitations

**Implementation:** `src/module_plots.rs::{cell_stats,welch_test,p_label,write_module_plots_named}`. This chapter documents the code after the supplied CI and Bonferroni patches; it does not claim an independent validation study.

## Observation and estimand

For module \\(j\\), population/phase \\(k\\), and condition \\(h\\), each cell \\(c\\) contributes one **uncorrected** module mean \\(A_{cj}\\). The comparison estimates

\\[
\Delta_{jk}=\bar A_{jk,\mathrm{case}}-\bar A_{jk,\mathrm{control}}.
\\]

Thus a positive value indicates higher average **log1p expression** in case cells. The effect is *not* a log fold change of counts, and it is not a difference of background-corrected module scores.

## Confidence intervals

Within each condition, the code computes the sample SD using \\(n-1\\) and draws

\\[
\mathrm{CI}_{95\%}(\bar A)=\bar A\pm1.96\frac{s}{\sqrt n}.
\\]

These are normal-approximation **individual-group** intervals, not Welch confidence intervals for \\(\Delta\\), and not simultaneous Bonferroni-adjusted intervals. For fewer than two cells, SD and CI are unavailable.

## Two-sided Welch test

For case and control, let \\(\bar A_1,\bar A_0\\) be group means, \\(s_1,s_0\\) their sample SDs, and \\(n_1,n_0\\) their cell counts. The implementation calculates

\\[
V_1=\frac{s_1^2}{n_1},\quad V_0=\frac{s_0^2}{n_0},\qquad
 t=\frac{\bar A_1-\bar A_0}{\sqrt{V_1+V_0}},
\\]

\\[
\nu=\frac{(V_1+V_0)^2}{V_1^2/(n_1-1)+V_0^2/(n_0-1)},\qquad
p_{\rm raw}=2\Pr(T_\nu\ge |t|).
\\]

The survival probability is calculated with `statrs::distribution::StudentsT`. A test is unavailable (`NA`) if either group has fewer than two cells, the pooled standard-error variance is nonpositive/nonfinite, or the degrees of freedom are invalid. This is a **two-sided, unequal-variance, cell-level** test.

## Bonferroni family and annotations

For each trajectory output, the code counts nonzero module IDs \\(J\\) and plotted population/phase categories \\(K\\), then sets

\\[
m=J\,K,\qquad p_{\rm Bonferroni}=\min(1,\;m\,p_{\rm raw}).
\\]

**Untestable module × category combinations are still counted** in \\(m\\), conservatively. In an output with eight modules and three phases, \\(m=24\\). The factor is determined by that **output's** modules and displayed categories; it is not a hardcoded 24. If an optional developmental-order output contains a different set of categories, its family can differ. This implementation does **not** automatically combine separate analysis runs into a single multiplicity family.

SVG labels show the **adjusted** p-value in scientific notation, prefixed with `*` if \\(p_{\rm Bonferroni}<0.01\\), `**` if \\(p_{\rm Bonferroni}<0.001\\), and no star otherwise. `p=NA` means no valid test. The plot legend reports the factor and its module × category components. The thresholds are intentionally **not** the conventional 0.05/0.01 pair.

The TSV records `mean_log1p`, `sd`, `se`, `ci95_low`, `ci95_high`, `n_cells`, `case_n`, `control_n`, `delta_mut_minus_wt`, `welch_t`, `welch_df`, `p_raw`, `p_bonferroni`, and `bonferroni_factor`. The field `delta_mut_minus_wt` is named for the TP53 comparison, but mathematically means **configured case minus configured control**, whatever their labels.

## What the p-value cannot establish

1. **Pseudoreplication:** cells from one patient, library, culture or sample are not independent biological replicates. A small cell-level p-value is not a patient-level inference. Biological replication and hierarchical/pseudobulk analysis are needed for generalizable claims.
2. **Post-selection inference:** modules are discovered using observed case–control effects in the training data. Testing the *same* data with those selected modules can make nominal p-values anti-conservative. Bonferroni across plotted modules/phases **does not correct the earlier selection step**. Independent validation on held-out samples is preferable.
3. **Dependence and scale:** modules and phases can be correlated; Bonferroni remains conservative under arbitrary dependence **for a prespecified family of valid p-values**. Means of log1p values are not count fold changes.
4. **Small effects:** with many cells, statistically convincing differences can still be biologically negligible. Always inspect \\(\Delta\\), group means, coverage, sample origin and trajectory shape.
5. **No causal proof:** neither a STRING association nor a gene-set match establishes a TP53 regulatory mechanism or functional rescue.
