# Module trajectories and plots

**Implementation:** `src/gene_model.rs::plot_frozen_modules`, `src/module_plots.rs::write_module_plots_named`.

Each plotted x-position is an **annotation category**, not a numerical time point. For each frozen nonzero module \\(j\\), category \\(k\\) and condition \\(h\\), the plotted ordinate is the mean across contributing cells of their uncorrected module-gene means \\(A_{cj}\\). Blue denotes reference/control; red denotes case. An absent condition is a gap, not a zero or an imputed reference.

The error bars are **approximately 95% normal confidence intervals for the cell-level mean**, \\(\bar A\pm1.96\,s/\sqrt n\\), not gene-to-gene dispersion or multiplicity-adjusted confidence intervals. When \\(n<2\\), the code cannot produce a finite SD/CI. Counts under the plot reflect observed cells in the respective category. Adjusted p-values and stars come from the separate Welch/Bonferroni procedure in [Statistics and limitations](statistical_analysis.md).

The x-axis annotation sequence is **cell counts → corrected p-value → category label → optional cell-type annotation**. `--label-angle 0` makes category labels horizontal; supported angles are 0–85 degrees. `--count-label-layout` supports `horizontal` and `rotated`.

Population ordering follows explicitly requested `--co-variant-clusters`, or relevant H5AD categorical/developmental metadata; a separate AI-suggested developmental order may generate `module_trajectories_developmental.tsv` and `module_*_trajectory_developmental.svg`. Lines between categories are visual guides; without a justified biological order, they do **not** establish a developmental trajectory. The analysis-wide correction in each output is computed from its modules and displayed categories; changing the set of categories can change the correction factor.

**Outputs:** `module_trajectories.tsv` and `module_<id>_trajectory.svg` (and optional `_developmental` counterparts). TSV rows are **per module × category × condition**; the same comparison's p-value and effect are therefore repeated in the two condition rows. `genes` is the number of module genes, while `n_cells` counts observations in that condition. See [Statistics and limitations](statistical_analysis.md) for the column definitions.
