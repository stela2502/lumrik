# Scientific overview

TraConMech begins with a nominated gene (for example TP53) and a biological comparison such as mutation-bearing versus wild-type cells. It searches the supplied STRING association network for partners of the nominated gene, connects partner names to regulatory gene sets from the configured Ommverse/MSigDB data, and assembles candidate target genes. These associations nominate hypotheses; they do **not** prove physical binding, direction of regulation, or causation.

For each candidate gene, TraConMech measures expression in case and control cells **within the same annotated population**. It then groups genes that show similar *patterns of relative expression change* across populations. A model created from those groups stores module membership and population-level reference centroids. Applying that model to another dataset uses the original modules without refitting them.

The principal trajectory figures answer a straightforward question: **within a given cell population, do cells in the two conditions differ in the average expression of the module's genes?** Each cell contributes one observation (the mean of its available module-gene log1p values). The figures show group means and approximate 95% confidence intervals. Two-sided Welch tests compare groups, and Bonferroni correction accounts for all module-by-population comparisons in that trajectory output.

A separate AddModuleScore-style calculation subtracts an expression-matched control-gene background from each cell's module-gene mean. That score is exported for interpretation and visualization; it is **not** the quantity tested in the trajectory figures.

The current inference is descriptive/exploratory at the **cell** level. Cells sampled from one patient or experiment are not independent biological replicates; corrected p-values therefore do not establish a patient-level effect or a rescue mechanism. Read [Statistics and limitations](statistical_analysis.md) before interpreting small effects.
