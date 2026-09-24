# sc-analysis AI Contract

## Role

Lightweight single-cell analysis primitives used by Lumrik: sparse input, QC, normalization, dimensional reduction, clustering, and reporting.

## Preserve raw meaning

Normalization/embedding/clustering must not overwrite raw counts or make transformed values appear to be molecule counts.

## Determinism

Where algorithms permit deterministic seeding/order, keep outputs reproducible for identical inputs and parameters.

## Scope

This crate analyzes an already-defined cell-by-feature dataset. It should not silently redefine upstream canonical-cell or molecule-provenance contracts.

## Initial population partitioning

Initial groups are deterministic geometric micro-patches in PCA space, not biological clusters.
Recursively split a patch at the median of its highest-variance informative PCA dimension until
each patch contains at most `max(50, ceil(1% of retained cells))` cells. Informative PCA
dimensions are the smallest leading set explaining at least 90% of the variance captured by
the computed PCs (minimum two dimensions when available). Biological grouping remains the
subsequent recomputed mean-expression Pearson merge; spatial patch IDs must not be presented
as final biological populations.

## Immune-receptor segment handling

Immunoglobulin/T-cell-receptor V, D, and J segment genes are receptor-identity features,
not ordinary transcriptional-state features. They MUST remain in raw input/count meaning,
and raw UMI cell filtering MUST still count them. They MUST NOT be eligible as variable
genes for PCA. When explicit pre-normalization VDJ exclusion is enabled, remove only V/D/J
segment genes from the surviving normalized expression matrix; constant-region genes such
as Ighm/Igha/Ighg*, Igkc/Iglc*, Trac/Trbc* and Jchain remain ordinary expression features.

## UMAP reporting

UMAP scatter panels are qualitative spatial views and should not draw coordinate axes.
Continuous overlays use a low-to-high yellow -> red -> blue scale with a numeric color
legend; integer-valued overlays display integer legend endpoints. Categorical UMAP
overlays retain an explicit category legend.

## Final-cluster reporting

Final biological cluster abundance is part of the interpretation and MUST be visible alongside
cluster identity. Report cell count and retained-cell fraction in cluster summaries and in the
final-cluster UMAP legend. The final-cluster legend belongs in a dedicated area outside the UMAP
coordinate system so it never obscures cells.

