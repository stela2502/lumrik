# sc-analysis

`sc-analysis` is Lumrik’s Rust crate for single-cell gene expression analysis.
Its `sc-analyse gex` command takes raw UMI matrices from one or more samples
through QC, PCA, graph clustering, marker analysis, and SVG reporting. It
exports an `analysis.h5ad` with normalized expression and cell annotations.

## How the analysis works

1. **Load one or more samples.** Read 10x-style Matrix Market (MEX) folders,
   Norn sample folders with an exonic matrix, or H5AD files with raw integer
   counts in `X`. For multiple samples, align genes by loaded feature name
   across the union of features and prefix barcodes with the sample name.
   ([Input and sample merging](src/bioinformatics.rs), [MEX reader](src/mex.rs))
2. **Filter features and cells.** Apply an optional `--target-gene-list` first.
   Exclude mitochondrial (`MT-`) and ribosomal (`RPS`/`RPL`) genes from
   downstream expression. Add `--exclude-vdj-before-normalization` to also
   exclude immunoglobulin and T-cell receptor V, D, and J segments while
   retaining constant-region genes. The UMI threshold uses counts from the
   surviving genes. IG/TR V(D)J segments are always excluded from PCA feature
   selection. ([QC](src/qc.rs), [Filtering and PCA preparation](src/bioinformatics.rs))
3. **Equalize UMI depth for PCA.** Keep cells with at least
   `--min-umi-count` surviving UMIs (800 by default). Rescale their counts to
   that depth with seeded stochastic rounding and an exact-total correction.
   Apply `log1p`, select the highest-variance eligible genes, standardize each
   selected gene, and fit PCA with Linfa.
   ([Resampling](src/normalize.rs), [PCA](src/pca.rs))
4. **Scale expression for interpretation.** Separately compute
   `log1p(10,000 × gene UMIs / surviving cell UMIs)` from the surviving
   original counts. Use this matrix for expression export, marker statistics,
   expression overlays, and Pearson merging. ([Normalization](src/normalize.rs))
5. **Partition PCA space.** Recursively split cells at the median of the
   locally highest-variance retained PC. Splits consider the leading PCs
   covering 90% of captured variance, with at least two when available. By
   default, each patch holds at most the larger of 50 cells or 1% of retained
   cells rounded up. These patches index the neighbor search.
   ([Spatial partitioning](src/partition.rs))
6. **Build a sparse neighbor graph.** Search the patch hierarchy in parallel
   for each cell’s exact Euclidean nearest neighbors in PCA space. Bounding
   boxes prune branches that cannot contain a closer cell, and neighbors can
   span patch boundaries. Symmetrize the links and weight them by distance
   relative to each cell’s local neighbor scale. ([Neighbor search](src/leiden.rs))
7. **Find Leiden communities.** Run seeded Leiden with the RB configuration
   quality function on the weighted PCA graph. Tune `--leiden-neighbors` and
   `--leiden-resolution` to control the graph and community granularity.
   ([Leiden](src/leiden.rs))
8. **Merge similar communities.** Merge the eligible pair with the highest
   Pearson correlation between mean profiles across all surviving genes.
   Recompute the merged profile before choosing the next pair. Stop at the
   correlation cutoff (`--merge-pearson`) or request a cluster count with
   `--target-clusters`. Both modes honor `--max-cluster-fraction`. The merge
   history and cell-level crosswalk record the transitions from PCA patch to
   initial and final community. ([Pearson merging](src/cluster.rs))
9. **Report cluster statistics.** Compare each final cluster with the rest
   using Wilcoxon rank-sum tests. Report adjusted p-values, fold changes,
   detection fractions, and expression quantiles. The output also includes
   cluster sizes, UMI distributions, and cluster-quality diagnostics.
   ([Statistics](src/stats.rs), [Reports](src/report.rs))
10. **Plot the results.** Compute a 2-D UMAP from the PCA coordinates and write
    SVGs for initial and merged Leiden labels, sample membership when available,
    UMI depth, and QC fractions. Also export PCA explained variance and UMI
    rank plots. ([UMAP](src/embedding.rs), [Plots](src/report.rs))

## How to get it

Follow the [Lumrik build instructions](../../README.md#building) for the
workspace setup. From the root of your Lumrik checkout, build the CLI:

```sh
cargo build --release -p sc-analysis --bin sc-analyse
```

Run a joint analysis using a MEX folder and an H5AD count matrix:

```sh
./target/release/sc-analyse gex \
  --data sample-a/mex sample-b.h5ad \
  --sample-names sample-a sample-b \
  --out sc_analysis
```

For one input, pass its path to `--data` and omit `--sample-names`.
For all options, run `./target/release/sc-analyse gex --help`.
