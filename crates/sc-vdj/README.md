# sc-vdj 0.2.0

Posterior antigen-receptor analyzer for Lumrik/Nelrune.

`sc-vdj` is deliberately run **after** normal expression mapping. Its two required evidence sources are:

1. the retained Nelrune BAM, carrying the original cell/UMI identity; and
2. the complete per-cell Nelrune expression matrix.

The analyzer uses the same GTF + genome as the mapping run to build an explicit IG/TR germline reference with physical locus geometry.

## What it returns

For every cell:

- aggregate germline support from **all BAM reads/molecules for that cell**;
- best V, optional D, J and C identities and a compact `V(D)JC` notation;
- explicit recombination stage (`DJ`, `VJ`, or `VDJ`);
- V position as physical distance and normalized locus position, never inferred from gene nomenclature;
- locus-resolved germline/sterile-transcription candidates at configurable spatial resolution (default 64 bins) plus exact aligned intervals;
- RAG activity from RAG1/RAG2 expression;
- an explainable B- or T-lineage recombination-development score based on GEX marker evidence plus locus transcription;
- a reversible `evidence_code`: high bits contain confidence, low bits encode the selected program and marker-presence rationale.

For every sample, `write_reports()` **always** writes all seven receptor loci, including zero-only loci.

## Output files

- `vdj_cell_summary.tsv`
- `vdj_rearrangements.tsv`
- `vdj_sterile_spatial.tsv`
- `vdj_sterile_intervals.tsv`
- `vdj_development_rationale.tsv`
- `vdj_sample_summary.tsv`

The spatial sterile-transcription table is not limited to proximal/middle/distal labels. Those are derived summaries only. Every cell/locus carries a configurable binned profile, while exact aligned reference blocks are retained in the interval table. CIGAR `N` introns are excluded from covered blocks rather than being counted as transcriptional coverage.

## Heavy chain

IGH/TRB/TRD are represented as two-stage loci. A posterior call can therefore distinguish `DJ` evidence from completed `VDJ` evidence. IGK/IGL/TRA/TRG use the simpler `VJ` stage.

## Germline assignment

The fast k-mer index is only a candidate selector. Final segment support is rescored using local alignment and aggregated over every relevant BAM read for the cell. D segments fall back to exhaustive local scoring if short D segments seed poorly.

## Nelrune integration

The crate intentionally does **not** guess Nelrune's QNAME encoding or expression-matrix internals.

- Implement `BamIdentityResolver` using Nelrune's existing QNAME parser, then call `read_bam()`.
- Implement `ExpressionMatrix` directly for Nelrune/scdata's expression object. No matrix conversion is necessary.
- `LongTsvExpression` exists only as a convenient adapter for testing (`cell<TAB>gene<TAB>value`).

This keeps cell/UMI identity and expression semantics owned by Nelrune rather than duplicating them inside `sc-vdj`.

## Important interpretation boundary

The spatial BAM signal is stored as *sterile/germline-transcription evidence*. A read that contains a convincing V-J rearrangement is excluded from that evidence. However, partial rearranged reads can still align germline-like. The output therefore preserves raw bins/intervals so later chain-specific sterile-transcript rules (promoter/J-C splice structure, known sterile exons, allele-specific rules) can become stricter without reparsing the BAM.
