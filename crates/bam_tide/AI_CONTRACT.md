# Lumrik AI Engineering Contract

This file documents architectural invariants that MUST be preserved when
modifying Lumrik. Do not infer alternative semantics from implementation
details. If code appears to violate this contract, investigate before changing
the contract.

## Molecule origin is biological state

Mapped BAM QNAMEs retain molecule origin.

- `G` = gene-expression (GEX) molecule
- `V` = VDJ/custom-capture molecule

These classes MUST NOT be collapsed merely because both map to the genome.

### Canonical GEX cell calling

Canonical GEX cells are defined from GEX (`G`) molecules ONLY.

VDJ (`V`) molecules:
- MUST NOT be quantified by `bam_tide` GEX quantification
- MUST NOT contribute to GEX expression matrices or UMI totals
- MUST NOT contribute to the barcode-rank distribution used for GEX cell calling
- MUST NOT influence the GEX knee/cell threshold
- are quantified by the dedicated `sc-vdj` crate

The final canonical-cell set derived from GEX is subsequently used to filter
the GEX quantified data.

## Splice classification

`MatchClass` values have distinct biological meanings.

- `ExactJunctionChain`: exonic GEX evidence
- `Compatible`: exonic GEX evidence
- `JunctionMismatch`: exonic GEX evidence; observed splice structure is absent
  from the transcript annotation
- `Intronic`: POSITIVE evidence that sequence occupies an annotated intron
- `Incompatible`: transcript model cannot explain the alignment; NOT synonymous
  with intronic
- `OverhangTooLarge`: outside accepted transcript-boundary tolerance

Failure to fit an exon model MUST NOT automatically become `Intronic`.

Transcript-end annotations are imprecise. The effective terminal-overhang
tolerance MUST NOT be below 100 bp.

## Quantification versus cell calling

Do not confuse these operations.

Quantification may retain evidence that is intentionally excluded from the
population used to CALL canonical GEX cells.

Cell calling is a biological selection operation, not simply a summary of every
molecule present in the expression matrix.

## Health server

The health server is part of the correctness/debugging contract, not cosmetic UI.

During BAM quantification it should expose meaningful quantification state,
including:
- BAM records processed
- explicit GEX (`G`) provenance count; the existing `Candidate molecules` server
  slot is used for this during BAM-only quantification and MUST exclude `V` and
  legacy/undecodable records
- canonical GEX cells once called
- exonic/intronic cells, genes and UMIs
- transcript match-class counts

During BAM-only nelrune quant, the health server should display only metrics that this execution actually measures. FASTQ processing/routing metrics belong to full nelrune runs and should be omitted or "NA" rather than displayed as zero.

## General rule

Do not "simplify" biologically distinct states merely because they currently
share a data structure.

If an invariant in this file appears wrong, STOP and ask before changing it

## Splice-mismatch diagnostics

When diagnosing splice matching, preserve signed donor and acceptor displacement
relative to the nearest annotated junction. Do not reduce all junction mismatches
to a single counter: small systematic +/-1 or +/-2 bp shifts are diagnostically
different from genuinely novel junctions. Use bounded histograms rather than
storing per-read offsets.
