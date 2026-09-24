# scdata AI Contract

## Role

`scdata` stores sparse single-cell feature/UMI observations and exports deterministic matrices. It is a data structure/selection layer, not the authority for primer chemistry or splice biology.

## UMI semantics

UMI counts are molecule counts after the caller's intended deduplication semantics. Do not convert a threshold contract from `>= N` to `> N`, or vice versa, without an explicit API/contract change and tests at the boundary value.

## Filtering

Filtering cells and filtering features are distinct operations. Selecting canonical cells MUST NOT accidentally erase the unfiltered feature universe or silently redefine which observations are biologically valid.

Unfiltered outputs must remain genuinely unfiltered with respect to canonical-cell selection. Filtered outputs must be derived from the explicit selected-cell set.

## Determinism

Sparse matrix export, barcode ordering, feature ordering, and identifiers should remain deterministic for identical input. Do not introduce hash-iteration-dependent output ordering.

## Separation of concerns

`scdata` should not infer GEX versus VDJ provenance from feature names. Provenance must be decided upstream and only appropriate observations inserted/selected.
