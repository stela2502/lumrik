# sc-te AI Contract

## Role

`sc-te` performs single-cell transposable-element quantification from multimapping alignments.

## Multimapping is signal

TE reads are expected to multimatch. Do not apply unique-mapping assumptions from ordinary gene quantification unless explicitly justified.

## Cell identity

Preserve upstream cell/UMI identity and deduplication semantics. TE assignment policy may differ from GEX gene assignment, but must not redefine the cell barcode itself.

## Reference provenance

Repeat annotation/reference version is part of the interpretation of TE counts and should remain traceable in outputs.
