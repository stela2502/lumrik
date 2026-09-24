# onehot_dna AI Contract

## Role

Fixed-length one-hot/IUPAC DNA representations for barcode, primer, and small-sequence matching.

## Two distinct semantics

Strict fixed barcode encoding and IUPAC biological-sequence compatibility are intentionally different. Do not collapse them.

Strict `OneHot<N>` treats non-ACGT as mismatch/zero. IUPAC-aware sequence representation preserves possibility masks, where compatibility is set intersection.

`OneHot<N>` uses four bits/base in `u128`; therefore `N <= 32` is structural.
