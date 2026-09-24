# snp-index AI Contract

## Role

`snp-index` indexes known variants and matches read evidence to loci with single-cell integration.

## Evidence semantics

Reference, alternate, ambiguous, and insufficient-coverage evidence are distinct states. Absence of alternate evidence is NOT automatically wild-type evidence.

Anchor/base-quality/coverage requirements are part of variant-call semantics and must remain explicit.

## Coordinates and alleles

Reference build, chromosome naming, coordinate convention, REF/ALT orientation, and strand handling must remain consistent between index construction and read matching.
