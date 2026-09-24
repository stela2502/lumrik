# clonomap AI Contract

## Role

ClonoMap infers mutation-aware substructure within large receptor clones using scalable geometric methods (PCA/MST), not full phylogenetic reconstruction.

## Interpretation

Geometric proximity/MST edges are analysis structure, not automatically evolutionary ancestry. Do not label them as phylogenetic relationships without an explicit model supporting that interpretation.

## Mutation-aware input

Reference/germline context and mutation representation are part of the model. Do not mix incompatible receptor/reference coordinate systems silently.

## Scale

Algorithms should remain viable for very large clones; avoid accidental all-pairs/quadratic replacements without an explicit reason.
