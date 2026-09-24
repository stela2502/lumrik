# sc-vdj AI Contract

## Role

`sc-vdj` owns single-cell V(D)J evidence, receptor reconstruction, recombination inference, rescue, and AIRR-compatible output. VDJ-origin (`V`) molecules belong here rather than in GEX quantification.

## Evidence first

Preserve per-cell receptor evidence and distinguish direct evidence, rescued/fragment evidence, germline/reference support, and inferred recombination state. Do not turn an inference into apparent direct evidence.

## Recombination and clone identity

V/D/J/C calls and junction/CDR3 evidence have different confidence and biological roles. Do not merge receptor families solely because segment labels look similar.

When reconstructing receptors, support multiple recombinations per locus where evidence supports them; do not force one receptor per locus by construction.

## AIRR output

AIRR fields must reflect the reconstructed biological state consistently, including locus, V/D/J/C calls, productivity, frame/stop status, junction metrics, and sequence fields. Missing/uncertain evidence should remain missing/uncertain rather than being fabricated to satisfy a schema.

## Performance

The crate is expected to operate on large BAMs. Chunking/batching may change implementation strategy but MUST NOT change per-cell evidence semantics or make results depend on chunk boundaries.
