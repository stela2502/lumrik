# gtf_splice_index AI Contract

## Role

This crate is the annotation/model layer for splice-aware genomic matching. It is BAM-independent and operates on genomic blocks using 0-based, half-open coordinates.

## Match semantics

A match class must describe what was positively established about the relationship between a read and a transcript model. Do not use a class merely as a convenient fallback.

- `ExactJunctionChain`: read junction chain exactly matches the transcript junction chain.
- `Compatible`: read blocks/junctions are consistent with the transcript model.
- `JunctionMismatch`: exonic structure is present but one or more observed splice junctions are absent from that transcript model.
- `Intronic`: positive overlap with an annotated intronic interval.
- `Incompatible`: the transcript model cannot explain the alignment; this is distinct from positive intronic evidence.

Failure of exon fitting MUST NOT by itself imply `Intronic`.

## Transcript boundaries

Transcript 5'/3' boundaries and poly(A)-proximal annotation are not exact experimental boundaries. Terminal exon matching must support biological end overhangs and must not reject ordinary reads merely because annotation ends are slightly short. Lumrik currently requires an effective minimum terminal-overhang tolerance of 100 bp.

Boundary handling must work for both strands and for single-exon transcripts. Do not derive biological 5'/3' semantics from block order alone without considering transcript orientation/model boundaries.

## Separation of concerns

This crate classifies transcript/read relationships. It does NOT decide whether a class is counted as GEX, intronic RNA, VDJ, or discarded by a downstream assay. Downstream crates must make those policy decisions explicitly.
