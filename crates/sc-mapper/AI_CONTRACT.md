# sc-mapper AI Contract

## Role

`sc-mapper` is the external-mapper integration layer (for example STAR/minimap2/BWA). It streams reads to mapper processes and returns alignments while preserving Lumrik molecule identity.

## Identity preservation

Read/QNAME normalization must preserve the information required to recover cell, UMI, quality, and grammar provenance downstream. Do not replace corrected Lumrik identity with an original/raw identifier when writing BAM output.

## Mapper boundary

Mapper-specific process/FIFO details belong here; biological quantification does not. A mapper backend change must not alter cell/UMI/provenance semantics.

## Streaming

FIFO/stdin/stdout orchestration must be safe for large datasets and must not require loading the complete dataset into memory.
