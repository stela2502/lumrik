# nelrune AI Contract

## Role

`nelrune` orchestrates single-cell preparation, mapping, GEX quantification, feature-tag integration, and reporting. It wires specialized crates together; it must not blur their biological boundaries.

## Modality ownership

GEX quantification/canonical GEX cell calling uses GEX-origin molecules. VDJ reconstruction belongs to `sc-vdj`. Nelrune must preserve primer/grammar provenance when routing data between stages.

Do not make canonical GEX cell calling depend on VDJ abundance, feature-tag abundance, or another modality's molecule distribution.

## Canonical cells

Once canonical GEX cells are called, publish that count immediately to progress/server state and use the same authoritative set for downstream filtered outputs. Do not maintain competing hidden definitions of "retained cell".

An explicit user-provided `--min-umi-count` is authoritative cell calling. It MUST bypass sc-beacon entirely; automatic knee calling is used only when no explicit cutoff is supplied.

## Health server

The health server is a correctness/debugging instrument, not cosmetic UI.

Only display metrics meaningful for the active stage/mode. In BAM-only quantification, FASTQ throughput/routing counters that were never populated should not be presented as meaningful zeros.

Expose enough live biological/algorithmic counters to diagnose distributional failures before completion, including expression totals, canonical cell count when available, and transcript match-class populations.

A final-stage transition must not erase already-known useful measurements.

## Orchestration

Nelrune should call crate APIs according to their contracts rather than reimplementing primer detection, splice semantics, sparse-matrix semantics, or VDJ reconstruction locally.
