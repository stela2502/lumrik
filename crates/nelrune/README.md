# nelrune

> **Lumrik crate.** Nelrune is part of [Lumrik](../../README.md) and is distributed as part of the Lumrik workspace. See the [Lumrik README](../../README.md#license) and root `LICENSE` for the workspace licensing terms.

## What this crate does

Nelrune is Lumrik's integrated single-cell processing pipeline. It connects primer/read-structure detection, cell and UMI extraction, feature routing, molecule deduplication, external mapping and BAM quantification without turning every stage into another full intermediate dataset.

## Why Nelrune exists

The main goal is not to replace excellent specialized aligners. It is to make sure they only do work that is still biologically useful.

A sequencer can produce tens or hundreds of millions of reads, but many reads can be rejected or classified before genome alignment. A PCR duplicate does not become more informative because STAR maps it again. A sample-tag or guide read does not become useful genomic evidence because STAR is allowed to fail on it. Primer/barcode failures do not need to enter the mapper at all.

Nelrune therefore performs the cheap, chemistry-aware decisions first and sends only accepted genomic molecules to the external mapper. This reduces wasted CPU time, mapper I/O and temporary data, while preserving the specialized mapping quality of STAR, minimap2 or BWA.

This is also a sustainability decision: **the cheapest computation is the computation we can prove we do not need to perform.** On large sequencing runs, avoiding unnecessary mapping is more useful than merely making unnecessary mapping slightly faster.

The intended flow is:

```text
FASTQ / normalized reads
        |
        v
primer + cell/UMI detection
        |
        +---- feature/sample-tag/guide reads ---> feature processing
        |
        +---- duplicates / unusable reads ------> stop
        |
        v
accepted genomic molecules
        |
        v
STAR / minimap2 / BWA
        |
        v
BAM processing + quantification
```

## Binaries

- `nelrune` — the production integration binary.

V(D)J reconstruction is implemented in the separate [`sc-vdj`](../sc-vdj/README.md) crate and exposed through `nelrune-vdj`. Keeping it downstream means the primary mapper BAM can be reused as receptor evidence rather than requiring Nelrune's core FASTQ/mapping orchestration to contain receptor-specific logic.

## Library use

Nelrune is primarily an integration binary. Its orchestration code is also the reference implementation for composing [`sc_primer`](../sc_primer/README.md), [`bam_tide`](../bam_tide/README.md), [`sc-mapper`](../sc-mapper/README.md), [`scdata`](../scdata/README.md), `mapping_info` and `lumrik-status`.

## Detailed documentation

## Build

From the Lumrik workspace root:

```bash
cargo build --release --bin nelrune
```

## Illumina example

```bash
cargo run --release --bin nelrune -- \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --chemistry <CHEMISTRY> \
  --mapper minimap2 \
  --mapper-index reference.mmi \
  --mapper-threads 8 \
  --index reference.splice.idx \
  --threads 8 \
  --outpath nelrune_out
```

## ONT example

```bash
cargo run --release --bin nelrune -- \
  --bam dorado.bam \
  --chemistry <CHEMISTRY> \
  --mapper minimap2 \
  --mapper-index reference.mmi \
  --mapper-threads 8 \
  --index reference.splice.idx \
  --threads 8 \
  --outpath nelrune_out
```

Use the exact `sc_primer` chemistry/options appropriate for the experiment.

## Live health server

The server binds to `0.0.0.0` so it is reachable through the host/node network.
Nelrune prints the externally useful URL using, in order:

1. `--health-hostname`
2. `SLURMD_NODENAME`
3. `hostname -f`
4. `HOSTNAME`
5. `localhost`

Default port: `8787`.

Endpoints:

- `/` live dashboard
- `/health` simple `OK` probe
- `/status` JSON state

On clusters where the compute-node hostname used by your browser differs from
`SLURMD_NODENAME`, pass it explicitly with `--health-hostname`.

## Output

- `nelrune.log` — stage/progress log plus the original normalizer and quantifier `MappingInfo` reports
- `nelrune-report.txt` — final quantification report, including permanent cell-calling/accounting diagnostics
- `exonic/`, `intronic/`, and optional SNP output directories — **filtered canonical-cell matrices**
- `raw/exonic/`, `raw/intronic/`, and optional SNP directories — **all observed barcode evidence (>=1 UMI) before canonical cell calling**
- `filtered/exonic/`, `filtered/intronic/`, and optional SNP directories — **canonical called cells only**
- mapper BAM only when `--bam-out` is supplied; otherwise the temporary mapper BAM is removed after quantification

### Filtered and unfiltered matrices are both part of the output contract

Nelrune must preserve the pre-cell-calling GEX evidence. The normal top-level
matrices contain the canonical called cells used downstream; the matching
`unfiltered/` matrices contain every observed barcode with at least one UMI.
The unfiltered matrices are intentionally retained so that barcode-rank/cell
calling can be audited, alternative callers can be tested, and evidence is not
lost when a caller is unexpectedly stringent. **Do not remove the unfiltered
export as an output-size optimization or replace it with filtered-only output.**

`nelrune-report.txt` records the number of exonic/intronic barcodes before the
cutoff, UMI totals, counts above standard UMI thresholds, the cell-calling
method/cutoff, retained cells, and cells not called.

## Deliberate minimalism

This version does not add a second summary hierarchy. `RunProgress` owns only
live orchestration state; `MappingInfo` remains the source of truth for
normalization/quantification counters and report strings.
