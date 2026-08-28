# nelrune

Minimal Lumrik orchestration binary for a complete single-cell analysis pass:

1. normalize ONT/Dorado BAM **or** Illumina R1/R2 input with `bam_tide`
2. stream normalized molecules into `sc-mapper`
3. drain mapper results into a temporary BAM
4. quantify that BAM with `bam_tide::BamCollector`
5. write sparse quantification output
6. write the authoritative `MappingInfo` report
7. expose live orchestration state through the health server

Nelrune deliberately does not duplicate normalizer or quantifier counters. The
normalizer and quantifier reports are printed and written using their existing
`MappingInfo` display strings.

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
- `nelrune-report.txt` — final quantification `MappingInfo`
- `exonic/`, `intronic/`, and optional SNP output directories from `QuantData::write`
- mapper BAM only when `--bam-out` is supplied; otherwise the temporary mapper BAM is removed after quantification

## Deliberate minimalism

This version does not add a second summary hierarchy. `RunProgress` owns only
live orchestration state; `MappingInfo` remains the source of truth for
normalization/quantification counters and report strings.
