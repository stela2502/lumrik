# Ommverse

Ommverse is Lumrik's integrated genome-to-protein reference model. It builds a compact, serialised reference for one genome assembly from UCSC-hosted genome, gene and UniProt annotation, then adds reference regulatory annotation such as FANTOM5 enhancers and the compact union of ENCODE4 protein-binding regions when those sources are supported for the assembly.

Large experimental evidence sets are deliberately kept in the source cache rather than copied wholesale into the `.ommverse` file. For example, the full ENCODE4 TF rPeak BigBed remains available under `Chromatin/ENCODE4/` while Ommverse stores only its merged candidate binding-region union.

The normal workflow has two phases:

1. build the base genome/reference models;
2. enrich one or many base models with InterPro protein annotations in a single streaming pass.

## Build

From the Lumrik repository root:

```bash
cargo build -r -p ommverse
```

The examples below assume the release binary:

```bash
OMMVERSE=target/release/ommverse
```

## 1. Build a base model

`ommverse build --assembly` downloads the sources known to Ommverse, caches them, builds the model and writes a build summary beside it.

For hg38:

```bash
$OMMVERSE build \
    --assembly hg38 \
    --cache /data1/UCSC \
    --out /data1/UCSC/hg38/hg38.ommverse
```

The cache layout is assembly-centred:

```text
/data1/UCSC/hg38/
├── Genome/
├── Genes/
├── Protein/
├── Chromatin/
│   ├── FANTOM5/
│   └── ENCODE4/
├── hg38.ommverse
└── hg38.ommverse-build-summary.yaml
```

Existing files are reused. Re-running the same command therefore rebuilds from the cached source set and only fetches required files that are missing.

The build health server is enabled by default on port 8787. Its URL is printed when the build starts. The standard controls are:

```text
--health-port <PORT>
--health-hostname <HOSTNAME>
--no-health-server
```

To build the reference assemblies currently used by Lumrik, run the same operation for each assembly:

```bash
for assembly in hg38 mm39 ce11 danRer11 dm6; do
    $OMMVERSE build \
        --assembly "$assembly" \
        --cache /data1/UCSC \
        --out "/data1/UCSC/$assembly/$assembly.ommverse" || exit 1
done
```

Source availability is assembly-dependent. Optional annotation sources are skipped when Ommverse does not support them for that assembly; the genome/gene sources required to construct the base model must be available.

### Build from an already downloaded reference tree

If the UCSC-style source tree already exists and no fetching is wanted:

```bash
$OMMVERSE build \
    --reference /data1/UCSC/hg38 \
    --out /data1/UCSC/hg38/hg38.ommverse
```

`--reference` and `--assembly` are mutually exclusive.

## What is already in the base model?

The base build includes the genome/transcript model and the UniProt information imported from the UCSC UniProt tracks. The build summary reports, among other things, mapped proteins, linked transcripts, protein feature records, chromatin elements and binding-union regions.

InterPro is a separate enrichment step. A base file named `hg38.ommverse` is therefore intentionally distinct from the InterPro-enriched `hg38.interpro.ommverse` described below.

## 2. Download the InterPro bulk annotation once

The bulk importer consumes the InterPro files directly; it does not require running InterProScan for every Ommverse protein. InterPro distributes these data from its current-release download area.

Choose a shared cache location, for example:

```bash
mkdir -p /data1/UCSC/InterPro
cd /data1/UCSC/InterPro

wget -c https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list
wget -c https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt
```

`protein2ipr.dat.gz` is streamed and can be large. `entry.list` supplies the InterPro entry vocabulary. `ParentChildTreeFile.txt` is optional to the CLI, but supplying it preserves the InterPro parent/child relationships and is recommended for the full model.

## 3. Add the protein dimension to one model

For a single base index, `--out` is required:

```bash
$OMMVERSE ingest-interpro \
    --index /data1/UCSC/hg38/hg38.ommverse \
    --protein2ipr /data1/UCSC/InterPro/protein2ipr.dat.gz \
    --entry-list /data1/UCSC/InterPro/entry.list \
    --parent-child-tree /data1/UCSC/InterPro/ParentChildTreeFile.txt \
    --out /data1/UCSC/hg38/hg38.interpro.ommverse
```

The importer matches InterPro records by UniProt accession to proteins already present in Ommverse. It reports matching records, features added, duplicates, malformed records and unknown entries.

The original base model is not overwritten by this command.

## 4. Bulk-enrich all models

This is the preferred route after building several assemblies. `--index-dir` recursively finds base `*.ommverse` files and deliberately ignores files already ending in `.interpro.ommverse`.

```bash
$OMMVERSE ingest-interpro \
    --index-dir /data1/UCSC \
    --protein2ipr /data1/UCSC/InterPro/protein2ipr.dat.gz \
    --entry-list /data1/UCSC/InterPro/entry.list \
    --parent-child-tree /data1/UCSC/InterPro/ParentChildTreeFile.txt
```

The large `protein2ipr.dat.gz` file is streamed **once for all discovered Ommverse indices**. For every input such as:

```text
/data1/UCSC/hg38/hg38.ommverse
```

Ommverse writes:

```text
/data1/UCSC/hg38/hg38.interpro.ommverse
```

The InterPro importer also exposes the same health-server controls as the base builder and reports progress for the shared stream and each target index.

This two-stage design is intentional: rebuilding genome/reference annotation does not require repeatedly streaming the global InterPro mapping, while one InterPro pass can enrich every assembly after their base models are ready.

## Sanity checks

Inspect a gene and its mapped proteins:

```bash
$OMMVERSE gene \
    --index /data1/UCSC/hg38/hg38.interpro.ommverse \
    TP53
```

Inspect one UniProt accession, including its feature list:

```bash
$OMMVERSE protein \
    --index /data1/UCSC/hg38/hg38.interpro.ommverse \
    P04637
```

Include the reconstructed protein sequence with:

```bash
$OMMVERSE protein \
    --index /data1/UCSC/hg38/hg38.interpro.ommverse \
    P04637 \
    --sequence
```

If the base build reports proteins but unexpectedly reports zero protein features, do not treat the InterPro step as a repair for a broken base import: the UCSC/UniProt protein import and the later InterPro enrichment are separate layers and should be diagnosed separately.

## Reproducibility and source cache

Ommverse keeps downloaded source files instead of deleting them after import. The source tree is therefore both a cache and the provenance/evidence reservoir for later Lumrik stages.

In particular, large state-dependent experimental resources do not automatically become static biological truth inside Ommverse. The cached source may contain substantially richer information than the compact `.ommverse` representation. Later evidence-generating code can query those sources directly when it needs the original observations.
