# Lumrik

**High-performance single-cell sequencing tools in Rust**

Lumrik is a Rust workspace for processing, mapping, quantifying, and analysing single-cell sequencing data. It is built around a simple principle: large sequencing datasets should be processed as streams whenever possible, with compact representations and clearly separated components rather than repeatedly materialising large intermediate datasets.

The workspace contains reusable crates for primer and barcode detection, read normalisation, external mapper integration, gene and SNP quantification, sparse single-cell data storage, feature mapping, and guide assignment.

At the centre of the workspace is **nelrune**, Lumrik's end-to-end single-cell processing pipeline.

> **Status:** Lumrik is under active development. Interfaces and output formats may still change.

---

## nelrune

`nelrune` connects the Lumrik components into a complete sequencing workflow.

A typical Illumina run looks approximately like:

```
FASTQ R1/R2
     │
     ▼
sc_primer / bam_tide
├── chemistry / primer detection
├── cell barcode extraction
├── UMI extraction
└── molecule deduplication
     │
     ├──────────────► fast feature matching
     │                ├── sample tags
     │                ├── HTO
     │                └── CRISPR guides
     │
     ▼
streaming mapper
├── STAR
├── minimap2
└── BWA
     │
     ▼
BAM stream
     │
     ▼
bam_tide
├── exonic counts
├── intronic counts
└── optional SNP-aware quantification
     │
     ▼
sparse single-cell output
```

ONT/Dorado BAM input is also supported through the normalisation layer.

The important part is that these stages are designed to cooperate as a pipeline. Reads do not need to be converted into a succession of enormous intermediate representations before the next component can begin working.

---

## Performance

Lumrik is designed for datasets containing tens of millions of sequencing reads while keeping memory consumption bounded wherever possible.

During development testing on **28 August 2026**, nelrune processed a real BD Rhapsody dataset using STAR with four mapper threads on a workstation with only **31 GiB of RAM**.

Observed during the run:

* more than **40 million reads** processed without memory growth proportional to read count;
* approximately **14,000–18,000 reads/second** end-to-end during the observed portions of the run;
* STAR itself occupied approximately **23 GiB resident memory** with the genome index loaded;
* nelrune continued normalisation, feature classification, mapping-result processing, and molecule accounting under the remaining memory constraint;
* the same development workload had previously exhausted the machine at approximately 8 million reads before removal of unnecessary per-read state retention.

These numbers are **development observations, not a controlled benchmark**. Performance depends strongly on sequencing chemistry, reference genome, mapper, storage, compression, CPU, feature references, and enabled analysis stages.

They do, however, demonstrate an important design property: processing additional reads does not inherently require retaining all previously processed read metadata in memory.

For production whole-genome STAR workloads, substantially more than 32 GiB RAM is recommended. The test above intentionally operated very close to the hardware limit.

---

## Workspace

Lumrik is composed of small crates that can also be used independently.

### `nelrune`

End-to-end orchestration for single-cell sequencing analysis.

It connects input normalisation, feature detection, external mapping, BAM processing, sparse quantification, reporting, and the live health server.

### `sc_primer`

Chemistry-aware primer, barcode, UMI, adapter, and insert detection.

The grammar supports sequencing-system-specific read structures rather than hard-coding one platform into the rest of the pipeline.

It is used for both Illumina and long-read preprocessing.

### `bam_tide`

Read normalisation and BAM quantification.

Responsibilities include:

* Illumina FASTQ processing;
* ONT/Dorado BAM processing;
* cell/UMI handling;
* molecule deduplication;
* splice-aware gene quantification;
* optional SNP-aware processing;
* sparse quantification output.

### `sc-mapper`

Streaming interface to external sequence aligners.

Currently designed around:

* STAR;
* minimap2;
* BWA.

Mapper input and BAM output are streamed so that mapping can proceed while nelrune continues processing sequencing reads.

### `gtf_splice_index`

Compact genomic annotation/index support used for splice-aware gene assignment.

### `snp-index`

Reference-genome and VCF-backed SNP support.

Provides genomic reference access and aligned-read refinement for allele-aware quantification without duplicating the genome for individual read jobs.

### `fast_tag_mapper`

Fast matching of reads against small feature reference sets such as:

* sample tags;
* hashtag oligonucleotides;
* CRISPR guide features.

### `scdata`

Sparse single-cell data structures used by Lumrik analysis components.

### `sc-beacon`

Ambient-aware CRISPR guide assignment from single-cell guide-count data.

The caller supports single- and multi-guide assignments and models background guide signal rather than relying only on a winner-takes-all threshold.

### Supporting crates

The workspace also contains focused utility crates including:

* `read-tag-table`
* `mapping_info`
* `int_to_str`
* `onehot_dna`

---

## Building

Lumrik currently targets modern Rust.

Clone the repository and build the complete workspace:

```
cargo build --release --workspace
```

To build only nelrune:

```
cargo build --release --bin nelrune
```

The resulting binary is located under:

```
target/release/nelrune
```

External mappers are not bundled. Install the mapper required for your workflow separately.

---

## Running nelrune

The exact command depends on sequencing chemistry and mapper configuration.

A simplified Illumina example is:

```
nelrune \
    --r1 sample_R1.fastq.gz \
    --r2 sample_R2.fastq.gz \
    --chemistry <CHEMISTRY> \
    --mapper star \
    --mapper-index /path/to/star_index \
    --mapper-threads 8 \
    --index reference.splice.idx \
    --threads 8 \
    --outpath nelrune_out
```

A simplified ONT example is:

```
nelrune \
    --bam dorado.bam \
    --chemistry <CHEMISTRY> \
    --mapper minimap2 \
    --mapper-index reference.mmi \
    --mapper-threads 8 \
    --index reference.splice.idx \
    --threads 8 \
    --outpath nelrune_out
```

Use:

```
nelrune --help
```

for the options supported by the current build.

---

## Live progress

nelrune includes a lightweight health/progress server for monitoring long-running analyses.

By default it uses port:

```
8787
```

Available endpoints include:

```
/          live dashboard
/health    health probe
/status    machine-readable status
```

This is particularly useful for long analyses running on workstations or compute nodes where the terminal itself does not need to remain the primary progress display.

---

## Design principles

### Stream first

Sequencing datasets are large enough already. Components should not create another complete copy of the dataset merely to pass information to the next stage.

### Keep ownership clear

Long-lived data belongs to the component responsible for it. Per-read jobs should contain read-specific state, not copies of global reference structures.

### Bound transient memory

Buffers should be large enough for efficient processing, but not so large that batching itself becomes the dominant memory consumer.

### Preserve information

Multimapping, splice structure, molecule identity, feature assignments, and allele observations should not be prematurely collapsed when later analysis may require them.

### Separate mechanisms

Primer detection, mapping, quantification, sparse storage, and statistical calling are separate crates. nelrune orchestrates them rather than reimplementing them.

### Measure on real data

Synthetic tests are essential for correctness, but performance decisions are also tested against real single-cell sequencing datasets.

---

## Development status

Lumrik is currently research software under active development.

The project already contains substantial working implementations for:

* Illumina single-cell normalisation;
* ONT/Dorado input;
* BD Rhapsody read structures;
* external streaming mapping;
* STAR, minimap2, and BWA integration;
* feature/sample-tag matching;
* exon/intron quantification;
* SNP-aware quantification;
* sparse single-cell data;
* CRISPR guide assignment;
* live progress reporting.

Additional sequencing systems, analysis modes, validation, documentation, and performance work are ongoing.

Do not assume that command-line interfaces or file formats are stable between development versions.

---

## Testing

Run the workspace tests with:

```
cargo test --workspace
```

For performance-sensitive code, release builds should always be used:

```
cargo test --release --workspace
```

Individual crates can be tested separately, for example:

```
cargo test -p sc_primer
cargo test -p sc-mapper
cargo test -p bam_tide
cargo test -p sc-beacon
```

Some integration tests require external tools or larger reference/test datasets and may therefore be ignored during a normal test run.

---

## Why Rust?

Single-cell sequencing pipelines combine several workloads that benefit from Rust:

* parsing very large compressed sequencing files;
* compact barcode and UMI representations;
* parallel sequence processing;
* sparse numerical data;
* high-throughput hashing;
* safe concurrency;
* streaming between processes;
* predictable ownership of large data structures.

Lumrik uses Rust not merely as a wrapper around existing bioinformatics tools, but as the implementation language for the data-intensive parts surrounding them.

External aligners such as STAR remain extremely specialised and highly optimised tools. Lumrik integrates them rather than attempting to replace them.

---

## License

Lumrik is developed by **Stefan Lang**.

The workspace is available under the **GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later)**.

A commercial licensing option is also available. See `Commercial.md` for details.

Third-party components retain their respective licenses. See `THIRD_PARTY_NOTICES.md`.

---

## Project philosophy

Lumrik is intended to make sophisticated single-cell sequencing analysis possible without turning every new assay into another monolithic pipeline.

The long-term goal is a collection of efficient, interoperable building blocks where sequencing chemistry, mapping strategy, feature detection, quantification, and downstream analysis can evolve independently — while nelrune provides a practical route through them for complete datasets.

