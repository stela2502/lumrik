# Introduction

`bam_tide` is a collection of high-performance genomics and single-cell analysis tools written in Rust.

The project was created out of frustration with modern bioinformatics software stacks that often require:

- multiple incompatible Python environments
- Conda dependency chains
- R installations
- fragile workflow glue
- extensive preprocessing steps
- and large memory footprints despite mediocre runtime performance

The goal of `bam_tide` is not to re-implement entire ecosystems.  
The goal is to provide a focused toolkit that performs the expensive and repetitive parts of sequencing analysis:

- faster
- with fewer dependencies
- with reproducible behavior
- and with a simple deployment model

In most cases, installation should be as simple as copying a small set of binaries into a workflow or HPC environment.

---

# Core Philosophy

The project is built around several design principles.

## 1. Performance matters

Modern sequencing datasets are enormous.

Single-cell and long-read experiments routinely contain:

- millions of reads
- large sparse matrices
- expensive BAM traversal operations
- repeated parsing of the same annotation structures

Many existing pipelines repeatedly move large intermediate datasets through:

- Python
- R
- shell glue
- serialization layers
- and temporary file conversions

This creates unnecessary overhead.

`bam_tide` attempts to keep critical operations:

- compiled
- parallel
- streaming where possible
- and memory-efficient

Rust alone does not magically make software fast, but it enables:

- predictable memory usage
- efficient threading
- low-overhead parsing
- and safer systems-level optimization

---

## 2. Monolithic deployment over fragmented environments

A major project goal is operational simplicity.

Many sequencing workflows currently depend on combinations of:

- Python packages
- Conda environments
- R packages
- system libraries
- external scripts
- and workflow-specific wrappers

These environments are often difficult to reproduce and expensive to maintain on HPC systems.

`bam_tide` instead aims for:

- standalone binaries
- minimal runtime dependencies
- reproducible command-line interfaces
- and easy integration into Nextflow or shell-based workflows

The preferred deployment model is:

```bash
cargo build -r --target x86_64-unknown-linux-musl

find target/x86_64-unknown-linux-musl/release \
  -maxdepth 1 \
  -type f \
  -executable \
  -exec cp {} pipeline/bin/ \;

```

and run.

This produces portable static binaries that can be copied directly into
workflow environments, HPC systems, or container images without requiring
Conda, Python package stacks, or large runtime environments.

---

## 3. Streaming and direct processing

Where possible, tools are designed to operate directly on:

- BAM files
- FASTQ streams
- sparse matrices
- annotation indexes
- and compressed inputs

without large temporary conversion steps.

The project strongly prefers:

- direct processing
- piped execution
- and reusable binary indexes

over repeated re-parsing of large text-based formats.

---

## 4. Practical single-cell and ONT workflows

The toolkit is especially focused on:

- ONT long-read sequencing
- single-cell quantification
- transcriptome projection
- SNP-aware quantification
- and barcode-aware processing

Several tools are designed to work together as one integrated workflow rather than isolated utilities.

The toolkit currently contains the following command-line applications.

---

## `bam-ont-normalizer`

Preprocess ONT BAM files into transcript-compatible single-cell reads.

The tool detects and extracts:

- cell barcodes
- UMIs
- adapter structures
- poly-T sequences
- read orientation

and generates cleaned FASTQ output together with read-tag tables that can later be used during quantification.

Typical use cases:

- ONT single-cell preprocessing
- barcode extraction
- cassette detection
- read normalization
- preparation for transcriptome alignment

See: [bam-ont-normalize](tools/bam-ont-normalize.md)

---

## `bam-transcriptome-to-genome`

Project transcriptome-aligned reads back into genomic coordinates.

This tool converts transcript-space alignments into genome-space BAM files while preserving:

- splice structures
- read orientation
- transcript provenance
- and optional SNP compatibility

Typical use cases:

- transcriptome-first ONT workflows
- genome-coordinate visualization
- SNP-aware downstream quantification
- genomic BAM generation from transcript alignments

See: [bam-transcriptome-to-genome](tools/bam-transcriptome-to-genome.md)

---

## `bam-quant`

Single-cell quantification engine for transcript and gene expression.

`bam-quant` operates directly on BAM files using a compiled splice index and can optionally integrate:

- external read-tag tables
- SNP-aware quantification
- transcript-level quantification
- intronic splitting
- and multi-BAM merged analysis

The tool produces sparse matrix outputs compatible with downstream single-cell workflows.

Typical use cases:

- ONT single-cell quantification
- transcript quantification
- SNP side-channel counting
- barcode-aware quantification
- large-scale matrix generation

See: [bam-quant](tools/bam-quant.md)

---

## `gtf-splice-index`

Compile GTF/GFF annotations into a reusable binary splice index.

The generated index is optimized for rapid transcript and splice matching during quantification and projection workflows.

Typical use cases:

- transcript annotation indexing
- splice-aware matching
- reusable workflow preprocessing

See: [gtf-splice-index](tools/gtf-splice-index.md)

---

## `bam-coverage`

Generate coverage tracks from BAM files.

The tool is designed as a fast alternative to common BAM coverage utilities and supports filtering behavior similar to `deeptools bamCoverage`.

Typical use cases:

- genome browser visualization
- coverage generation
- BigWig and bedGraph workflows
- QC and inspection

See: [bam-coverage](tools/bam-coverage.md)

---

## `bw-compare`

Compare and combine BigWig coverage tracks.

This utility supports efficient signal comparison workflows for downstream visualization and analysis.

Typical use cases:

- coverage normalization
- differential track comparison
- signal subtraction and scaling

See: [bw-compare](tools/bw-compare.md)

---

## `bam-subset-tag`

Subset BAM files using tag-aware filtering.

Useful for extracting:

- cells
- barcodes
- read groups
- or experimental subsets

from large sequencing datasets.

Typical use cases:

- debugging
- targeted re-analysis
- workflow prototyping
- barcode-specific extraction

See: [bam-subset-tag](tools/bam-subset-tag.md)

---

## `bam-read-tag-stats`

Summarize and inspect external read-tag tables.

Provides statistics for:

- cell barcode usage
- UMI distributions
- pair frequencies
- and barcode complexity

Typical use cases:

- QC of ONT normalization
- barcode inspection
- workflow debugging
- read-tag validation

See: [bam-read-tag-stats](tools/bam-read-tag-stats.md)

---

## `sum_up_ont_tab_result`

Utility helper for summarizing ONT intermediate table outputs.

Mainly intended for:

- debugging
- workflow inspection
- and experimental preprocessing analysis

See: [sum_up_ont_tab_result](tools/sum_up_ont_tab_result.md)

---

## `bam-quant-testdata`

Generate synthetic or controlled test datasets for quantification validation.

Useful for:

- integration testing
- reproducibility checks
- benchmarking
- and development validation

See: [bam-quant-testdata](tools/bam-quant-testdata.md)

---


# Intended Use Cases

The project is particularly suited for:

- HPC environments
- reproducible pipelines
- large sequencing experiments
- long-read transcriptomics
- rapid iterative analysis
- and environments where deployment simplicity matters

---

# Current Status

`bam_tide` is an actively evolving research-oriented toolkit.

Interfaces and output formats may still evolve, especially in experimental ONT and SNP-processing components.

The project currently prioritizes:

- correctness
- reproducibility
- performance
- and workflow integration

over strict backward compatibility during early development.
