[![Rust](https://github.com/stela2502/bam_tide/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/bam_tide/actions/workflows/rust.yml)

# bam_tide

`bam_tide` is a collection of high-performance genomics and single-cell analysis tools written in Rust.

The project focuses on:

- ONT long-read processing
- single-cell quantification
- transcriptome-first workflows
- splice-aware analysis
- SNP-aware quantification
- and scalable HPC deployment

The toolkit was created to reduce dependence on large fragmented bioinformatics software stacks involving:

- Conda environments
- Python dependency chains
- R installations
- workflow glue code
- and repeated parsing of large annotation files

Instead, `bam_tide` aims to provide:

- standalone compiled binaries
- reproducible workflows
- streaming-friendly processing
- and efficient large-scale analysis

---

# Main Features

- ONT barcode and UMI extraction
- transcriptome-to-genome projection
- splice-aware quantification
- SNP side-channel quantification
- BAM coverage generation
- BigWig comparison and validation
- barcode-aware BAM subsetting
- reusable serialized splice indexes

---

# Included Tools

| Tool | Description |
|---|---|
| `bam-ont-normalizer` | Extract and normalize ONT single-cell molecules |
| `bam-transcriptome-to-genome` | Convert transcriptome BAMs back to genomic coordinates |
| `bam-quant` | Splice-aware single-cell quantification |
| `gtf-splice-index` | Build reusable serialized splice indexes |
| `bam-coverage` | Generate genomic coverage tracks |
| `bw-compare` | Compare BigWig coverage tracks |
| `bam-subset-tag` | Split BAMs using barcode/tag values |
| `bam-read-tag-stats` | QC and summarize read-tag tables |
| `bam-quant-testdata` | Generate quantification test datasets |

Detailed documentation for the software suite is available at [bam_tide documentation](https://stela2502.github.io/bam_tide/)

---

# Installation

## System Dependencies

`bam_tide` depends on HTSlib through `rust-htslib`.

Ubuntu/Debian dependencies:

```bash
sudo apt-get update

sudo apt-get install -y \
    libhts-dev \
    libbz2-dev \
    liblzma-dev \
    libz-dev \
    pkg-config
```

---

# Building

Standard development build:

```bash
cargo build --release
```

---

## Static MUSL Build

Recommended for HPC deployment and workflow portability:

```bash
cargo build --release \
    --target x86_64-unknown-linux-musl
```

This produces portable static binaries suitable for:

- HPC systems
- workflow deployment
- Singularity/Apptainer containers
- and reproducible pipelines

---

# Deploying Binaries Into Pipelines

Example deployment workflow:

```bash
find target/x86_64-unknown-linux-musl/release \
  -maxdepth 1 \
  -type f \
  -executable \
  -exec cp {} pipeline/bin/ \;
```

This copies all compiled executables into a workflow `bin/` directory.

---

# Documentation

The project documentation is written using `mdBook`.

## Install mdBook

```bash
cargo install mdbook
```

---

## Build Documentation

```bash
mdbook build book
```

Generated HTML documentation:

```text
book/book/index.html
```

---

## Live Documentation Server

```bash
mdbook serve book
```

Default local URL:

```text
http://localhost:3000
```

---

# Example ONT Workflow

```text
Dorado BAM
  ↓
bam-ont-normalizer
  ↓
normalized FASTQ
  ↓
transcriptome alignment
  ↓
bam-transcriptome-to-genome
  ↓
genome-coordinate BAM
  ↓
bam-quant
  ↓
sparse matrices + SNP matrices
```

---

# Nextflow Pipeline

The tools are designed to integrate directly into Nextflow workflows.

Example companion pipeline:

https://github.com/stela2502/sc-ont-mut-quant/

The pipeline implements:

- ONT preprocessing
- transcriptome alignment
- transcriptome-to-genome projection
- splice-aware quantification
- SNP-aware matrix generation

using the `bam_tide` toolchain.

---

# Testing

Run tests:

```bash
cargo test
```

Verbose testing:

```bash
cargo test --verbose
```

---

# GitHub CI

The repository uses GitHub Actions for:

- building
- dependency validation
- and automated testing

Current CI workflow installs HTSlib dependencies and runs:

```bash
cargo build --verbose
cargo test --verbose
```

---

# Design Philosophy

The project emphasizes:

- performance
- reproducibility
- deployment simplicity
- low runtime overhead
- and scalable long-read analysis

The goal is not to replace all bioinformatics tooling, but to accelerate and simplify the computationally expensive parts of sequencing workflows.

---

# Current Status

`bam_tide` is an actively evolving research-oriented toolkit.

Interfaces and output formats may still evolve, especially in experimental ONT and SNP-processing components.

The current focus is:

- correctness
- reproducibility
- workflow integration
- and scalable performance
