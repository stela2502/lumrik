# sc_primer

A grammar-driven Rust crate for detecting, validating, correcting and generating
single-cell primer structures for 10x Genomics, BD Rhapsody, Illumina and ONT workflows.

## Features

- 10x Genomics chemistry presets
- BD Rhapsody chemistry presets
- Custom grammar definitions
- Exact whitelist validation
- Unique 1-mismatch fuzzy rescue
- ONT multimer detection
- Reverse-complement detection
- Primer synthesis and regeneration
- Cell barcode translation tables
- Reusable clap CLI integration

## Installation and availability

`sc_primer` is currently intended to be installed directly from GitHub.

The crate embeds large 10x whitelist resources, especially the 3M whitelist used
for modern 10x chemistries. These data files make the package too large for a
comfortable crates.io release, so the supported dependency style is currently a
Git dependency.

```toml
[dependencies]
sc_primer = { git = "https://github.com/stela2502/bam_tide", package = "sc_primer" }
```

For local development:

```toml
[dependencies]
sc_primer = { path = "../sc_primer" }
```

The Rust package registry itself is called **crates.io**, but this crate is not
currently distributed there because of the embedded whitelist size.

## Core Concepts

```text
Chemistry
    ↓
Grammar
    ↓
PrimerDetector
    ↓
PrimerMatch
```

Built-in chemistries are convenience presets built on top of the grammar engine.

## Supported Chemistries

### 10x Genomics

- tenx-three-prime-v1
- tenx-three-prime-v2
- tenx-three-prime-v3
- tenx-three-prime-v4
- tenx-five-prime
- tenx-multiome-arc-v1

### BD Rhapsody

- bd-v1
- bd-v2-96
- bd-v2-384

## Quick Start

```rust
use sc_primer::{Chemistry, PrimerDetector};

let detector =
    PrimerDetector::from_chemistry(
        Chemistry::TenxThreePrimeV3
    )?;

let hits = detector.detect_all(seq, qual)?;
```

## Detection APIs

### detect_first()

Use for Illumina-style reads where only one primer is expected.

```rust
let hit = detector.detect_first(seq, qual)?;
```

### detect_all()

Use for ONT reads that may contain multiple concatenated molecules.

```rust
let hits = detector.detect_all(seq, qual)?;
```

### generate()

Generate a target primer from a ReadTagRecord.

```rust
let (target_cell, primer) =
    detector.generate(&record)?;
```

## Grammar Language

Supported operations:

```text
FIXED:<SEQ>:mm=<N>

CELL:<LEN>
TENX_CELL:<VERSION>
BD_CELL:<VERSION>

UMI:<LEN>

POLYT:min=<N>

INSERT

SEARCH:<START>..<END>

SKIP:<LEN>

TAG:<LEN>

FEATURE:<LEN>
```

### Example 10x grammar

```text
FIXED:CTACACGACGCTCTTCCGATCT:mm=2
+TENX_CELL:3p-v3
+UMI:12
+POLYT:min=0
+INSERT
```

### Example BD grammar

```text
FIXED:ATAGGAAACTCATGGT:mm=2
+BD_CELL:v2.384
+POLYT:min=0
+INSERT
```

## PrimerMatch

Extract data from a detected primer:

```rust
let cell = hit.get_cell(seq, qual)?;
let umi = hit.get_umi(seq, qual)?;
let insert = hit.get_insert(seq, qual)?;
```

## Barcode Validation

### 10x

Workflow:

```text
exact whitelist lookup
    ↓
unique 1-mismatch rescue
    ↓
fail
```

### BD Rhapsody

Workflow:

```text
C1/C2/C3 lookup
    ↓
unique 1-mismatch rescue
    ↓
cell id reconstruction
```

## ONT Multimers

The detector supports reads containing multiple:

```text
primer + insert
primer + insert
primer + insert
```

structures within the same sequencing read.

## Reverse Complement Detection

Detected matches report:

```rust
Orientation::Forward
Orientation::ReverseComplement
```

## Cell Translation

The detector can translate invalid or foreign barcodes into valid target-system
barcodes and stores the mapping internally.

```rust
let translation =
    detector.primer_translation();
```

## CLI Usage

Example:

```bash
identify_primers   --chemistry tenx-three-prime-v3   --seq ACTG...
```

Custom grammar:

```bash
identify_primers   --primer-structure 'FIXED:AAA:mm=1+CELL:16+UMI:12+INSERT'   --seq ACTG...
```

## Benchmarks

Recent release-mode measurements:

| Benchmark | µs/hit |
|------------|---------:|
| 10x exact | ~0.835 |
| 10x fuzzy | ~0.850 |
| BD exact | ~0.995 |
| BD fuzzy | ~1.767 |

The BD fuzzy benchmark intentionally exercises fuzzy rescue on every barcode and
represents a worst-case workload.

## Error Handling

```rust
PrimerResult<T>
PrimerError
```

## Development

Run tests:

```bash
cargo test
```

Run benchmarks:

```bash
cargo test --release benchmark_detect_all -- --ignored --nocapture
```

Build documentation:

```bash
cargo doc --no-deps
```

## Design Goals

- Grammar-driven architecture
- Shared implementation across chemistries
- Fast exact lookup paths
- Conservative fuzzy rescue
- ONT multimer support as a first-class feature
- Minimal overhead abstractions
