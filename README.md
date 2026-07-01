# scdata

[![Rust](https://github.com/stela2502/scdata/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/scdata/actions/workflows/rust.yml)

A high-performance Rust library for constructing **sparse single-cell UMI count data**, with deterministic export to Matrix Market (10x-style) format.


`scdata` is built around a strict separation of concerns:

- data accumulation (UMI-aware)
- export (deterministic, index-driven)

It intentionally does **not** implement matrix operations.

---

## Overview

`scdata` is designed for building **production-grade single-cell pipelines** in Rust.

It provides:

- Incremental UMI-aware insertion
- Duplicate detection
- Efficient merging
- Deterministic export using a **canonical feature index**
- Matrix Market compatibility

---

## Core Concepts

### Scdata

Sparse container for single-cell data.

- Stores cells in 256 buckets (hash-partitioned)
- Tracks unique `(feature_id, UMI)` pairs
- Aggregates counts per feature

---

### CellData

Per-cell storage:

- `seen`: unique `(feature_id, UMI)` pairs
- `total_reads`: aggregated counts per feature


### Cell Identifiers

`scdata` stores cell identifiers internally as `u64` values.

If you want exported `barcodes.tsv.gz` entries to appear as DNA barcode strings, these `u64` values must encode the barcode sequence in the 2-bit format used by `int_to_str`.

During accumulation, `scdata` treats cell IDs simply as numeric identifiers. During export, those numeric values are rendered back to DNA strings using `int_to_str`.

#### Encoding DNA barcodes for use with `scdata`

```toml
[dependencies]
stela-int-to-str = "0.1"
```

```rust
use int_to_str::int_to_str::IntToStr;

let barcode = "ATGACTCTCAGCATGG";
let cell_id: u64 = IntToStr::new(barcode.as_bytes()).into_u64();
```

You can then use `cell_id` as the cell identifier in `Scdata`:

```rust
scdata.try_insert(&cell_id, gene_umi_hash, 0.0, &mut report);
```

#### Important note

`scdata` does not validate whether a given `u64` really represents a valid DNA barcode.

If you pass arbitrary numeric IDs, they will still be stored correctly, but barcode export will interpret them as 2-bit encoded DNA and render the corresponding sequence.

---

### FeatureIndex (critical)

Defines the **canonical feature space** and export order.

```rust
pub trait FeatureIndex {
    fn feature_name(&self, feature_id: u64) -> &str;
    fn feature_id(&self, name: &str) -> Option<u64>;
    fn ordered_feature_ids(&self) -> Vec<u64>;
    fn to_10x_feature_line(&self, feature_id: u64) -> String;
}
```

👉 This ensures stable row order across exports.

---

### MatrixValueType

Controls Matrix Market output:

- Integer (recommended)
- Real (limited support)

---

## Workflow

### 1. Create dataset

```rust
let mut data = Scdata::new(4, MatrixValueType::Integer);
```

### 2. Insert data

```rust
data.try_insert(&cell_id, gene_umi_hash, 0.0, &mut report);
```

### 3. Finalize for export

```rust
data.finalize_for_export(min_total_umis, &feature_index);
```

### 4. Export

```rust
data.write_sparse(&output_dir, &feature_index)?;
```

---

## Matrix Market Output

Produces standard 10x-style files:

- matrix.mtx.gz
- features.tsv.gz
- barcodes.tsv.gz

Feature order is derived from `FeatureIndex`, not data.

---

## Design Philosophy

- UMI-centric (not raw count matrix first)
- Deterministic export
- Explicit feature indexing
- Separation of:
  - data ingestion
  - finalization
  - export
- No hidden magic

---

## When to use scdata

✔ Rust-native single-cell pipelines  
✔ High-performance UMI counting  
✔ Deterministic reproducible outputs  
✔ Avoid Python/R dependencies

---

## Limitations

- Real-valued matrices are not a primary target
- Requires external feature index for export
- MatrixMarket import reconstructs integer counts only

---

## License

MIT or Apache-2.0

---

## Author

Stefan Lang