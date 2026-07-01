[![Rust](https://github.com/stela2502/mapping_info/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/mapping_info/actions/workflows/rust.yml)
# mapping_info

[![Rust](https://github.com/stela2502/mapping_info/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/mapping_info/actions/workflows/rust.yml)

A lightweight Rust helper library for collecting, timing, logging and
reporting mapping statistics in sequencing pipelines.

This crate is designed for:

- BAM processing tools
- Single-cell quantification workflows
- Multi-threaded pipelines
- Performance benchmarking
- Structured error reporting

---

# Overview

`MappingInfo` is a structured statistics collector that tracks:

- Read filtering categories
- Cell/gene assignment success
- PCR duplicates
- Multi-mappers
- Error types
- Per-category read logs
- Histogram counters
- Detailed timing measurements
- File I/O vs CPU vs subprocess time

The **recommended modern usage pattern** is centered around:

```rust
mapping_info.report("some_error_type");
```

---

# Installation

```toml
mapping_info = { git = "https://github.com/stela2502/mapping_info" }
```

---

# Basic Usage

## Create a MappingInfo instance

```rust
use mapping_info::MappingInfo;
use std::fs::File;

let logfile = File::create("run.log").ok();

let mut info = MappingInfo::new(
    Some(logfile),    // Optional<File>
    20.0,       // min_quality
    0           // max_reads (0 = unlimited)
);
```

---

# Error Reporting (Recommended)

Modern usage should rely on the flexible `report()` API.

```rust
info.report("low_quality");
info.report("invalid_barcode");
info.report("invalid_barcode");
```

Retrieve counts:

```rust
let n = info.get_issue_count("invalid_barcode");
```

Export as:

### CSV file

```rust
info.report_to_csv("errors.tsv");
```

### String (tab-separated)

```rust
println!("{}", info.report_to_string());
```

---

# Built-in Read Counters

MappingInfo tracks structured counters:

- `quality`
- `length`
- `n_s`
- `poly_a`
- `no_sample`
- `no_data`
- `ok_reads`
- `cellular_reads`
- `multimapper`
- `pcr_duplicates`

Additionally, you can register custom read types:

```rust
info.iter_read_type("expression reads");
info.iter_read_type("antibody reads");
```

And print formatted summaries:

```rust
println!(
    "{}",
    info.read_types_to_string(vec![
        "expression reads",
        "antibody reads",
        "sample reads"
    ])
);
```

---

# Timing Measurements

MappingInfo contains structured timing support.

## Categories

- `single_processor_time`
- `multi_processor_time`
- `file_io_time`
- `subprocess_time`
- total runtime

## Example

```rust
info.start_counter();

// do work

info.stop_single_processor_time();
```

Split elapsed time:

```rust
let (h, m, s, ms) = info.elapsed_time_split();
```

Human-readable timing summary:

```rust
println!("{}", info.program_states_string());
```

---

# Logging Functionality

MappingInfo supports structured logging to file.

## Create with log file

```rust
let file = File::create("run.log").ok();
let mut info = MappingInfo::new(file, 20.0, 0);
```

## Write arbitrary log line

```rust
info.write_to_log("Pipeline started".to_string());
```

## Log report table

```rust
info.log_report();
```

## Periodic Progress Logging

Designed to work with `indicatif::ProgressBar`:

```rust
use indicatif::ProgressBar;

let pb = ProgressBar::new(100);

info.total += 1;
info.log(&pb);
```

This logs every `split` reads (default: 1,000,000).

---

# Histogram Support

```rust
info.iterate_hist(3);
```

Internally maintains a fixed-size histogram vector.

---

# Merging Results

Designed for multi-threaded aggregation:

```rust
info.merge(&other_info);
```

This merges:

- Counters
- Error reports
- Histograms
- Read logs

---

# Full Summary Output

```rust
let summary = info.summary(
    expression_umis,
    antibody_umis,
    sample_umis
);

println!("{}", summary);
```

Outputs:

- Read breakdown
- Filter reasons
- UMI statistics
- PCR duplicates
- Detailed timing
- Automatically writes to log file (if enabled)

---

# When To Use mapping_info

✔ BAM parsing tools\
✔ Single-cell quantifiers\
✔ High-performance Rust pipelines\
✔ Benchmarking Rust vs Python tools\
✔ Reproducible run summaries\
✔ HPC logging

---

# Design Goals

- Zero heavy dependencies
- Deterministic aggregation
- Thread-safe merge model
- Explicit timing separation
- Human-readable summaries
- Machine-readable exports

---

# Author

Stefan Lang\
Bioinformatics / Single-Cell Systems
