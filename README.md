[![Rust](https://github.com/stela2502/int_to_str/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/int_to_str/actions/workflows/rust.yml)
[![Rust](https://github.com/stela2502/int_to_str/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/int_to_str/actions/workflows/rust.yml)

# int_to_str (v0.2.0)

A lightweight Rust library for encoding DNA sequences into compact 2-bit representations.

This crate allows efficient conversion between DNA strings and integer-based representations (`u16`, `u32`, `u64`, `u128`), while also supporting slicing and low-level sequence manipulation.

> ⚠️ The implementation prioritizes performance and compactness. It is designed for practical use in bioinformatics workflows rather than being a fully polished general-purpose library.

---

# Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
int_to_str = { git = "https://github.com/stela2502/int_to_str.git", version = "0.2.0" }
```

---

# Basic Usage

Encode a DNA sequence into the internal 2-bit representation:

```rust
use int_to_str::int_to_str::IntToStr;

let tool = IntToStr::new(
    b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACTACCTCTGGGCCTGGGATAC"
);
```

You can:

* Convert to integers (`u16`, `u32`, `u64`, `u128`)
* Slice sequences efficiently
* Convert back to string representation

---

# Binary Usage

The crate includes a small CLI tool:

```bash
int_to_str <sequence or integer>
```

### Integer → DNA

```bash
int_to_str 12343
Integer input: 12343
→ Sequence: TCTAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

### DNA → Integer

```bash
int_to_str TCTAAAT
Sequence input: TCTAAAT
→ Sequence: 12343
```

⚠️ Only the **significant bits** correspond to real DNA.
The binary does not track sequence length for you, so trailing `A`s (00) may appear.

---

# Optional Feature: Alignment (v0.2.0)

Version **0.2.0** introduces an optional `alignment` feature that adds sequence comparison utilities directly on top of the `IntToStr` representation.

## Enable the feature

```toml
[dependencies]
int_to_str = { git = "https://github.com/stela2502/int_to_str.git", version = "0.2.0", features = ["alignment"] }
```

or:

```bash
cargo build --features alignment
```

---

## What it adds

When enabled, the following capabilities become available:

### Exact overlap detection

* `best_exact_overlap(...)`
* `exact_overlap_with_offset(...)`

Finds the best positional overlap between two sequences, including:

* overlap length
* mismatch count
* relative offset

---

### Mismatch-aware comparison

* `overlap_mismatches(...)`

Allows fast and deterministic comparison of overlapping regions without full alignment.

---

### Needleman–Wunsch alignment

* `needleman_wunsch(...)`
* `needleman_wunsch_distance(...)`

Provides global alignment and a normalized distance metric suitable for:

* validating overlaps
* rescuing near-matching sequences
* comparing short DNA fragments

---

## Example

```rust
use int_to_str::IntToStr;

let a = IntToStr::new(b"AAGCAGTGGTATCAACGC");
let b = IntToStr::new(b"TGGTATCAACGCAGAGTAA");

// Exact overlap detection
if let Some(hit) = a.best_exact_overlap(&b, 10, 0.0) {
    println!("Overlap: {:?}", hit);
}

// Alignment-based similarity
let dist = a.needleman_wunsch_distance(&b);
println!("Distance: {}", dist);
```

---

## When to use this feature

Enable `alignment` if you need:

* robust overlap detection (e.g. primer reconstruction)
* mismatch-aware sequence comparison
* alignment-based validation of sequence relationships

You likely **do not need it** if you only use:

* encoding / decoding
* integer conversion
* slicing operations

---

## Design Notes

* The alignment functionality is **feature-gated** to keep the core crate lightweight.
* All operations work directly on the **2-bit encoded representation** (no string conversion overhead).
* Optimized for **short sequences (≤32 bp)**, making it ideal for:

  * primers
  * k-mers
  * short-read analysis

---

# Final Notes

This crate is intentionally minimal and pragmatic.
If you need full-scale bioinformatics tooling, you should likely combine it with other libraries.

If you need **fast, compact DNA handling with optional alignment support**, this crate is exactly that.
