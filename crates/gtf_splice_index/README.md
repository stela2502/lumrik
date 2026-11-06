# gtf_splice_index

> **Lumrik crate.** This crate is part of [Lumrik](../../README.md) and is distributed as part of the Lumrik workspace. See the [Lumrik README](../../README.md#license) and root `LICENSE` for the workspace licensing terms.

## What this crate does

Fast streaming GTF/GFF parsing and spatial/splice-aware genome annotation indexing. The parser is also used for annotations containing transposable-element records.

## Binaries

- `transcripts_at_position`

## Library use

Build a `SpliceIndex` from GTF/GFF and query it with `SplicedRead` objects; the parser/index types are the implementation truth used by downstream Lumrik crates.

## Novel splice-event names and frameshift annotation

During `bam_tide` novel-junction discovery, each distinct splice event is
represented by one name of the form:

`<Gene>*<chromosome>:<donor>-<acceptor>*del[_frameshift]`

For example: `TP53*chr17:7675994-7676272*del_frameshift` is the **format** of a
predicted frameshifting deletion; it is **not** a verified consequence for this
particular TP53 event. Coordinates use 0-based half-open genomic boundaries,
independent of transcription strand. A single junction can skip multiple exons
or part of an exon; it does not require full-length read coverage.

`gtf_splice_index::Transcript::coding_bases_removed_by_junction()` calculates
the intersection of the affected genomic interval with the reference
transcript's **coding exons** (CDS span clipped to annotated exons). Introns
are never counted. The `_frameshift` suffix is appended only if coding bases
are removed and their number is not divisible by three. A missing suffix does
**not** prove an in-frame protein: the reference may lack CDS annotation, or
the event may be noncoding. No translation or protein-effect prediction is
performed.

For speed, the first reference transcript used to create a new junction event
provides the CDS annotation and representative exon model. Subsequent reads
with the same gene/chromosome/junction reuse that event ID; other reference
isoforms may have different coding consequences. The suffix is therefore a
**reference-dependent prediction**, not a universally established effect.

## Detailed documentation



A flexible, streaming **GTF/GFF3 parser** plus a **splice-aware transcript index**
for fast matching of spliced reads against transcript models.

This crate is designed to sit between:

- **Alignment layer** (e.g. BAM → spliced blocks)
- **Annotation layer** (GTF/GFF transcripts)

and provide fast, deterministic transcript matching.

Coordinates are **0-based, half-open**: `[start, end)`.


## Installation

You can add `gtf_splice_index` directly from GitHub.

### Using `cargo add` (recommended)

```bash
cargo add gtf_splice_index --git https://github.com/stela2502/lumrik
```

This will add an entry like this to your `Cargo.toml`:

```toml
[dependencies]
gtf_splice_index = { git = "https://github.com/stela2502/lumrik" }
```

---

### Pin to a specific commit or branch (optional)

For reproducibility, you may want to pin a commit or branch.

**Specific commit:**

```toml
[dependencies]
gtf_splice_index = { git = "https://github.com/stela2502/lumrik", rev = "COMMIT_HASH" }
```

**Specific branch:**

```toml
[dependencies]
gtf_splice_index = { git = "https://github.com/stela2502/lumrik", branch = "main" }
```

---

### Local development (path dependency)

If you’re developing both crates together:

```toml
[dependencies]
gtf_splice_index = { path = "../gtf_splice_index" }
```

---


## Main workflow

1. Build a `SpliceIndex` from a GTF/GFF file.
2. Convert an aligned read into a `SplicedRead`.
3. Call `match_transcripts()` to obtain best transcript matches.

This is the intended **public entry point** of the crate.

---

## Core exported types

From the crate root:

```rust
pub use types::{RefBlock, SplicedRead, Strand};
```

### `SplicedRead`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicedRead {
    pub chr_id: usize,
    pub strand: Strand,
    pub blocks: Vec<RefBlock>,
    finalized: bool,
}
```

A `SplicedRead` is simply a list of aligned reference blocks.

Typically these come from a BAM CIGAR string.

---

## Building the index from a GTF/GFF

High-level API:

```rust
use gtf_splice_index::AnnotationBuilder;

// 1 Mb bins: good default for mammalian genomes
let index = AnnotationBuilder::new(1_000_000)
    .build_from_path("genes.gtf")
    .unwrap();
```

### Bin size parameter

The value passed to `AnnotationBuilder::new()` controls the **genomic bin size**
used for transcript bucketing.

Recommended values:

| Genome type | Bin size |
|-------------|----------|
| Human / mouse | 1,000,000 |
| Fly / yeast | 100,000 |
| Bacteria | 10,000–50,000 |

For most users: **1,000,000 is a safe default**.

---

## Streaming parser (advanced usage)

If you need low-level access to annotation records:

```rust
use std::fs::File;
use std::io::BufReader;
use gtf_splice_index::annotation::io::AnnotationReader;

let file = File::open("genes.gtf").unwrap();
let reader = BufReader::new(file);

let rdr = AnnotationReader::new(reader);
for rec in rdr.records() {
    let rec = rec.unwrap();
    println!("{} {}-{}", rec.seqname, rec.start0, rec.end0);
}
```

Most users should use `AnnotationBuilder` instead.

---

## MatchOptions: controlling what “match” means

Matching behavior is configured using `MatchOptions`.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchOptions {
    /// If true, require read blocks to be on a compatible strand.
    pub require_strand: bool,

    /// If true, require the read to have the exact same splice junction chain as the transcript.
    pub require_exact_junction_chain: bool,

    /// Maximum allowed 5′ overhang (bp). If exceeded -> OverhangTooLarge.
    pub max_5p_overhang_bp: u32,

    /// Maximum allowed 3′ overhang (bp). If exceeded -> OverhangTooLarge.
    pub max_3p_overhang_bp: u32,

    /// Allowed sequencing micro-gap inside exons.
    /// Gaps ≤ this size are ignored when computing junctions.
    /// Larger gaps become real junctions and may lead to JunctionMismatch.
    pub allowed_intronic_gap_size: u32,
}
```

### Semantics

- Overhangs are **always measured** and returned.
- Options only control whether a read is considered valid.
- Small internal gaps can be tolerated to avoid false junction mismatches.

---

## Matching a read against transcripts

The main API:

```rust
let hits = index.match_transcripts(&read, opts);
```

### Complete example

```rust
use gtf_splice_index::{
    AnnotationBuilder,
    RefBlock,
    SplicedRead,
    Strand,
    MatchOptions,
};

// 1) Build index
let index = AnnotationBuilder::new(1_000_000)
    .build_from_path("genes.gtf")
    .unwrap();

// 2) Build spliced read
let mut read = SplicedRead::new(
    1,
    Strand::Plus,
    vec![
        RefBlock::new(100, 150),
        RefBlock::new(200, 250),
    ],
);
read.finalize();

// 3) Configure matching
let opts = MatchOptions {
    require_strand: true,
    require_exact_junction_chain: false,
    max_5p_overhang_bp: 10,
    max_3p_overhang_bp: 10,
    allowed_intronic_gap_size: 5,
};

// 4) Match
let hits = index.match_transcripts(&read, opts);

for hit in hits {
    println!(
        "transcript={} class={:?} over5={} over3={}",
        hit.transcript.primary_name().unwrap_or("unknown"),
        hit.hit.class,
        hit.hit.overhang_5p_bp,
        hit.hit.overhang_3p_bp,
    );
}
```

---

## Match classes (conceptual overview)

Typical classifications:

- `ExactJunctionChain`  
  Read junction chain equals transcript junction chain.

- `Compatible`  
  Read junctions are a subset of transcript junctions.

- `JunctionMismatch`  
  Read contains junction(s) not present in transcript.

- `Intronic`  
  Read overlaps transcript span but includes intronic sequence.

- `OverhangTooLarge`  
  End overhang exceeds configured threshold.

- `StrandMismatch`  
  Strand incompatible (if required).

- `NoOverlap`  
  No genomic overlap with transcript.

The returned `MatchHit` always includes measured 5′ and 3′ overhangs.

---

## Typical integration pattern

1. Parse BAM record → build `SplicedRead`
2. Call `index.match_transcripts()`
3. Select best hit(s) by match class and overhangs
4. Assign transcript or gene label

---

## Targeted exon-skipping quantification (ONT use case; proposed workflow)

### Why this mode is useful

A splice-acceptor mutation can cause an exon to be skipped. For example, to
investigate **TP53 exon 5**, the biological question is whether individual ONT
molecules support **exons 4→5→6 (inclusion)** or **exons 4→6 (skipping)**.
Counting independent 4→5 and 5→6 junction features would allow a single long
read spanning 4→5→6 to contribute to both features. Instead, represent the
**two competing exon chains** as two synthetic transcript models and quantify
molecules assigned to those models.

This is intentionally **gene-specific**. The index should contain only the
requested gene and exon event, not all genes or all possible splice junctions.
It is especially useful when cell identifiers and UMIs are available in BAM
tags and downstream analysis needs per-cell counts in scData/H5AD.

### Intended usage

1. Choose a **specific reference transcript** for TP53; exon numbers are
   transcript-dependent. Resolve exons 4, 5, and 6 from the same GTF and
   record the reference transcript ID and genome assembly.
2. Construct a small GTF containing **two synthetic transcript records** with
   the same strand and actual exon coordinates:
   - `TP53_E4_5_6`: exon 4, exon 5, exon 6 (inclusion).
   - `TP53_E4_6`: exon 4, exon 6 (skipping).
3. Build a `SpliceIndex` from that GTF using `AnnotationBuilder`.
4. Match ONT alignments against the two models. Count only **positive,
   discriminating splice evidence**; do not count unspliced reads classified
   `Compatible` as proof of either isoform.
5. Reuse the downstream BAM quantifier's cell-barcode/UMI handling to produce
   per-cell inclusion and skipping counts. Reads with insufficient or
   conflicting evidence should remain unassigned, not silently counted twice.

The two synthetic transcripts are **event models**, not claims that these are
complete biological TP53 transcript isoforms. An ONT read may also contain
splice junctions outside the event.

### Existing matcher behavior and limitations

`MatchOptions::require_exact_junction_chain = true` makes the matcher return
`ExactJunctionChain` when a read's **entire** observed junction chain equals
that of the synthetic transcript. The `Compatible` class is still possible
for reads without observed junctions; downstream quantification **must filter
for `ExactJunctionChain`** if using this strict mode.

Strict equality can exclude otherwise informative ONT reads with additional
upstream or downstream junctions. A future **event-chain matching** mode should
instead test the local junction chain around the targeted exon while allowing
extra junctions outside the event, and require unique assignment to one model.
This event-local matching mode and automatic synthetic-GTF generation are
**design proposals, not currently implemented CLI features**.

Large `max_5p_overhang_bp` and `max_3p_overhang_bp` values (for example,
10,000 bp) relax **outer transcript boundary** checks only; they do not relax
internal splice-junction coordinates. They also do not by themselves guarantee
that additional ONT exonic blocks outside a synthetic model will pass the
existing exon-block compatibility check. Validate this behavior before using
long-read counts for biological interpretation.

**Suggested output:** one sparse cell × {`TP53_E4_5_6`, `TP53_E4_6`}
count matrix, with event/transcript metadata and optional read counts alongside
UMI-deduplicated molecule counts. The skipping fraction is
`skip / (include + skip)` when the denominator is nonzero; otherwise it is
missing, not zero.

---

## Performance characteristics

- Streaming annotation parser
- Binned transcript index for fast candidate lookup
- Junction-based transcript filtering
- Deterministic matching results

---

## License

Add your license here (e.g. MIT / Apache-2.0).
