# bam-ont-normalizer

## Overview

`bam-ont-normalizer` converts raw Oxford Nanopore (ONT) Dorado BAM reads into clean, 10x-style FASTQ molecules.

It reconstructs canonical single-cell molecules of the form:

```
adapter + cell barcode (CB) + UMI + polyT + transcript
```

Each detected molecule is emitted as a separate FASTQ record.

This tool is designed as a preprocessing step for pipelines such as `wf-single-cell`, which assume well-formed 10x read structure.

This statement is currently checked...

---

## Motivation

ONT reads frequently contain:

* multiple concatenated molecules per read
* reverse-complement orientation
* truncated adapters
* noisy sequence around the barcode cassette

As a result, many reads are discarded during barcode extraction.

`bam-ont-normalizer` addresses this by:

* detecting partial adapter matches
* enforcing strict CB/UMI geometry
* normalizing each molecule into canonical orientation
* splitting multi-cassette reads into multiple FASTQ entries

---

## Input / Output

### Input

* Dorado BAM (mapped or unmapped)
* Requires query sequence and qualities

### Output

1. **Normalized FASTQ**

Each record:

```
@<original_read_id>/mol<N>
<adapter+CB+UMI+polyT+transcript>
+
<quality>
```

2. **Molecule metadata TSV (optional)**

Per molecule:

* original read ID
* molecule index
* orientation
* CB / UMI (raw)
* quality slices
* coordinates
* polyT statistics

---

## Core Algorithm

For each BAM record:

1. Extract sequence and quality
2. Evaluate both orientations:

   * forward
   * reverse complement
3. Scan for adapter suffix matches:

   * allow mismatches
   * allow truncated prefix
4. Infer full adapter position from suffix match
5. Enforce strict structure:

```
adapter_end → CB(16) → UMI(12) → polyT
```

6. Validate polyT:

   * minimum T count within window
   * exact expected position
7. Define molecule boundaries:

   * from inferred adapter start
   * to next cassette or read end
8. Emit one FASTQ record per valid cassette

---

## Key Design Principles

### 1. Strict CB/UMI Geometry

Cell barcode and UMI positions are **not fuzzy matched**.

Reads with indels or shifts in CB/UMI are discarded.

> Rationale: incorrect CB/UMI is worse than losing the read as it ould not be counted later on anyhow.

### 2. Flexible Adapter Matching

Adapter detection tolerates:

* truncated prefix
* limited mismatches
* ONT sequencing noise

Adapter matches are used only to **anchor structure**, not as final output.

### 3. PolyT as Structural Anchor

PolyT detection ensures correct placement of CB/UMI relative to transcript.

### 4. Multi-cassette Splitting

Reads containing multiple molecules are split:

```
read → mol1, mol2, mol3, ...
```

---

## Example Statistics

Typical output summary:

```
reads processed        : 10,496
reads w/ cassette      : 5,311 (50.6%)
reads w/o cassette     : 5,185

molecules emitted      : 5,365
mean molecules/read    : 0.51

multi-cassette reads   : 52 (0.5%)
forward / reverse      : ~50 / 50
```

Interpretation:

* ~50% of ONT reads contain usable 10x structure
* multi-molecule reads are rare (~1%)
* orientation is balanced

---

## Usage

### Basic

```
bam-ont-normalizer \
  --bam input.bam \
  --out normalized.fastq.gz \
  --tags molecule_tags.tsv
```

### High-performance streaming

This tool is horribly inefficient at the moment.
Only if the assumptions can be proven it will be developed further.

```
bam-ont-normalizer \
  --bam input.bam \
  --out - \
  --no-tags \
| pigz -p 8 -1 > normalized.fastq.gz
```

---

## Important Limitations

* No barcode whitelist correction
* CB/UMI must match exact expected positions
* Indels in CB/UMI → read is discarded
* Assumes 10x 3' chemistry structure
* Not designed for full-length transcript protocols

---

## When to Use

Use this tool when:

* working with ONT single-cell data
* Dorado BAM contains raw read sequence
* downstream pipeline expects strict 10x structure

Do **not** use when:

* CB/UMI already extracted reliably
* working with Illumina data
* reads lack polyT tail

---

## Citation

If you use this tool, please cite via the repository:

```
Lang, S. bam-ont-normalizer. GitHub.
https://github.com/stela2502/bam_tide
```

---

## Status

This is an experimental preprocessing tool.

Biological validation should be performed downstream (e.g., cell recovery, UMI counts, clustering consistency).

---

## Summary

`bam-ont-normalizer` reconstructs usable single-cell molecules from noisy ONT reads by combining:

* flexible adapter detection
* strict structural enforcement
* per-molecule normalization

It is intended as a pragmatic bridge between ONT data and 10x-oriented analysis pipelines.

