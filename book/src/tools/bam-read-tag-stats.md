# bam-read-tag-stats

`bam-read-tag-stats` summarizes external read-tag tables produced by:

- `bam-ont-normalizer`
- or compatible preprocessing workflows

The tool provides lightweight QC and exploratory statistics for:

- cell barcode usage
- UMI complexity
- barcode/UMI pair frequencies
- and filtering behavior

It is primarily intended for validation and inspection of preprocessing workflows before quantification.

---

# Typical Use Cases

- QC of ONT preprocessing
- barcode complexity inspection
- UMI distribution analysis
- debugging read-tag generation
- threshold tuning
- validation of preprocessing pipelines

---

# Basic Usage

```bash
bam-read-tag-stats \
  --read-tag-table molecule_tags.tsv
```

Compressed input is supported automatically:

```bash
bam-read-tag-stats \
  --read-tag-table molecule_tags.tsv.gz
```

---

# Multiple Tables

Multiple read-tag tables may be analyzed together:

```bash
bam-read-tag-stats \
  --read-tag-table run1.tsv.gz \
  --read-tag-table run2.tsv.gz
```

This is useful for:

- lane merges
- repeated sequencing
- multi-run experiments
- and workflow comparisons

---

# Filtering Parameters

The tool applies lightweight filtering to reduce noise from weak barcode/UMI observations.

---

## Minimum pair count

```bash
--min-pair-count 2
```

Minimum number of observations required for one:

```text
cell barcode + UMI
```

combination.

---

## Minimum UMIs per cell

```bash
--min-cell-umis 3
```

Minimum number of surviving UMIs required for a barcode after pair filtering.

---

# Input Format

The expected input is a TSV file containing:

- read IDs
- cell barcodes
- UMIs
- and optional quality/provenance columns

Typical source:

```text
bam-ont-normalizer
```

Compressed `.gz` files are supported automatically.

---

# Configurable Columns

The parser supports configurable column names.

Examples:

| Option | Meaning |
|---|---|
| `--rt-read-id-column` | Read identifier column |
| `--rt-cell-column` | Cell barcode column |
| `--rt-umi-column` | UMI column |
| `--rt-original-read-id-column` | Original read provenance |

This allows compatibility with external preprocessing pipelines.

---

# Example

```bash
bam-read-tag-stats \
  --read-tag-table molecule_tags.tsv.gz \
  --min-pair-count 3 \
  --min-cell-umis 5
```

---

# Typical Workflow Position

```text
ONT preprocessing
  ↓
read-tag table generation
  ↓
bam-read-tag-stats
  ↓
QC and threshold tuning
  ↓
bam-quant
```

---

# Current Status

The tool is intended as a lightweight QC and debugging utility.

Future versions may include:

- barcode correction statistics
- whitelist matching summaries
- UMI complexity metrics
- collision estimation
- sequencing saturation estimates
- and additional cell-level QC metrics
