# bam-coverage

`bam-coverage` generates genomic coverage tracks from BAM files.

The tool is designed as a lightweight and high-performance alternative to Python-based coverage exporters such as:

```text
deeptools bamCoverage
```

It supports:

- coverage binning
- multiple normalization strategies
- SAM flag filtering
- mapping-quality filtering
- and BigWig or bedGraph export workflows

depending on the selected output format and build configuration.

---

# Typical Use Cases

- genome browser visualization
- coverage track generation
- QC and inspection
- ChIP-seq style signal generation
- RNA-seq coverage visualization
- benchmarking against deeptools

---

# Basic Usage

```bash
bam-coverage \
  --bam sample.bam \
  --outfile sample.bw
```

---

# Example

```bash
bam-coverage \
  --bam sample.bam \
  --outfile sample.bw \
  --normalize cpm \
  --width 50
```

---

# Input Requirements

The input BAM does not need to be:

- sorted
- indexed
- or coordinate-ordered

The tool processes BAM records directly.

---

# Coverage Normalization

Supported normalization modes:

| Mode | Description |
|---|---|
| `not` | Raw read counts |
| `rpkm` | Reads Per Kilobase per Million mapped reads |
| `cpm` | Counts Per Million mapped reads |
| `bpm` | Bins Per Million mapped reads |
| `rpgc` | Reads Per Genomic Content (1× genome coverage) |

Example:

```bash
--normalize cpm
```

---

# Bin Width

Coverage is calculated in fixed-width bins.

Example:

```bash
--width 50
```

Smaller bins produce higher-resolution tracks but increase output size and runtime.

---

# Mapping Quality Filtering

Reads below the specified MAPQ threshold are ignored.

Example:

```bash
--min-mapping-quality 20
```

---

# SAM Flag Filtering

The tool supports deeptools-style SAM flag filtering.

---

## Excluding Reads

Example:

```bash
--sam-flag-exclude 2816
```

This excludes:

- secondary alignments
- QC-failed reads
- supplementary alignments

A read is excluded if:

```text
(read_flag & mask) != 0
```

---

## Including Reads

Example:

```bash
--sam-flag-include 64
```

This keeps only reads matching:

```text
read1
```

A read is included only if:

```text
(read_flag & mask) == mask
```

---

# Duplicate and Supplementary Reads

Optional inclusion flags:

| Option | Meaning |
|---|---|
| `--include-secondary` | Include secondary alignments |
| `--include-supplementary` | Include supplementary alignments |
| `--include-duplicates` | Include duplicate-marked reads |

By default these are excluded to match common coverage-generation conventions.

---

# Read1-Only Coverage

```bash
--only-r1
```

restricts coverage generation to read1 alignments.

Useful for:

- paired-end single-cell workflows
- duplicate reduction
- and compatibility with certain downstream analyses

---

# Example Comparison Workflow

```text
BAM
  ↓
bam-coverage
  ↓
BigWig
  ↓
bw-compare
  ↓
validation report
```

---

# Performance Notes

The tool is optimized for:

- streaming BAM traversal
- low-overhead coverage accumulation
- large genomic datasets
- and high-throughput workflows

Compared to Python-based coverage generation, the main performance improvements typically come from:

- reduced serialization overhead
- direct memory management
- and compiled execution

rather than algorithmic shortcuts.

---

# Typical Workflow Context

```text
aligned BAM
  ↓
bam-coverage
  ↓
BigWig / bedGraph
  ↓
genome browser visualization
```

---

# Example Full Workflow

```bash
bam-coverage \
  --bam aligned.bam \
  --outfile aligned.bw \
  --normalize cpm \
  --width 50 \
  --min-mapping-quality 20 \
  --sam-flag-exclude 2816
```
