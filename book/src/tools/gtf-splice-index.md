# gtf-splice-index

`gtf-splice-index` builds and inspects serialized splice indexes derived from transcript annotations.

The generated index is a compact binary representation of:

- transcripts
- exons
- splice structures
- gene metadata
- transcript metadata
- and genomic interval bins

The index is designed for fast splice-aware matching during:

- quantification
- transcript projection
- and long-read processing workflows

within the `bam_tide` ecosystem.

---

# Typical Use Cases

- preprocessing transcript annotations
- building reusable splice indexes
- transcript-aware quantification
- splice-aware matching
- transcriptome projection workflows

---

# Basic Workflow

```text
GTF/GFF
  ↓
gtf-splice-index build
  ↓
serialized splice index
  ↓
bam-quant
bam-transcriptome-to-genome
```

---

# Building an Index

```bash
gtf-splice-index build \
  --annotation gencode.v49.annotation.gtf \
  --index gencode.v49.annotation.gtf.dat
```

---

# Inspecting an Index

```bash
gtf-splice-index stats \
  --index gencode.v49.annotation.gtf.dat
```

This prints summary statistics describing the serialized index contents.

---

# Input Formats

The builder supports:

- GTF
- GFF
- GFF3

annotation files.

---

# Output

The generated index is written as a serialized binary file:

```text
gencode.v49.annotation.gtf.dat
```

This file can later be loaded directly by:

- `bam-quant`
- `bam-transcriptome-to-genome`
- and related workflows

without repeatedly parsing the original annotation.

---

# Example

```bash
gtf-splice-index build \
  --annotation stringtie.gff \
  --index stringtie.gff.dat
```

---

# Annotation Customization

The builder supports configurable annotation keys for compatibility with different GTF/GFF conventions.

Examples include:

| Option | Meaning |
|---|---|
| `--gene-id-key` | Gene ID attribute |
| `--gene-name-key` | Gene name attribute |
| `--transcript-id-key` | Transcript ID attribute |
| `--transcript-name-key` | Transcript name attribute |
| `--parent-key` | GFF3 exon→transcript linkage |
| `--exon-feature-type` | Features treated as exons |

This improves compatibility with:

- Gencode
- Ensembl
- StringTie
- and custom transcript annotations

---

# Bin Width

```bash
--bin-width 1000000
```

controls genomic interval binning used internally for efficient overlap lookup.

The default value works well for most workflows.

Advanced users may tune this depending on:

- genome size
- transcript density
- and workload characteristics

---

# Performance Notes

Building the splice index is usually a one-time preprocessing step.

The serialized index is optimized for:

- fast loading
- repeated reuse
- low-overhead transcript lookup
- and splice-aware matching

during downstream analysis.

---

# Typical Workflow Context

```text
annotation.gtf
  ↓
gtf-splice-index build
  ↓
annotation.gtf.dat
  ↓
bam-quant
```

or:

```text
annotation.gtf
  ↓
gtf-splice-index build
  ↓
annotation.gtf.dat
  ↓
bam-transcriptome-to-genome
```

---

# Important Notes

## Transcript IDs must remain consistent

The transcript IDs stored in the splice index must match:

- transcriptome FASTA identifiers
- transcriptome BAM reference names
- and downstream workflow annotations

Inconsistent transcript naming is one of the most common causes of projection and quantification failures.

---

## Reusable binary format

The serialized index is intended to replace repeated parsing of large text-based annotations during high-throughput workflows.

This significantly improves startup time and scalability for large sequencing datasets.
