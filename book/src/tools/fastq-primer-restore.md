# primer-restore

`primer-restore` reconstructs complete single-cell sequencing reads from normalized insert FASTQ files.

The tool restores the primer-derived barcode structure that was previously removed during normalization.

For each FASTQ record:

* the read ID is looked up in a read-tag table
* cell barcode information is recovered
* UMI information is recovered
* sample-tag information is recovered (if present)
* a primer sequence is regenerated using `sc_primer`
* the reconstructed primer is prepended to the insert sequence
* a complete FASTQ record is emitted

The resulting FASTQ can be processed as if the original single-cell primer structure had never been removed.

---

# Typical Use Cases

* restoring 10x-style reads after normalization
* restoring BD Rhapsody reads after normalization
* rebuilding primer structure for downstream tools
* recreating test datasets
* validating primer extraction workflows
* round-trip testing of normalization pipelines

---

# Typical Workflow

```text
ONT / Illumina reads
  ↓
bam-ont-normalizer
  ↓
normalized FASTQ
+
read_tags.tsv.gz
  ↓
primer-restore
  ↓
restored FASTQ
```

---

# Basic Usage

```bash
primer-restore \
  --fastq normalized.fastq.gz \
  --read-tags read_tags.tsv.gz \
  --out restored.fastq.gz \
  --chemistry bd-v2-384
```

---

# Input Files

The tool requires:

* a normalized FASTQ file
* a read-tag table generated during normalization

Example:

```bash
primer-restore \
  --fastq normalized.fastq.gz \
  --read-tags molecule_tags.tsv.gz \
  --out restored.fastq.gz
```

---

# Read Tag Table

The read-tag table stores the information required to reconstruct the original primer structure.

Typical fields include:

| Field       | Description                   |
| ----------- | ----------------------------- |
| read_id     | FASTQ read identifier         |
| raw_cb      | Cell barcode                  |
| quality_cb  | Cell barcode qualities        |
| raw_umi     | UMI sequence                  |
| quality_umi | UMI qualities                 |
| orientation | Original molecule orientation |
| status      | Extraction status             |

The exact columns depend on the normalizer version.

---

# Primer Reconstruction

Primer generation is performed through `sc_primer`.

The selected chemistry determines:

* barcode layout
* UMI layout
* adapter structure
* primer sequence
* chemistry-specific formatting

Examples:

```bash
--chemistry tenx-three-prime-v3
```

```bash
--chemistry bd-v2-384
```

```bash
--chemistry bd-v2-96
```

---

# Error Handling

By default the tool operates in strict mode.

If a FASTQ record cannot be matched to a read-tag entry:

```text
read skipped
error recorded
```

If primer reconstruction fails:

```text
read skipped
error recorded
```

The program exits with a non-zero return code whenever restoration errors occurred.

This prevents accidental generation of incomplete datasets.

---

# Allow Partial Output

Optional relaxed mode:

```bash
--allow-partial
```

In this mode:

* successfully restored reads are written
* failed reads are skipped
* warnings are emitted
* execution completes successfully

Useful for debugging and recovery workflows.

---

# Example

```bash
primer-restore \
  --fastq normalized.fastq.gz \
  --read-tags read_tags.tsv.gz \
  --out restored.fastq.gz \
  --chemistry tenx-three-prime-v3
```

---

# Output

The output FASTQ contains:

```text
restored primer
+
restored cell barcode
+
restored UMI
+
original insert
```

All qualities are reconstructed and preserved where available.

---

# Summary Report

At completion a restoration summary is printed.

Example:

```text
primer-restore summary
======================
FASTQ reads processed      : 1000000
FASTQ reads restored       : 999872
reads without tag record   : 103
primer generation failures : 25
```

If restoration errors occurred:

```text
WARNING: primer-restore did not restore all reads.
Output FASTQ may be incomplete.
```

---

# Performance Notes

The tool streams FASTQ records directly.

Memory usage is dominated by:

* the read-tag lookup table
* primer reconstruction metadata

The FASTQ itself is processed sequentially and is not loaded entirely into memory.

---

# Typical Workflow Position

```text
Original sequencing reads
  ↓
normalization
  ↓
normalized FASTQ
+
read tags
  ↓
primer-restore
  ↓
restored FASTQ
  ↓
single-cell analysis
```
