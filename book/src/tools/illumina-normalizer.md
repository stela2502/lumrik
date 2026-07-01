# illumina_normalizer

`illumina_normalizer` converts paired-end Illumina single-cell sequencing reads into a normalized FASTQ representation suitable for mapping and downstream quantification.

The tool detects primer-derived barcode structures using `sc_primer`, extracts cell barcode and UMI information, writes a mapper-facing FASTQ file, and records all recovered metadata in a read-tag table.

Unlike ONT data, Illumina single-cell libraries typically contain exactly one molecule per FASTQ pair. Therefore only the first valid primer match is used.

The resulting normalized FASTQ is designed for:

* mapping
* transcript quantification
* fusion detection
* mutation analysis
* custom downstream workflows

while preserving all cell barcode and UMI information in a sidecar file.

The normalized FASTQ is intended to be mapped using a conventional aligner such as minimap2 or STAR. The resulting BAM file can then be quantified using [bam-quant](bam-quant.md), which consumes the read-tag table to recover cell barcode, UMI, sample-tag, and feature-tag information while generating single-cell count matrices.


---

# Typical Use Cases

* 10x Genomics single-cell RNA-seq
* BD Rhapsody gene expression
* BD sample-tag experiments
* CRISPR screening libraries
* feature barcode experiments
* custom primer layouts
* normalization prior to mapping

---

# Typical Workflow

```text
Raw Illumina FASTQ
  ↓
illumina_normalizer
  ↓
normalized FASTQ
+
read_tags.tsv.gz
  ↓
mapper
  ↓
bam-quant
```

---

# Basic Usage

```bash
illumina_normalizer \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --out normalized.fastq.gz \
  --read-tags molecule_tags.tsv.gz \
  --chemistry tenx-three-prime-v3
```

---

# What Gets Extracted

The normalizer detects:

* cell barcodes
* UMIs
* sample tags
* feature tags
* insert sequence

depending on the selected chemistry.

Recovered barcode information is written to the read-tag table while the biological insert becomes the emitted FASTQ record.

---

# Standard 10x and BD Usage

Most Illumina single-cell libraries place barcode and UMI information in R1 and the biological insert in R2.

Example:

```bash
illumina_normalizer \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --out normalized.fastq.gz \
  --read-tags molecule_tags.tsv.gz \
  --primer-read r1 \
  --insert-read r2 \
  --chemistry bd-v2-384
```

Workflow:

```text
R1 → barcode extraction
R2 → mapper FASTQ
```

---

# Same-Read Primer and Insert Layouts

Some assays contain both primer structure and biological insert in the same sequencing read.

For these cases:

```bash
--insert-read detected
```

Example:

```bash
illumina_normalizer \
  --r1 input_R1.fastq.gz \
  --r2 input_R2.fastq.gz \
  --insert-read detected
```

The insert returned by the primer detector becomes the emitted FASTQ record.

---

# Primer Detection

Primer detection is provided by `sc_primer`.

Supported presets include:

| Chemistry            | Description               |
| -------------------- | ------------------------- |
| tenx-three-prime-v1  | 10x Chromium 3' v1        |
| tenx-three-prime-v2  | 10x Chromium 3' v2        |
| tenx-three-prime-v3  | 10x Chromium 3' v3 / v3.1 |
| tenx-three-prime-v4  | 10x Chromium 3' v4        |
| tenx-five-prime      | 10x Chromium 5'           |
| tenx-multiome-arc-v1 | 10x Multiome ARC          |
| bd-v1                | BD Rhapsody v1            |
| bd-v2-96             | BD Rhapsody 96-cell       |
| bd-v2-384            | BD Rhapsody 384-cell      |

Example:

```bash
--chemistry tenx-three-prime-v3
```

---

# Custom Primer Structures

Custom primer grammars can be supplied directly.

Example:

```bash
--primer-structure "CELL(16)UMI(12)POLYT(20+)INSERT"
```

When supplied, the custom grammar overrides any chemistry preset.

---

# Reverse Complement Detection

Some datasets contain reads in mixed orientation.

Optional reverse-complement searching:

```bash
--detect-reverse-complement
```

Both orientations are evaluated and the best valid primer match is selected.

---

# Sample Tags

Built-in BD sample-tag detection:

```bash
--species human
```

or

```bash
--species mouse
```

The corresponding sample-tag panel is loaded automatically.

---

# Feature Tags

Additional feature barcode references can be supplied as FASTA.

Example:

```bash
--tags feature_tags.fa
```

Typical applications:

* antibody-derived tags
* CRISPR guides
* hashing tags
* custom feature panels

---

# Read Tag Table

A metadata table is written for every accepted molecule.

Example fields:

| Column           | Description               |
| ---------------- | ------------------------- |
| read_id          | Emitted FASTQ identifier  |
| original_read_id | Original FASTQ identifier |
| orientation      | Match orientation         |
| raw_cb           | Cell barcode              |
| quality_cb       | Cell barcode qualities    |
| raw_umi          | UMI sequence              |
| quality_umi      | UMI qualities             |
| status           | Detection status          |

This table is typically consumed by downstream quantification tools.

---

# Minimum Insert Length

Very short inserts can be discarded.

Example:

```bash
--min-insert-len 30
```

Default:

```text
20 bp
```

---

# Example: BD Rhapsody

```bash
illumina_normalizer \
  --r1 BD_R1.fastq.gz \
  --r2 BD_R2.fastq.gz \
  --out normalized.fastq.gz \
  --read-tags read_tags.tsv.gz \
  --chemistry bd-v2-384 \
  --species mouse
```

---

# Example: 10x Chromium

```bash
illumina_normalizer \
  --r1 R1.fastq.gz \
  --r2 R2.fastq.gz \
  --out normalized.fastq.gz \
  --read-tags read_tags.tsv.gz \
  --chemistry tenx-three-prime-v3
```

---

# Output Files

Typical output:

```text
normalized.fastq.gz
read_tags.tsv.gz
```

The FASTQ contains only accepted biological inserts.

The read-tag table contains barcode and molecule metadata.

---

# Performance Notes

The normalizer streams paired FASTQ files directly.

Memory usage is dominated by:

* feature-tag indices
* sample-tag indices
* primer lookup structures

The FASTQ files themselves are not loaded into memory.

Threading accelerates:

* FASTQ decompression
* FASTQ compression
* molecule processing

using the `--threads` parameter.

---

# Typical Workflow Position

```text
Illumina FASTQ
  ↓
illumina_normalizer
  ↓
normalized FASTQ
+
read_tags.tsv
  ↓
mapping
  ↓
bam-quant
  ↓
single-cell count matrix
```

---

# Related Tools

```text
illumina_normalizer
    ↓
bam-quant

illumina_normalizer
    ↓
primer-restore

bam-ont-normalizer
    ↓
bam-quant
```

`primer-restore` can reconstruct the original primer-derived barcode structure from the normalized FASTQ and read-tag table when needed.
