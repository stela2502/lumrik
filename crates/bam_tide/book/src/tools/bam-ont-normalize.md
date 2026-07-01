# bam-ont-normalizer

`bam-ont-normalizer` converts noisy ONT/Dorado BAM reads into normalized single-cell FASTQ molecules.

The tool searches each read in both orientations for a 10x-style cassette structure:

```text
adapter + cell barcode + UMI + polyT + transcript
```

Detected molecules are extracted, orientation-normalized, and written as independent FASTQ records suitable for transcriptome alignment workflows.

In addition to FASTQ output, the tool generates a TSV sidecar file containing:

- extracted cell barcodes
- UMIs
- quality strings
- molecule coordinates
- orientation
- and extraction status

The normalizer does not perform barcode correction or whitelist matching. All CB/UMI values are emitted exactly as observed in the read.

---

# Typical Use Cases

- ONT single-cell preprocessing
- Dorado BAM normalization
- barcode and UMI extraction
- transcriptome-first ONT workflows
- long-read single-cell preprocessing

---

# Basic Usage

```bash
bam-ont-normalizer \
  --bam dorado.bam \
  --out normalized.fastq.gz \
  --tags molecule_tags.tsv
```

---

# Typical Workflow

```text
Dorado BAM
  ↓
bam-ont-normalizer
  ↓
normalized FASTQ
  ↓
transcriptome alignment
  ↓
bam-transcriptome-to-genome
  ↓
bam-quant
```

---

# Output FASTQ Structure

Each detected molecule becomes one FASTQ entry:

```text
original_read_name/mol<N>
```

The emitted sequence is normalized into the expected transcript orientation.

Reads without detectable cassettes are not emitted but remain included in the final summary statistics.

---

# Read-Tag TSV Output

The sidecar TSV contains one row per emitted molecule.

Typical columns include:

- read ID
- original read ID
- orientation
- raw cell barcode
- barcode qualities
- raw UMI
- UMI qualities
- extraction status

This table can later be reused during quantification using:

```bash
bam-quant --read-tag-table
```

---

# Adapter and Cassette Detection

The tool searches for:

```text
adapter + CB + UMI + polyT
```

using configurable matching rules.

Default adapter:

```text
CTACACGACGCTCTTCCGATCT
```

Default structure:

| Component | Default Length |
|---|---|
| Cell barcode | 16 bp |
| UMI | 12 bp |

---

# Example

```bash
bam-ont-normalizer \
  --bam dorado.bam \
  --out normalized.fastq.gz \
  --tags molecule_tags.tsv \
  --threads 8 \
  --gzip-level 1
```

---

# Important Parameters

## Adapter Matching

```bash
--min-adapter-match 13
```

Minimum adapter suffix length required for cassette detection.

Useful for noisy ONT reads where only part of the adapter is present.

---

## PolyT Detection

```bash
--poly-t-min 10
--poly-t-window 14
```

Control expected polyT validation behavior after:

```text
adapter + CB + UMI
```

---

## Minimum Transcript Length

```bash
--min-transcript-len 20
```

Minimum transcript sequence length required after the polyT region.

Shorter molecules are counted but not emitted.

---

# Compression

FASTQ output is gzip-compressed by default.

Compression level:

```bash
--gzip-level 1
```

is recommended for large ONT datasets because it significantly improves throughput.

Plain FASTQ output may be enabled using:

```bash
--no-gzip
```

---

# Performance Notes

The tool is optimized for:

- streaming BAM processing
- ONT-scale datasets
- and large preprocessing workloads

Threading currently accelerates:

- BAM decompression
- BAM reading
- and FASTQ compression

depending on storage throughput.

---

# Important Notes

## No barcode correction

The tool intentionally does not perform:

- whitelist correction
- barcode correction
- UMI collapsing
- or cell calling

All sequences are emitted exactly as observed.

This keeps preprocessing reproducible and transparent.

---

## Orientation normalization

Both read orientations are scanned automatically.

Detected molecules are emitted in normalized transcript orientation regardless of original read direction.

---

# Example Full Workflow

```bash
bam-ont-normalizer \
  --bam dorado.bam \
  --out normalized.fastq.gz \
  --tags molecule_tags.tsv \
  --threads 8

minimap2 \
  -ax map-ont transcriptome.fa normalized.fastq.gz \
  | samtools view -b -o mapped_tx.bam

bam-transcriptome-to-genome \
  --bam mapped_tx.bam \
  --gtf transcripts.gtf \
  --genome GRCh38.fa.gz \
  --out mapped_genome.bam

bam-quant \
  --bam mapped_genome.bam \
  --read-tag-table molecule_tags.tsv \
  --index gencode.v49.annotation.gtf.dat \
  --outpath quant_out
```
