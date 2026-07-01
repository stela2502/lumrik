# bam-quant

`bam-quant` is the core quantification engine of the `bam_tide` toolkit.

The tool quantifies single-cell BAM alignments against a compiled splice index and produces sparse matrix outputs compatible with downstream single-cell analysis workflows.

The quantifier supports:

- gene-level quantification
- transcript-level quantification
- ONT long-read workflows
- multi-BAM merged analysis
- external read-tag tables
- and optional SNP-aware quantification

The tool is designed for direct integration into high-throughput sequencing pipelines and large-scale HPC workflows.

---

# Typical Use Cases

- ONT single-cell quantification
- transcript-level counting
- barcode-aware quantification
- SNP-aware expression analysis
- transcriptome-first workflows
- large sparse matrix generation

---

# Basic Usage

```bash
bam-quant \
  --bam mapped.bam \
  --index gencode.v49.annotation.gtf.dat \
  --outpath quant_out
```

---

# Transcript Quantification

```bash
bam-quant \
  --bam mapped.bam \
  --index gencode.v49.annotation.gtf.dat \
  --quant-mode transcript \
  --outpath quant_out
```

Supported modes:

- `gene`
- `transcript`

---

# Multi-BAM Quantification

Multiple BAM files may be quantified together in one run.

This is useful for:

- re-sequencing runs
- lane merges
- repeated ONT sequencing
- or workflow partitioning

Example:

```bash
bam-quant \
  --bam sample1.bam sample2.bam sample3.bam \
  --index gencode.v49.annotation.gtf.dat \
  --outpath merged_quant
```

All counts are accumulated into one shared output matrix.

---

# External Read-Tag Tables

Optional external read-tag tables may be supplied.

These tables map:

```text
read_id → cell barcode + UMI
```

and are useful when preprocessing pipelines preserve read names during alignment.

Example:

```bash
bam-quant \
  --bam sample1.bam sample2.bam \
  --read-tag-table sample1.tags.tsv sample2.tags.tsv \
  --index gencode.v49.annotation.gtf.dat \
  --outpath quant_out
```

BAM files and read-tag tables are paired by argument order.

---

# SNP-Aware Quantification

Optional SNP quantification can be enabled using:

```bash
--vcf variants.vcf
```

This generates additional sparse matrices containing:

- reference allele counts
- alternative allele counts

per cell.

SNP quantification requires:

```bash
--genome
```

because reads are refined against the genomic reference sequence.

Example:

```bash
bam-quant \
  --bam mapped.bam \
  --index gencode.v49.annotation.gtf.dat \
  --genome GRCh38.fa.gz \
  --vcf variants.vcf \
  --outpath quant_out
```

---

# Output Structure

Outputs are written in sparse matrix format compatible with many downstream single-cell workflows.

Typical outputs include:

```text
matrix.mtx.gz
features.tsv.gz
barcodes.tsv.gz
```

Additional SNP matrices may also be generated.

---

# Splice Matching Behavior

The quantifier performs splice-aware transcript matching using the compiled splice index.

Optional strictness settings include:

| Option | Description |
|---|---|
| `--require-strand` | Require strand-compatible matching |
| `--require-exact-junction-chain` | Require exact splice junction chains |
| `--max-5p-overhang-bp` | Allowed 5′ overhang |
| `--max-3p-overhang-bp` | Allowed 3′ overhang |
| `--allowed-intronic-gap-size` | Allowed sequencing error gap |

These settings are mainly intended for advanced tuning and experimental workflows.

---

# Intronic Splitting

Optional intronic splitting:

```bash
--split-intronic
```

attempts to separate intronic from exonic assignments.

This mode is currently considered experimental for many real sequencing datasets.

---

# Performance Notes

`bam-quant` is designed for:

- large BAM files
- parallel execution
- streaming quantification
- and memory-efficient sparse matrix generation

Threading is controlled using:

```bash
--threads 16
```

depending on available compute resources.

---

# Typical Workflow

```text
ONT BAM
  ↓
bam-ont-normalizer
  ↓
transcriptome alignment
  ↓
bam-transcriptome-to-genome
  ↓
bam-quant
  ↓
sparse matrices
```

---

# Important Inputs

## Splice Index

The splice index must first be generated using:

```text
gtf-splice-index
```

Example:

```bash
gtf-splice-index \
  --gtf gencode.v49.annotation.gtf \
  --out gencode.v49.annotation.gtf.dat
```

---

## Genome Refinement

When:

```bash
--genome
```

is supplied, reads may be refined against the genomic reference before quantification.

This improves SNP-aware processing and transcript matching in some workflows.

Genome refinement may be disabled explicitly:

```bash
--no-genome-refine
```

---

# Example Full Workflow

```bash
bam-quant \
  --bam results/*.genomic.bam \
  --read-tag-table results/*.tags.tsv \
  --index gencode.v49.annotation.gtf.dat \
  --genome GRCh38.p14.genome.fa.gz \
  --vcf SNPs.vcf \
  --quant-mode transcript \
  --threads 16 \
  --min-cell-counts 1 \
  --outpath merged_quant
```
