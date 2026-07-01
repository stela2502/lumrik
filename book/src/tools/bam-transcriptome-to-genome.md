# bam-transcriptome-to-genome

`bam-transcriptome-to-genome` converts transcriptome-aligned BAM records back into genomic coordinates.

The tool projects transcript-space alignments into genome-space using:

- a transcript annotation (GTF/GFF)
- and the genome FASTA used to construct the transcriptome reference

This is especially useful for transcriptome-first workflows where reads are initially aligned against transcript sequences for improved sensitivity or simplified long-read alignment.

The resulting BAM can then be:

- visualized in genome browsers
- used for SNP-aware analysis
- processed by genome-coordinate tools
- or quantified using downstream workflows

---

# Typical Use Cases

- transcriptome-first ONT workflows
- projection of spliced transcript alignments
- genome-coordinate visualization
- SNP-aware downstream analysis
- preparation for genome-based quantification

---

# Basic Usage

```bash
bam-transcriptome-to-genome \
  --bam transcriptome.bam \
  --gtf transcripts.gtf \
  --genome GRCh38.fa.gz \
  --out genomic.bam
```

---

# Typical Workflow

```text
normalized FASTQ
  ↓
transcriptome alignment
  ↓
transcriptome BAM
  ↓
bam-transcriptome-to-genome
  ↓
genome-coordinate BAM
  ↓
bam-quant
```

---

# Input Requirements

## Transcript IDs must match

The transcript IDs present in the annotation must match the transcript names used in the transcriptome BAM.

Example:

```text
ENST00000335137
```

must appear consistently in:

- the BAM reference names
- the transcriptome FASTA
- and the annotation

---

## Matching genome FASTA

The supplied genome FASTA must match the genome build used to generate the transcriptome reference.

Chromosome naming should also be compatible.

Example mismatches:

```text
chr1
```

vs.

```text
1
```

may require chromosome mapping support or matching references.

---

# Output

The tool produces a genome-coordinate BAM preserving:

- splice structure
- transcript orientation
- transcript provenance
- and genomic positioning

Spliced transcript alignments are converted into genomic splice-aware CIGAR strings.

---

# Example

```bash
bam-transcriptome-to-genome \
  --bam mapped_tx.bam \
  --gtf stringtie.gff \
  --genome GRCh38.p14.genome.fa.gz \
  --out mapped_genome.bam \
  --threads 8
```

---

# Performance Notes

Large transcriptome BAM files can be processed efficiently using multiple threads:

```bash
--threads 16
```

depending on:

- storage throughput
- BAM compression
- and filesystem performance

Threading primarily accelerates:

- BAM reading
- BAM decompression
- BAM writing
- BGZF compression

---

# Quiet Mode

```bash
--quiet
```

reduces logging output while still reporting failed or skipped-read statistics.

Useful for large production workflows.

---

# Typical Workflow Context

This tool is mainly intended for workflows such as:

```text
ONT reads
  ↓
bam-ont-normalizer
  ↓
transcriptome alignment
  ↓
bam-transcriptome-to-genome
  ↓
bam-quant
```

rather than direct genomic alignment workflows.

---

# Important Notes

## Transcriptome-first workflows

The tool assumes reads were aligned against transcript sequences rather than genomic chromosomes.

It reconstructs genomic coordinates by projecting transcript exon structures back onto the genome.

---

## Long-read compatibility

The tool is designed for:

- ONT long reads
- spliced transcript alignments
- transcriptome-derived mapping workflows
- and multi-exon transcript structures

including noisy long-read sequencing data.

---

# Example Full Workflow

```bash
bam-ont-normalizer \
  --bam dorado.bam \
  --out normalized.fastq.gz \
  --tags molecule_tags.tsv

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
