# BAM Dataset Summary

**File:** `/home/med-sal/NAS/David.Bryder/UMI_Test/GMLP_S1_R2_001.AGATAG.fastq.gz.sorted.bam`  
**Generated:** 2026-02-04

---

## Overview

This BAM file contains alignments generated with **HISAT (v0.1.5-beta)** and processed with **samtools (v1.19.2)**.  
The file is **coordinate-sorted** and indexed.

It represents a **single-end sequencing experiment** (no paired reads present).

---

## Reference Genome

The BAM is aligned against a mouse-like reference genome with the following chromosomes:

- Autosomes: chr1–chr19  
- Sex chromosomes: chrX, chrY  
- Mitochondrial: chrM  

Example lengths:

- chr1: 195,471,971 bp  
- chrX: 171,031,299 bp  
- chrM: 16,299 bp  

Total references: **22 chromosomes**.

---

## Read Statistics

| Metric                   | Value                   |
| ------------------------ | ----------------------- |
| Total reads              | **20,178,657**          |
| Mapped reads             | **13,074,466 (64.79%)** |
| Primary alignments       | **16,286,433**          |
| Secondary alignments     | **3,892,224**           |
| Supplementary alignments | **0**                   |
| Duplicates               | **0**                   |
| Paired reads             | **0**                   |
| Properly paired reads    | **0**                   |

---

## Chromosome Coverage (Top 15 by mapped reads)

| Chromosome | Mapped reads |
| ---------- | ------------ |
| chr11      | 1,041,870    |
| chr1       | 945,718      |
| chr2       | 902,454      |
| chr5       | 865,417      |
| chr7       | 839,191      |
| chr4       | 767,100      |
| chr3       | 689,622      |
| chr9       | 683,826      |
| chr6       | 673,502      |
| chr17      | 662,542      |
| chr8       | 625,175      |
| chr13      | 571,932      |
| chr14      | 526,552      |
| chr10      | 525,797      |
| chr15      | 479,465      |

---

## Interpretation

This dataset can be described as:

> **A coordinate-sorted, single-end BAM file (~20 million reads) aligned with HISAT to a mouse-like reference genome, with ~65% of reads mapped and broad genome-wide coverage.**

It is suitable for:

- Coverage profiling
- Bin-based signal aggregation (bedGraph / BigWig)
- Benchmarking coverage tools
- UMI or mutation density analysis (if tags are present)

---

## Provenance

Alignment pipeline:

- Aligner: **hisat v0.1.5-beta**
- SAM/BAM processing: **samtools v1.19.2**

---

## Notes

- No paired-end information present.
- No duplicate-marked reads.
- No supplementary alignments.
- Secondary alignments are present (~3.9M), indicating multi-mapping reads.

---

*Generated automatically using `describe_bam.sh`.*
