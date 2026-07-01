# bam-subset-tag

`bam-subset-tag` splits a BAM file into multiple output BAM files based on BAM tag values.

Each supplied list file defines one subset output BAM.

The tool is particularly useful for:

- cell barcode extraction
- experimental partitioning
- targeted re-analysis
- debugging workflows
- and generation of smaller test datasets

---

# Typical Use Cases

- extracting selected cell populations
- splitting BAMs by barcode groups
- isolating experimental subsets
- debugging single-cell workflows
- generating reduced-size BAMs for testing

---

# Basic Usage

```bash
bam-subset-tag \
  --bam input.bam \
  --list selected_cells.txt \
  --prefix subset
```

---

# Multiple Output Groups

Multiple list files may be supplied.

Each list file creates one output BAM.

Example:

```bash
bam-subset-tag \
  --bam input.bam \
  --tag CB \
  --list tumor_cells.txt \
  --list immune_cells.txt \
  --prefix subsets/sample
```

This generates separate BAM files for:

- tumor cells
- immune cells

based on barcode membership.

---

# Input List Format

Each list file contains one accepted tag value per line.

Example:

```text
AAACCCAAGGTT
AAACCCATCGTA
AAACGGTTCGAT
```

---

# BAM Tags

The tool can match any two-character BAM tag.

Common examples:

| Tag | Meaning |
|---|---|
| `CB` | Corrected cell barcode |
| `CR` | Raw cell barcode |
| `UB` | Corrected UMI |
| `UR` | Raw UMI |

Default:

```text
CR
```

---

# Writing Unmatched Reads

Optional unmatched-read output:

```bash
--write-unmatched
```

This creates an additional BAM file containing reads that did not match any supplied tag list.

Useful for:

- QC
- debugging
- completeness checks
- and workflow validation

---

# Example

```bash
bam-subset-tag \
  --bam sc_reads.bam \
  --tag CB \
  --list selected_cells.txt \
  --prefix subsets/cells \
  --threads 8
```

---

# Output Structure

Outputs are written using the supplied prefix.

Example:

```text
subsets/cells_selected_cells.bam
```

and optionally:

```text
subsets/cells_unmatched.bam
```

---

# Performance Notes

The tool streams the BAM file directly and performs tag matching in memory.

Threading currently accelerates:

- BAM decompression
- BGZF reading
- BGZF writing

but does not parallelize downstream filtering logic.

---

# Typical Workflow Position

```text
Single-cell BAM
  ↓
bam-subset-tag
  ↓
subset BAMs
  ↓
targeted analysis or debugging
```
