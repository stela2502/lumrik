# bw-compare

`bw-compare` compares two BigWig coverage tracks and reports per-chromosome and global differences.

The tool was originally designed to validate the Rust implementation of `bam-coverage` against existing Python-based workflows such as `deeptools bamCoverage`.

Both BigWig files are binned using the same bin width and compared position-by-position.

The resulting report summarizes:

- absolute signal differences
- variance
- RMSE
- maximum deviation
- and Pearson correlation

between the two tracks.

---

# Typical Use Cases

- validating `bam-coverage` output
- comparing coverage generation methods
- regression testing after algorithm changes
- QC of new normalization or filtering strategies
- benchmarking replacement workflows

---

# Basic Usage

```bash
bw-compare \
  --python-bw python.bw \
  --rust-bw rust.bw
```

This compares both tracks using the default bin width of 50 bp.

---

# Example With Explicit Parameters

```bash
bw-compare \
  --python-bw sample_python.bw \
  --rust-bw sample_rust.bw \
  --bin-width 50 \
  --eps 0.00001 \
  --outfile comparison.txt
```

---

# Output Statistics

The report contains one line per chromosome and a final `TOTAL` summary line.

Reported columns:

| Column | Description |
|---|---|
| `n_over_eps` | Number of bins where `|python - rust| > eps` |
| `frac_n_over_eps` | Fraction of bins exceeding epsilon |
| `mean_abs` | Mean absolute difference |
| `var_abs` | Variance of absolute differences |
| `rmse` | Root mean squared error |
| `max_abs` | Maximum absolute difference |
| `pearson_rho` | Pearson correlation between tracks |

---

# Important Notes

## Matching bin width

Both BigWig files must have been generated using the same bin width.

Example:

```bash
bam-coverage --bin-width 50
```

must be compared against another track generated with:

```bash
--bin-width 50
```

otherwise differences are expected.

---

## Interpretation

A successful comparison usually shows:

- very small `mean_abs`
- very small `rmse`
- high `pearson_rho`
- and very few bins above epsilon

Small floating point differences are expected between implementations.

---

# Command Line Options

| Option | Description |
|---|---|
| `--python-bw` | Reference BigWig file |
| `--rust-bw` | BigWig file to validate |
| `--bin-width` | Bin width used during coverage generation |
| `--eps` | Difference threshold used for counting bins as different |
| `--outfile` | Optional report output file |

---

# Example Workflow

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
