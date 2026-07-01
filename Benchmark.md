# bam_tide Benchmark Summary

## Overview

We compared **bam_tide (Rust)** against **deepTools bamCoverage (Python)** using two datasets:

1. **GL000220.1** (very small chromosome fragment)
2. **pbmc_1k** (realistic medium-sized dataset)

For each dataset, multiple SAM flag exclusion scenarios were tested:

* 0
* 256
* 512
* 1024
* 2048

Metrics collected:

* Absolute differences between output BigWig tracks
* Pearson correlation
* Runtime
* Peak memory usage

---

## 1. GL000220.1 (Small reference)

### Key results

| Metric                  | Python   | Rust        |
| ----------------------- | -------- | ----------- |
| Runtime                 | ~46–50 s | ~4.2–4.4 s  |
| Peak RAM                | ~73 MB   | ~7.3–7.8 MB |
| Pearson correlation     | 1.0000   | 1.0000      |
| Max absolute difference | 0        | 0           |

### Interpretation

* **Rust is ~10× faster** than Python.
* **Rust uses ~10× less memory**.
* Outputs are **numerically identical** (correlation 1.0).

This indicates:

* The Rust implementation is **functionally equivalent**.
* Differences are negligible and likely due to:

  * Floating-point rounding

---

## 2. pbmc_1k (Realistic dataset)

### Key results (flag 0 representative)

| Metric                  | Python  | Rust   |
| ----------------------- | ------- | ------ |
| Runtime                 | 721.8 s | 91.3 s |
| Peak RAM                | 190 MB  | 232 MB |
| Pearson correlation     | 1.0     | 1.0    |
| Max absolute difference | 0       | 0      |

### Interpretation

**Performance**

* Rust is **~7–8× faster**.
* Memory usage:

  * Python: ~180–190 MB
  * Rust: ~227–232 MB
* Rust uses slightly more RAM here, likely due to:

  * Not depending on sorted bam files - collecting all data into memory

**Accuracy**

* Perfect


---

## Cross-dataset Summary

| Property                 | GL000220.1 | pbmc_1k        |
| ------------------------ | ---------- | -------------- |
| Speedup (Rust vs Python) | ~10×       | ~7–8×          |
| Memory (Rust vs Python)  | ~10× lower | ~20–30% higher |
| Correlation              | 1.0        | 1.0            |
| Fraction bins differing  | 0.0        | 0.0            |

---

## Main Conclusions

### 1. Speed

* Rust implementation is **consistently 7–10× faster**.
* This holds across small and realistic datasets.

### 2. Memory

* Small dataset: Rust is dramatically more efficient.
* Larger dataset: Rust uses slightly more RAM, but:

  * Still in the same order of magnitude
  * Likely tunable with streaming or chunked buffers

---

## Overall Assessment

**bam_tide** provides:

* Major speed improvements
* Comparable or better memory behavior
* Identical output to deepTools

This suggests it is a **valid high-performance replacement** for bamCoverage in many workflows.

---

