# sc-vdj Performance Notes

This document records observed performance of `sc-vdj` at selected Git revisions. These measurements are intended as engineering baselines rather than formal benchmarks.

## Revision `55076ba` — Batched parallel evidence processing

**Git commit:** `55076bae2ecb29e3f6fc5c251f0098deb165e424`

This revision introduced bounded, parallel processing of V(D)J evidence using `CellHash`.

Raw BAM evidence is collected in batches of **200,000 records**, partitioned by cell, compacted in parallel using Rayon, and merged into persistent per-cell receptor summaries. Raw evidence is discarded after each batch.

This solved the previously observed unbounded memory growth during BAM ingestion and restored effective multicore utilization.

### Test dataset

The run used the large NELRUNE V(D)J dataset that previously caused memory exhaustion during receptor evidence collection.

At completion of the initial evidence-collection stage:

| Metric                                | Observed value |
| ------------------------------------- | -------------: |
| Receptor-overlap records              |     16,831,205 |
| BAM cell IDs with receptor evidence   |      1,027,661 |
| Persistent compact receptor summaries |      1,819,343 |
| Receptor-evidence fragments           |     16,776,603 |
| Reference V(D)J segments              |            811 |
| Worker threads                        |              8 |

### Resource usage

During batched evidence collection, CPU utilization commonly reached approximately **500–800%** on 8 worker threads.

A representative observation during the later receptor-reconstruction stage was:

| Metric                 |                         Observed value |
| ---------------------- | -------------------------------------: |
| CPU utilization        |                                  ~794% |
| Resident memory        |                               ~6.8 GiB |
| Peak resident memory   |                               ~6.8 GiB |
| Runtime at observation |                        ~4 days 3 hours |
| Pipeline stage         | Reconstructing cell-specific receptors |

Importantly, the run remained alive and computationally active. The earlier catastrophic linear accumulation of raw BAM evidence had therefore been eliminated.

### Performance result

**The run did not complete in a practically acceptable time.**

The bottleneck had moved from evidence ingestion to receptor reconstruction.

The initial BAM scan accepted receptor evidence from **1,027,661 distinct BAM cell IDs** and generated **1,819,343 persistent receptor summaries**. Reconstruction was subsequently attempted across this enormous candidate population.

After approximately **4 days and 3 hours**, the process was still in:

> **Stage 2/3 — reconstructing cell-specific receptors**

The subsequent CDR3/constant-region BAM rescan had not yet started.

### Interpretation

Revision `55076ba` demonstrates that the batched `CellHash` architecture successfully addresses the original memory-scaling problem:

**Before batching**

```text
BAM records
    ↓
retain raw receptor evidence
    ↓
memory grows approximately with the number of routed reads
    ↓
OOM / impractical memory consumption
```

**At revision `55076ba`**

```text
BAM records
    ↓
200,000-record batch
    ↓
parallel per-cell compaction
    ↓
persistent CellHash receptor summaries
    ↓
discard raw batch
```

This reduced transient evidence memory dramatically while allowing effective parallel execution.

However, the benchmark exposed a second scaling problem: **cell selection**.

The BAM contained receptor-associated evidence for more than one million barcode-derived cell IDs. Most of these should not automatically become expensive receptor-reconstruction candidates when an expression-derived preliminary cell set is available.

Consequently, this revision should be considered an important intermediate performance baseline:

* **Raw evidence memory scaling:** solved.
* **Parallel evidence compaction:** effective.
* **CPU utilization:** good.
* **Persistent summary population:** too large.
* **Receptor reconstruction runtime:** unacceptable for an unfiltered BAM barcode population.

## Next optimization target

When an expression-derived cell set is supplied through `--exonic`, those barcodes can provide an early preliminary cell population.

The expected processing model is:

```text
Expression-derived preliminary cells
                ↓
          accepted Cell IDs
                ↓
BAM ──→ reject unrelated barcodes immediately
                ↓
       V(D)J evidence collection
                ↓
       receptor reconstruction
```

This should reduce both the persistent receptor-summary population and the reconstruction workload by orders of magnitude for datasets where expression-backed cell calls are available.

A separate V(D)J-native cell-calling strategy remains desirable for datasets where no expression-derived cell set exists, or where biologically useful receptor-bearing cells may not be represented in the expression-derived population.

Future performance measurements should be recorded against this revision to quantify the effect of cell selection independently from the batching and parallelization improvements established here.

