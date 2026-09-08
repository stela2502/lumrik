# sc-vdj

Single-cell V(D)J reconstruction for Lumrik/Nelrune.

`sc-vdj` reconstructs immunoglobulin and T-cell receptor rearrangements directly from single-cell BAM data. It is designed for transcriptomic single-cell datasets where receptor evidence may be fragmented across multiple reads and where no individual BAM record necessarily contains a complete V(D)J rearrangement.

The implementation is built around two requirements:

1. receptor reconstruction must combine evidence across fragments without treating a read name or UMI as the biological receptor itself;
2. memory use must depend primarily on the compact evidence retained for relevant cells, not on the number of BAM records processed.

The current implementation can use an expression-derived cell set from Nelrune as an early processing gate. A future V(D)J-native cell caller is intended to make expression-derived cell selection optional even for large BAMs.

## Current processing model

```text
VDJ reference
  |
  +-> VDJ index
       V / D / J / C segments
       germline sequence
       lookup structures

single-cell BAM
  |
  +-> optional preliminary cell gate
  |     |
  |     +-> --exonic supplied:
  |     |     use expression-called barcodes as allowed cells
  |     |
  |     +-> --exonic omitted:
  |           collect receptor evidence directly from BAM barcodes
  |
  +-> first BAM pass
  |     |
  |     +-> identify V/D/J/C-overlapping evidence
  |     +-> accumulate bounded raw-evidence batches
  |     +-> split batch by cell
  |     +-> compact cells in parallel
  |     +-> merge compact deltas into persistent CellHash
  |     +-> discard raw batch
  |
  +-> receptor reconstruction
  |     |
  |     +-> build cell-specific receptor sequence evidence
  |     +-> V/J anchoring
  |     +-> bounded D inference
  |     +-> junction reconstruction
  |     +-> identify candidate recombinations
  |
  +-> second BAM pass
  |     |
  |     +-> rescan wanted cells in bounded batches
  |     +-> rediscover receptor/CDR3 sequence
  |     +-> confirm constant-region evidence
  |     +-> retain compact fragment linkage
  |
  +-> output
        AIRR-compatible rearrangements
        receptor sequences
        cell/sample summaries
        development/debug summaries
```

## Evidence model

A BAM record is evidence for a receptor; it is not itself a receptor.

`BamFeatureEvidence` represents transient BAM-derived evidence, including sequence and mapping information required to update the receptor model. These objects are processed in bounded batches and are not retained for the lifetime of the analysis.

Within each batch, evidence is grouped by cell and converted into compact cell-local summaries. Cell-local compaction is parallelized, after which the compact deltas are merged into the persistent:

```text
CellHash<u64 -> CellEvidence>
```

Raw batch evidence is then discarded.

Persistent `CellEvidence` contains the information needed to reconstruct receptor sequence and segment support without retaining a copy of every contributing BAM record. Sequence evidence is represented through compact per-position observations rather than an ever-growing collection of raw reads.

This distinction is important for large single-cell BAMs: memory consumption should follow the amount of persistent receptor evidence rather than scale directly with BAM size.

## Cell selection

### Expression-derived preliminary cells

When `--exonic` is supplied, `nelrune-vdj` reads the barcode set from the exonic expression matrix and uses it as the preliminary allowed-cell population.

For a normal Nelrune result:

```bash
nelrune-vdj \
    --exonic <nelrune-output>/exonic \
    --bam <nelrune-output>/nelrune.mapper.bam \
    --index reference.vdjidx \
    --out vdj-output
```

`--exonic` may point to the exonic matrix directory or to a supported Nelrune output location containing it.

The barcode set is used as an **early BAM-ingestion gate**. Once a BAM record has been assigned to a cell barcode, records outside the allowed population are rejected before expensive V(D)J evidence processing and before persistent receptor state is allocated for that cell.

The expression-derived set is therefore not merely metadata or a reporting filter. It bounds the cell population entering receptor reconstruction.

The exonic barcode set does **not** pre-populate the receptor `CellHash`. A cell enters persistent V(D)J state only when receptor evidence for that cell is actually observed.

### BAM-derived mode

If `--exonic` is omitted, `sc-vdj` does not require expression support and can collect receptor evidence directly from the barcode population observed in the BAM.

This behavior is intentional. V(D)J evidence can be biologically useful even when a cell is absent from, or poorly represented in, an expression-derived cell call.

For very large single-cell BAMs, unrestricted BAM-derived processing can currently admit a very large background barcode population. The expression-derived gate is therefore the practical high-throughput mode when a suitable expression matrix is available.

## Planned V(D)J-native cell detection

The intended next step is an internal V(D)J cell-calling stage.

Instead of treating every barcode with any receptor overlap as a receptor-bearing cell, `sc-vdj` should use the receptor evidence itself to distinguish plausible cells from the much larger population of low-level/background barcode observations.

Conceptually:

```text
BAM barcode population
        |
        v
cheap receptor-evidence collection
        |
        v
VDJ-native cell scoring / calling
        |
        +-----------------------+
        |                       |
        v                       v
likely receptor cells       background
        |
        v
full receptor reconstruction
```

This should allow three useful operating modes:

```text
1. Expression-assisted
   --exonic supplied
   -> expression-called cells provide the preliminary gate

2. VDJ-native
   no expression matrix
   -> receptor evidence identifies the cells worth reconstructing

3. Combined
   expression calls + VDJ evidence
   -> retain expression-supported cells while allowing strong
      receptor-only evidence to recover additional cells
```

The important design principle is that expression support should remain useful without becoming a biological requirement for calling a V(D)J receptor.

The current no-`--exonic` collect-all behavior is retained partly so this V(D)J-native path can be developed without changing the underlying BAM evidence model.

## Receptor reconstruction

Receptor reconstruction operates on compact cell evidence rather than requiring a complete rearrangement to occur in one BAM record.

Evidence from multiple fragments can therefore contribute independently to:

* V support,
* D support,
* J support,
* constant-region support,
* receptor consensus sequence,
* junction structure,
* and final recombination support.

This is particularly important for single-cell RNA sequencing, where coverage of a receptor transcript can be fragmented and individual reads may cover only one part of the rearrangement.

The reconstruction logic performs V/J anchoring and bounded D inference and measures junction structure including deletions and P/N nucleotide contributions where supported by the reconstructed sequence.

`Recombination` is the canonical reconstructed biological result.

Constant-region evidence describes the observed receptor and can be used for receptor annotation and confirmation, but constant-region identity is deliberately excluded from the compact structural `HC:`/`LC:` recombination identifier.

## Second-pass confirmation

After candidate receptors have been reconstructed, `nelrune-vdj` performs a second BAM pass.

This pass is also bounded rather than accumulating all matching BAM records in a BAM-sized map.

It uses the reconstructed cell-specific receptor information to search relevant records for:

* receptor/CDR3 rediscovery,
* constant-region evidence,
* fragment linkage,
* and constant calls that can be confirmed or rescued from the additional evidence.

The second pass uses the same general memory principle as initial evidence collection:

```text
bounded BAM batch
    -> cell-local processing
    -> compact persistent result
    -> discard raw batch
```

## Parallelism and memory model

The expensive cell-local parts of evidence compaction and receptor reconstruction are parallelized with Rayon.

The intended scaling model is:

```text
raw BAM evidence:
    bounded

temporary batch memory:
    bounded by batch size

persistent memory:
    proportional to compact evidence for retained cells

cell-local computation:
    parallel
```

This avoids retaining all receptor-overlapping BAM records or all `(cell, QNAME)` combinations for the duration of the analysis.

For large datasets, restricting the preliminary population with `--exonic` currently has a much larger performance effect than simply adding threads, because it prevents background barcode IDs from reaching the expensive reconstruction stage at all.

## Live run status

`nelrune-vdj` provides a lightweight HTTP status server for monitoring long-running analyses.

By default, the dashboard is available on port `8787`:

```text
http://localhost:8787/
```

The machine-readable status endpoint is:

```text
http://localhost:8787/status
```

A compact terminal representation can be obtained with `jq`:

```bash
curl -s http://localhost:8787/status | jq -r '
.stage,
(.sections[] | "\n\(.title)", (.metrics[] | "  \(.label): \(.value)"))
'
```

This reports the current stage together with live metrics for:

* initial evidence collection,
* receptor reconstruction,
* CDR3/constant-region confirmation,
* stage and total execution time,
* and memory consumption.

A simple health endpoint is available at:

```bash
curl http://localhost:8787/health
```

Relevant command-line controls are:

```text
--health-port <PORT>
--health-hostname <HOSTNAME>
--no-health-server
```

For a run on a cluster compute node, the dashboard can be forwarded to the local machine with SSH:

```bash
ssh -N -L 8787:<compute-node>:8787 <cluster-login>
```

and then opened locally at:

```text
http://localhost:8787/
```

The status dashboard contains aggregate run statistics only; it is not intended as a per-cell receptor browser.

At successful completion, the same final run state is also written as:

```text
vdj-run-summary.yaml
vdj-report.html
```

The YAML file preserves numeric run counters for regression and workflow use. The HTML file is a portable static snapshot of the final dashboard and does not require the health server to remain running.

## V(D)J reference index

`sc-vdj` uses a dedicated V(D)J index containing the germline segments required for reconstruction.

An index can be generated from a genome and annotation with `vdj-index`:

```bash
vdj-index \
    --gtf annotation.gtf \
    --genome genome.fa \
    --out reference.vdjidx
```

The current index format is:

```text
LVDJIDX6
```

Existing indices written in older incompatible formats should be regenerated before use.

The index provides the V, D, J and constant-region germline sequence information used by both initial evidence detection and later receptor reconstruction.

## Compact recombination IDs

The compact recombination identifier represents the structural V(D)J rearrangement independently of its nucleotide sequence representation.

Version 2 encodes V/D/J segments as chain-local ordinals. The associated V(D)J index determines the number of hexadecimal digits required for each ordinal.

Common junction measurements in the range `0..14` use a single hexadecimal digit; larger values use the extended `Fxxxx` representation.

`vdj-decode` can decode the compact representation and also accepts version-1 identifiers using the earlier three-hex-digit global segment indices.

Constant-region identity is intentionally not part of the structural `HC:`/`LC:` identifier.

## Output

The primary biological output is the reconstructed `Recombination`.

`nelrune-vdj` additionally writes machine-readable and human-readable summaries describing the detected receptor population, including AIRR-compatible rearrangement output and sequence-oriented outputs used for inspection and downstream analysis.

Depending on the run and enabled outputs, these include information such as:

```text
vdj_rearrangements.tsv
vdj_rearrangement_candidates.fasta
vdj_rearrangement_candidate_proteins.fasta
vdj_cell_summary.tsv
vdj_sample_summary.tsv
vdj-mapping-info.txt
vdj-run-summary.yaml
vdj-report.html
```

The AIRR-oriented output describes reconstructed receptor calls rather than individual BAM alignments.

## Tests

The central synthetic regression is:

```text
tests/test-complex-vdj-detection.rs
```

It constructs an IGH receptor containing V/D/J coding-end deletions, V/J P additions, N1/N2 sequence, a short D segment, and constant-region continuation.

The receptor is deliberately fragmented across multiple BAM records: individual fragments map independently to V, D, J and C, and no BAM record contains the complete rearrangement.

The test therefore exercises the central assumption of `sc-vdj`: a receptor must be reconstructed from distributed single-cell evidence rather than requiring one read to describe the entire biological rearrangement.

The regression checks reconstruction of the receptor, V/D/J/C assignment, junction measurements, and compact structural-ID round-tripping.

Run the crate tests in release mode:

```bash
cargo check --release -p sc-vdj
cargo test --release -p sc-vdj
cargo test --release -p sc-vdj --test test-complex-vdj-detection -- --nocapture
```

## Development direction

The current architecture establishes the main bounded-memory reconstruction pipeline:

```text
bounded BAM evidence
    -> compact per-cell state
    -> parallel receptor reconstruction
    -> bounded confirmation pass
```

The next major architectural step is V(D)J-native cell detection.

That work should make the distinction between:

```text
barcode observed in BAM
```

and:

```text
cell with convincing receptor evidence
```

explicit before full receptor reconstruction.

Once that exists, `--exonic` can become one source of cell evidence rather than the only practical way of bounding a very large barcode population.

The long-term goal is therefore not to require expression-based cell calls, but to combine complementary evidence:

```text
expression evidence
        +
receptor-specific evidence
        |
        v
candidate single cells
        |
        v
V(D)J reconstruction
```

This keeps `sc-vdj` useful both as part of the complete Nelrune single-cell workflow and as a standalone receptor-analysis component for BAMs where expression-derived cell calls are unavailable or incomplete.

