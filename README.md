# Lumrik

**High-performance single-cell sequencing tools in Rust**

Lumrik is a Rust workspace for processing, mapping, quantifying, and analysing single-cell sequencing data. It is built around a simple principle: large sequencing datasets should be processed as streams whenever possible, with compact representations and clearly separated components rather than repeatedly materialising large intermediate datasets.

The workspace contains reusable crates for primer and barcode detection, read normalisation, external mapper integration, gene and SNP quantification, sparse single-cell data storage, feature mapping, and guide assignment.

At the centre of the workspace is **nelrune**, Lumrik's end-to-end single-cell processing pipeline.

> **Status:** Lumrik is under active development. Interfaces and output formats may still change.

---

## nelrune

`nelrune` connects the Lumrik components into a complete sequencing workflow.

A typical Illumina run looks approximately like:

```
FASTQ R1/R2
     │
     ▼
sc_primer / bam_tide
├── chemistry / primer detection
├── cell barcode extraction
├── UMI extraction
└── molecule deduplication
     │
     ├──────────────► fast feature matching
     │                ├── sample tags
     │                ├── HTO
     │                └── CRISPR guides
     │
     ▼
streaming mapper
├── STAR
├── minimap2
└── BWA
     │
     ▼
BAM stream
     │
     ▼
bam_tide
├── exonic counts
├── intronic counts
└── optional SNP-aware quantification
     │
     ▼
sparse single-cell output
```

ONT/Dorado BAM input is also supported through the normalisation layer.

The important part is that these stages are designed to cooperate as a pipeline. Reads do not need to be converted into a succession of enormous intermediate representations before the next component can begin working.

---

## Find experimental sequences before a full analysis

`find-sequences` is a lightweight sanity/debugging tool for asking whether known
experimental sequences are present in sequencing data before committing to a
full analysis. Give it a FASTA containing constructs, reporters, vector sequence,
custom-capture targets, primer targets, spike-ins, or other sequences of interest.
It reports the matching FASTA entry and interval and, when the library structure
provides them, the cell barcode and UMI.

```bash
find-sequences \
    --fastq R1.fastq.gz \
    --r2-fastq R2.fastq.gz \
    --chemistry bd-v2-384 \
    --fasta constructs.fa \
    --out sequence_hits.tsv
```

Raw ONT BAM is supported directly; BAM records do **not** need to be mapped:

```bash
find-sequences \
    --bam raw_ont.bam \
    --fasta constructs.fa \
    --out sequence_hits.tsv
```

Long FASTA entries are tiled internally using Lumrik's fast feature matcher, but
the output is expressed in the original FASTA coordinates:

```text
cell_id  umi  feature  feature_start  feature_end  strand  matched_sequence
```

Overlapping internal tiles are collapsed before reporting. For long-read data the
final console summary also reports how many reads contain each FASTA feature and
the maximum number of separable occurrences observed in one read. This makes the
tool useful for early construct/primer validation as well as for inspecting
reporters, vectors, repeated inserts, and unexpected sequence content.

---

## Performance

Lumrik is designed for datasets containing tens of millions of sequencing reads while keeping memory consumption bounded wherever possible.

During development testing on **28 August 2026**, nelrune processed a real BD Rhapsody dataset using STAR with four mapper threads on a workstation with only **31 GiB of RAM**.

Observed during the run:

* more than **40 million reads** processed without memory growth proportional to read count;
* approximately **14,000–18,000 reads/second** end-to-end during the observed portions of the run;
* STAR itself occupied approximately **23 GiB resident memory** with the genome index loaded;
* nelrune continued normalisation, feature classification, mapping-result processing, and molecule accounting under the remaining memory constraint;
* the same development workload had previously exhausted the machine at approximately 8 million reads before removal of unnecessary per-read state retention.

These numbers are **development observations, not a controlled benchmark**. Performance depends strongly on sequencing chemistry, reference genome, mapper, storage, compression, CPU, feature references, and enabled analysis stages.

They do, however, demonstrate an important design property: processing additional reads does not inherently require retaining all previously processed read metadata in memory.

For production whole-genome STAR workloads, substantially more than 32 GiB RAM is recommended. The test above intentionally operated very close to the hardware limit.

---

## Implementation crates

Lumrik is a workspace, not a single monolithic implementation. The crate READMEs below are the implementation-level documentation and should be treated as the source of truth for their respective components.

| Crate | Responsibility |
|---|---|
| [`nelrune`](crates/nelrune/README.md) | End-to-end streaming integration of the Lumrik components |
| [`sc_primer`](crates/sc_primer/README.md) | Chemistry/read-structure grammars, cell/UMI extraction and primer detection |
| [`bam_tide`](crates/bam_tide/README.md) | High-throughput FASTQ/BAM normalization and quantification |
| [`sc-mapper`](crates/sc-mapper/README.md) | Streaming STAR/minimap2/BWA process integration |
| [`gtf_splice_index`](crates/gtf_splice_index/README.md) | Fast GTF/GFF parsing, spatial annotation and splice-aware indexing |
| [`snp-index`](crates/snp-index/README.md) | SNP indexing and allele-aware read matching |
| [`fast_tag_mapper`](crates/fast_tag_mapper/README.md) | Fast matching of sample tags, HTOs, guides and other short features |
| [`scdata`](crates/scdata/README.md) | Sparse UMI-aware single-cell count accumulation and export |
| [`sc-beacon`](crates/sc-beacon/README.md) | Ambient-aware cellular feature/guide calling |
| [`sc-vdj`](crates/sc-vdj/README.md) | V(D)J indexing, receptor reconstruction and AIRR-compatible output |
| [`valkyrn`](crates/valkyrn/README.md) | Repertoire interpretation, structural clone analysis and structure prioritization |
| [`clonomap`](crates/clonomap/README.md) | Mutation-aware PCA/MST geometry and plots for large receptor clones |
| [`sc-te`](crates/sc-te/README.md) | Single-cell transposable-element analysis and multimapper resolution |
| [`read-tag-table`](crates/read-tag-table/README.md) | External read/cell/UMI tag tables |
| [`mapping_info`](crates/mapping_info/README.md) | Shared counters, timings and reports |
| [`onehot_dna`](crates/onehot_dna/README.md) | Compact fixed-length DNA matching |
| [`int_to_dna`](crates/int_to_dna/README.md) | Compact sequence/integer identifier conversion |
| [`lumrik-status`](crates/lumrik-status/README.md) | Live run-status HTTP server/dashboard |

### Main capabilities

The individual crates are deliberately reusable, but several capabilities are especially central to Lumrik as a whole:

- **Fast genome annotation parsing.** `gtf_splice_index` streams GTF/GFF annotation into compact spatial and splice-aware indexes. The parser is not restricted to conventional gene/transcript records and is also used for transposable-element annotation.
- **Fast FASTQ/BAM handling.** `bam_tide` provides the high-throughput normalization and BAM-processing layer used for Illumina and ONT workflows, including cell/UMI handling, molecule deduplication, splice-aware quantification and optional SNP-aware processing.
- **Grammar-driven read structures.** `sc_primer` keeps sequencing chemistry out of downstream analysis code. A chemistry describes a read structure, and multiple structures can be requested for a single FASTQ stream.
- **Streaming mapper integration.** `sc-mapper` feeds accepted molecules directly to STAR, minimap2 or BWA and consumes mapper output while the run is still progressing.
- **Integrated single-cell processing.** `nelrune` composes these pieces into one run so filtering, feature routing, deduplication, mapping and quantification cooperate rather than repeatedly materialising the dataset.
- **Specialized downstream analysis.** `sc-vdj`, `sc-te`, `snp-index` and `sc-beacon` build on the same core representations instead of each reimplementing FASTQ/BAM infrastructure.

---

## Building

Lumrik currently targets modern Rust.

Clone the repository and build the complete workspace:

```
cargo build --release --workspace
```

To build only nelrune:

```
cargo build --release --bin nelrune
```

The resulting binary is located under:

```
target/release/nelrune
```

External mappers are not bundled. Install the mapper required for your workflow separately.

---

## Running nelrune

The exact command depends on sequencing chemistry and mapper configuration.

A simplified Illumina example is:

```
nelrune \
    --r1 sample_R1.fastq.gz \
    --r2 sample_R2.fastq.gz \
    --chemistry <CHEMISTRY> \
    --mapper star \
    --mapper-index /path/to/star_index \
    --mapper-threads 8 \
    --index reference.splice.idx \
    --threads 8 \
    --outpath nelrune_out
```

A simplified ONT example is:

```
nelrune \
    --bam dorado.bam \
    --chemistry <CHEMISTRY> \
    --mapper minimap2 \
    --mapper-index reference.mmi \
    --mapper-threads 8 \
    --index reference.splice.idx \
    --threads 8 \
    --outpath nelrune_out
```

Use:

```
nelrune --help
```

for the options supported by the current build.

---

## Live progress

nelrune includes a lightweight health/progress server for monitoring long-running analyses.

By default it uses port:

```
8787
```

Available endpoints include:

```
/          live dashboard
/health    health probe
/status    machine-readable status
```

This is particularly useful for long analyses running on workstations or compute nodes where the terminal itself does not need to remain the primary progress display.

---

## Design principles

### Stream first

Sequencing datasets are large enough already. Components should not create another complete copy of the dataset merely to pass information to the next stage.

### Keep ownership clear

Long-lived data belongs to the component responsible for it. Per-read jobs should contain read-specific state, not copies of global reference structures.

### Bound transient memory

Buffers should be large enough for efficient processing, but not so large that batching itself becomes the dominant memory consumer.

### Preserve information

Multimapping, splice structure, molecule identity, feature assignments, and allele observations should not be prematurely collapsed when later analysis may require them.

### Separate mechanisms

Primer detection, mapping, quantification, sparse storage, and statistical calling are separate crates. nelrune orchestrates them rather than reimplementing them.

### Measure on real data

Synthetic tests are essential for correctness, but performance decisions are also tested against real single-cell sequencing datasets.

---

## Development status

Lumrik is currently research software under active development.

The project already contains substantial working implementations for:

* Illumina single-cell normalisation;
* ONT/Dorado input;
* BD Rhapsody read structures;
* external streaming mapping;
* STAR, minimap2, and BWA integration;
* feature/sample-tag matching;
* exon/intron quantification;
* SNP-aware quantification;
* sparse single-cell data;
* CRISPR guide assignment;
* live progress reporting.

Additional sequencing systems, analysis modes, validation, documentation, and performance work are ongoing.

Do not assume that command-line interfaces or file formats are stable between development versions.

---

## Testing

Run the workspace tests with:

```
cargo test --workspace
```

For performance-sensitive code, release builds should always be used:

```
cargo test --release --workspace
```

Individual crates can be tested separately, for example:

```
cargo test -p sc_primer
cargo test -p sc-mapper
cargo test -p bam_tide
cargo test -p sc-beacon
```

Some integration tests require external tools or larger reference/test datasets and may therefore be ignored during a normal test run.

---

## Why Rust?

Single-cell sequencing pipelines combine several workloads that benefit from Rust:

* parsing very large compressed sequencing files;
* compact barcode and UMI representations;
* parallel sequence processing;
* sparse numerical data;
* high-throughput hashing;
* safe concurrency;
* streaming between processes;
* predictable ownership of large data structures.

Lumrik uses Rust not merely as a wrapper around existing bioinformatics tools, but as the implementation language for the data-intensive parts surrounding them.

External aligners such as STAR remain extremely specialised and highly optimised tools. Lumrik integrates them rather than attempting to replace them.

---

## License

Lumrik is developed by **Stefan Lang**.

The workspace is available under the **GNU Affero General Public License v3.0 or later (AGPL-3.0-or-later)**.

A commercial licensing option is also available. See `Commercial.md` for details.

Third-party components retain their respective licenses. See `THIRD_PARTY_NOTICES.md`.

---

## Project philosophy

Lumrik is intended to make sophisticated single-cell sequencing analysis possible without turning every new assay into another monolithic pipeline.

The long-term goal is a collection of efficient, interoperable building blocks where sequencing chemistry, mapping strategy, feature detection, quantification, and downstream analysis can evolve independently — while nelrune provides a practical route through them for complete datasets.

