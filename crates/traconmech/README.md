# TraConMech

**Trace Context to Mechanism.**

TraConMech investigates whether biological observations are consistent with
candidate mechanisms represented by Lumrik's reference knowledge and external
evidence.

The first experiment is intentionally narrow: project expression-defined cell
populations onto Ommverse genomic architecture and test whether candidate
regions behave like differential transcriptional "light switches".

A missing experimental modality is not an error. TraConMech should continue
with the evidence that is available and report what could not be assessed.

## panc8 test data

`scripts/export_panc8_traconmech.R` exports the SeuratData `panc8` object to
`/tmp/panc8_traconmech` as a sparse Matrix Market matrix, ordered gene names,
and ordered cell metadata. R is used only to extract the existing dataset; the
experiment itself belongs in Rust.

## Investigation workspace: design and implementation status

TraConMech is growing into a scientific investigation workspace: compute expression
statistics, ask questions using attached evidence, receive both an explanation and
a structured proposal, inspect it, and approve changes before they affect the
working dataset. For example, a cell-identification question may propose new
`obs.cell_type` labels and a branching developmental graph. A later question
can propose an erythroid population for DBA/TP53 analysis and request additional
measurements. The original cluster identifiers and source H5AD remain intact.

Question templates are defined by **what they return**: an ordered column, a
new annotation, an analysis target, a gene set or a developmental graph. Each
result type maps to a named Rust callback; Ollama cannot run arbitrary code.
Changes to annotations and new analyses require validation and explicit user
approval. Every step records its evidence and provenance.

The planned web workspace has **Investigation**, **Data**, **Templates**,
**Evidence** and searchable **Help** tabs. Template, callback and evidence
provider documentation will come from one shared registry, so help can report
which functions work, which need configuration and which are still planned.
STRING, GSEA/MSigDB, ENCODE, InterPro and image support are architectural
targets, **not claims that those integrations are operational**.

**Implemented in this increment:** a typed registry with searchable metadata,
validated identifiers and machine-readable output:

```bash
cargo run -r -p traconmech -- --list-contracts
```

The interactive browser, file attachments, typed Ollama response validation,
callback execution and persistent approved annotations are **not implemented
yet**. See [AI_CONTRACT.md](AI_CONTRACT.md) for the normative design.

## Scientific methods handbook

The [TraConMech methods handbook](docs/README.md) documents the implemented
regulatory investigation, gene-model training and application, module scoring,
trajectories, Welch tests and Bonferroni correction. See the
[chapter index](docs/SUMMARY.md) for a book-style reading order.

The handbook distinguishes measured mean log1p module expression from
background-corrected AddModuleScore-style scores and describes important
limitations of cell-level inference.
