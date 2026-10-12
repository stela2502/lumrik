# TraConMech AI contract

## Purpose and boundaries

TraConMech is an evidence-driven scientific investigation workspace, not an autonomous source of biological facts. `sc-analysis` produces experimental measurements; Ommverse resolves biological identity and reference knowledge; Ollama interprets evidence and proposes next steps. Missing modalities must remain explicitly missing. Observations, reference associations, hypotheses and user-approved conclusions must never be conflated.

## Investigation lifecycle

1. Load a dataset without modifying the original; generate and register reproducible analysis artifacts (initially marker summaries and top-20 genes per population).
2. The investigator selects or edits a question template, attaches explicit evidence files, and submits the question. A template specifies required inputs, expected typed outputs and proposed callbacks.
3. Ollama returns both a human-readable explanation and machine-readable proposals. The response is *untrusted*: validate IDs, schemas, referenced evidence, permissible actions and output size before offering any action.
4. Display proposed changes with Accept / Edit / Reject. Only approved callbacks may mutate the *working* annotation state, initiate new analyses or export an annotated dataset. Preserve original Louvain labels; developmental ordering changes category order or creates new metadata, never relabels cluster IDs.
5. Record question, template version, evidence references, model, structured proposal, approval decision, callback outcome and generated artifacts. A follow-up question can reuse accepted outputs, but not unapproved suggestions as established facts.

## Typed results and callback contracts

A template is defined by its **return value**, not only its prompt. Initial result types:

- `column_order`: proposed ordered categorical levels for an existing `obs` column; validate each level against actual observed labels. May represent one lineage path; do not silently linearize branching trajectories.
- `annotation_column`: proposed new per-cell or per-population `obs` annotation; validate coverage, type and identifiers. Apply to working state only after approval.
- `analysis_target`: proposed population, condition, gene or comparison for an analysis; validate selections and launch a named registered analysis provider after approval.
- `gene_set`: proposed named genes and provenance; resolve identifiers through Ommverse and mark missing or ambiguous names.
- `developmental_graph`: proposed directed edges between actual population IDs, with evidence references and uncertainty. A graph may branch; do not equate inferred trajectories with lineage tracing.

A single question may return multiple typed proposals. Callback names refer to **registered Rust implementations**, never executable model-generated code or arbitrary shell commands. Validation and human approval are mandatory before writing data or starting consequential jobs. A callback must report success/failure and resulting artifacts; failures must not corrupt earlier accepted state.

## Evidence provider contract

Experimental providers calculate measurements from the user's dataset; reference providers return contextual knowledge. Each provider declares its identifier, description, accepted query schema, output schema, reference/index prerequisites, provenance, limitations and availability (`available`, `not_configured`, `planned`). Candidate sources include expression, SNP/variants, VDJ, genome/transcript models, protein domains (InterPro), STRING associations, GSEA/MSigDB gene sets, ENCODE regulation and images. **Listing a provider is not evidence that it is implemented.** Never invent measurements, regulatory edges, protein interactions or images.

Every generated evidence object records producer, method, input identity, parameters, timestamp/version and paths to underlying files. Keep full tables on disk and send selected, bounded excerpts to the model. Pseudo-bulk groups of cells are not independent biological replicates; patient/control inference must respect sample-level replication.

## Web interface and documentation registry

The web workspace is organized into Investigation (question/answer boxes and attachments), Data (artifacts and annotations), Templates (search/edit/use), Evidence (providers and availability), and Help (searchable reference and examples). The **same versioned registry** drives the template selector, model-facing tool descriptions, callback validation, contextual help and search. A template/provider/callback ships with its own human-readable documentation, inputs, outputs, limitations and example. Documentation must reflect implemented capability, not plans.

File attachments are explicitly selected for each question, limited to investigation-approved paths, previewed and size-bounded. A selected file is not automatically an instruction. Never allow model-generated paths to grant file access. Maintain conversation and evidence history across restarts; stream long-running model progress without blocking the web server. The existing `lumrik-status` server is monitoring infrastructure, not yet a conversation API; extend/reuse it deliberately rather than claiming the web workspace already exists.

## First implementation increment

Provide a typed, machine-readable registry of templates, result types, callbacks and evidence providers; validate the registry; expose a CLI command to print it as JSON and include explicit capability states. This increment **does not** yet implement the browser, Ollama structured response parsing, approvals, annotation mutation or reference integrations. Those must follow as separate tested increments.
