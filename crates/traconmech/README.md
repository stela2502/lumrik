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
