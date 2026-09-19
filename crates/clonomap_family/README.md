# ClonoMap

ClonoMap is Lumrik's receptor-family and mutation analysis crate. `clonomap_family`
reads reconstructed receptors from one or more `nelrune-vdj` output directories,
builds HC families using the same HC/CDR3 rules across the complete input pool,
validates mutation alignments, assigns light-chain clones, and writes cell-level
and family-level outputs.

## Joint multi-sample analysis

`clonomap_family` accepts one or more paths after `--vdj-out`:

```bash
clonomap_family \
  --vdj-out mouse1_bm/vdj_out mouse1_spleen/vdj_out \
  --out mouse1_clonomap \
  --reference-curator reference_curator.bin \
  --plots
```

There is deliberately no separate post-hoc merge algorithm. All receptors from
all inputs enter the normal ClonoMap family construction together. If receptors
from different inputs satisfy the ordinary HC/CDR3, HC/LC and mutation rules,
they share the same resulting family/clone. If they do not, they remain separate.

Each input receives a source label derived from its sample directory. Cell IDs
are internally scoped by that source so identical barcode strings from different
samples cannot collide. `cells.tsv` writes the source and original cell barcode
as separate columns.

When at least one family relationship is supported by more than one source,
ClonoMap additionally writes:

- `overlap_events.tsv` — machine-readable cross-source overlaps.
- `overlap_events.md` — compact human-readable summary.

The reported overlap levels are:

- `HC`: the same final HC family contains cells from multiple sources.
- `HC_LC`: the same HC family and LC clone contain cells from multiple sources.
- `HC_LC_MUTATION_SET`: cells from multiple sources share the same HC family,
  LC clone, and exact measured HC+LC mutation-coordinate set.

These files are derived summaries. `cells.tsv` remains the authoritative
cell-level output and is intended to carry the source/family/clone assignments
into later Seurat or Scanpy analysis.

A single input remains valid and behaves as before:

```bash
clonomap_family --vdj-out sample/vdj_out --out clonomap_out
```

## Reference curation

Reference-incompatible reconstructed sequence is handed to Lumrik's
`reference_curator`. ClonoMap decides biological model grouping; the curator
owns stable candidate identity, persistence, accumulated mapping evidence and
resolution state.

```bash
clonomap_family \
  --vdj-out sample/vdj_out \
  --out clonomap_out \
  --reference-curator reference_curator.bin
```

## Outputs

The main outputs are:

- `cells.tsv` — source-aware per-cell family, LC clone, mutation and plot metadata.
- `families.tsv` — final HC-family summaries and statistics.
- `reference_candidates.tsv` — reference-incompatible models seen by ClonoMap.
- `hc_unassigned.tsv` — receptors rejected from final HC membership.
- `plots/` — optional family SVGs with `--plots`.
- `overlap_events.tsv` / `overlap_events.md` — only when cross-source overlaps exist.

## Library

The crate also exposes the lower-level `ClonoMap`, family, PCA and MST types.
The geometric PCA/MST representation is descriptive: proximity is not by itself
proof of a parent/child phylogenetic relationship.

## Workspace development

```bash
cargo test -r -p clonomap
cargo build -r -p clonomap
```

## License

ClonoMap is part of Lumrik and uses the workspace AGPL-3.0-or-later license,
with Lumrik's commercial licensing option available separately.
