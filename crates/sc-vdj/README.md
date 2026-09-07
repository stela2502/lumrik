# sc-vdj

Clean single-cell V(D)J reconstruction for Lumrik/Nelrune.

## Architecture

```text
BAM
  -> index/       exact V/D/J/C overlap + germline reference
  -> cellrep/     CellHash<u64 -> CellEvidence>, whole BAM feature + mapper evidence
  -> recombination/  sequence assembly, V/J anchoring, bounded D inference, P/N/deletion measurement
  -> output/      AIRR-compatible TSV
  -> runner.rs    glue only
```

The sequence and its mapper results stay in one `BamFeatureEvidence`; sc-vdj does not create a later biological identity from QNAME or UMI strings.

`Recombination` is the canonical biological result. Constant-region evidence is allowed in the observed receptor but deliberately excluded from the structural `HC:`/`LC:` identifier.

## Compact recombination IDs

Version 2 keeps the original nucleotide-independent structural identity semantics but encodes V/D/J as chain-local ordinals. The matching VDJ index determines how many hex digits each ordinal needs. Common junction measurements 0..14 remain one hex digit; larger values retain the old `Fxxxx` escape. `vdj-decode` also accepts version-1 IDs using the previous three-hex-digit global segment indices.

## Tests

The important regression is `tests/test-complex-vdj-detection.rs`. It builds a synthetic IGH receptor with V/D/J coding-end deletions, V/J P additions, N1/N2, a short D, and Igha continuation. Four overlapping BAM fragments map independently to V, D, J and C; no BAM record contains the whole recombination. The test requires sc-vdj to reconstruct the receptor, recover the exact V/D/J/C calls and junction measurements, and round-trip the compact structural ID.

Run in release mode:

```bash
cargo check --release -p sc-vdj
cargo test --release -p sc-vdj
cargo test --release -p sc-vdj --test test-complex-vdj-detection -- --nocapture
```

The new index format is `LVDJIDX5`; regenerate existing `.vdjidx` files with `vdj-index` before using them with this crate.
