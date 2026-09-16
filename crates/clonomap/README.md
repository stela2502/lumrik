# ClonoMap

ClonoMap is Lumrik's mutation-aware geometric analysis crate for large receptor clones. It projects aligned receptor sequences into a consensus-relative mutation space, reduces that space with PCA, and builds a minimum-spanning tree (MST) over the resulting coordinates.

ClonoMap does **not** decide whether two receptors are members of the same biological clone. Clone identity belongs upstream. In the Lumrik V(D)J workflow, `sc-vdj` reconstructs receptors and structural recombination IDs, and Valkyrn defines the heavy-chain family context. ClonoMap then asks how sequence states are arranged *within* a supplied clone.

## Role in Lumrik

The current flow is:

```text
sc-vdj / nelrune-vdj
        |
        v
reconstructed receptors + HC:/LC: structural IDs
        |
        v
Valkyrn
  - conservative HC/LC family definitions
  - mutation classification
  - isotype composition
  - selection of large HC families
        |
        v
ClonoMap
  - consensus-relative sequence encoding
  - PCA
  - MST
  - large-clone geometry / plots
```

Valkyrn currently uses ClonoMap as an exploratory view for large IGH families. Before passing sequences to ClonoMap, each observed receptor is aligned back onto its reconstructed naive-recombination coordinate system. This gives ClonoMap equal-width sequences without asking it to redefine the clone.

The resulting PCA/MST is a visualization of sequence-state geometry. It is **not automatically a phylogeny**, and the sparsest PCA region is not automatically a proven biological ancestor.

## Library

The main API is `ClonoMap`:

```rust
use clonomap::ClonoMap;

let sequences = vec![
    "ACGTACGT".to_string(),
    "ACGTTCGT".to_string(),
    "ACGTACGA".to_string(),
];

let model = ClonoMap::new(sequences, 3, false)?;
println!("{} PCA rows", model.coords().nrows());
println!("{} MST edges", model.tree().len());
```

`OneHotEncoder` stores the supplied aligned sequences and encodes each position relative to the clone consensus. `PcaModel` performs dimensionality reduction and `MstTree` constructs the geometric tree.

## Binaries

The crate retains standalone tools that are useful outside Valkyrn:

- `clonomap` analyses one equal-width sequence set.
- `clonomap_batch` reads Change-O-style clone tables and analyses sufficiently large clones.
- `tree_viewer` renders Newick trees when the `plot` feature is enabled.

For new Lumrik V(D)J analyses, prefer Valkyrn as the entry point because it uses Lumrik's structural recombination identity rather than importing an external `clone_id`.

## Plotting

Plotting is optional at the ClonoMap crate level:

```bash
cargo build --release -p clonomap --features plot
```

Valkyrn enables the plotting feature because its current large-family test writes PCA and MST PNGs under:

```text
valkyrn/clonomap/<HC-family>/
├── coords.tsv
├── tree.tsv
├── rows.tsv
├── pca.png
└── mst.png
```

`valkyrn/clonomap/clonomap_summary.tsv` records which large families were analysed or skipped.

## Current interpretation boundary

ClonoMap is intentionally descriptive. A short MST edge means two supplied sequence states are close in the PCA representation; it does not by itself prove a parent/child relationship. Likewise, PCA clusters can reveal structure worth inspecting without establishing a formal clonal phylogeny.

The immediate use in Valkyrn is therefore pragmatic: find the large structural HC families, visualize their mutation geometry, and compare those patterns with light-chain partners, isotype states, and Valkyrn's mutation annotations.

## Workspace development

ClonoMap is a normal Lumrik workspace crate and inherits the workspace version, Rust edition, and AGPL-3.0-or-later license.

Build and test in release mode:

```bash
cargo check --release -p clonomap
cargo test --release -p clonomap
```

To validate the Valkyrn integration:

```bash
cargo check --release -p valkyrn
cargo test --release -p valkyrn
```

## License

ClonoMap is part of Lumrik and is licensed under the workspace license: **AGPL-3.0-or-later**, with Lumrik's commercial licensing option available separately.
