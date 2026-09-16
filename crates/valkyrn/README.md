# Valkyrn

> **Valkyrn turns reconstructed immune receptors into biological hypotheses.**

Valkyrn is Lumrik's repertoire-interpretation layer. `sc-vdj`/`nelrune-vdj` should remain the rigorous reconstruction engine: they decide which receptor sequences are supported by the sequencing evidence. Valkyrn starts *after* that point and asks what the reconstructed repertoire says biologically.

The first implementation deliberately focuses on questions that can be answered from Lumrik's own evidence-aware output rather than importing a second clonotyping model.

## What Valkyrn does now

Given a `nelrune-vdj` output directory, Valkyrn joins `vdj_calls.tsv` to `airr_rearrangements.tsv` through `lumrik_recombination_id` and performs several linked analyses.

### 1. Look hard at the cells

`valkyrn_cells.tsv` counts productive and total heavy/light rearrangements per cell and records whether a cell has a productive heavy/light pair. This keeps secondary and unproductive rearrangements visible instead of silently forcing one receptor per cell.

### 2. Find receptor families

Valkyrn does not invent a second clonotyping definition. `sc-vdj` already emits compact reversible structural recombination IDs. `HC:<HEX>` encodes chain-local V/D/J identity plus V/D/J trimming, retained D length, P-addition lengths, N1/N2 lengths and the P/N ambiguity flag. Constant-region identity and nucleotide substitutions are intentionally excluded. That makes the existing HC ID the primary heavy-chain clone/lineage key and allows somatic sequence changes to accumulate without changing the underlying recombination identity.

P/N decomposition is not always unique. `sc-vdj` marks such calls with `pn_alternative`. Valkyrn keeps exact `HC:<HEX>` identity as the normal rule, but ambiguous calls may be connected across different compact IDs only when V, D and J calls, CDR3 amino-acid length, and reconstructed naive-rearrangement length all agree. Such a family is marked with `~PNALT` in the family name; the original recombination IDs remain present in the detailed output. This fallback is intentionally much narrower than the old V/J + CDR3-distance clustering.

Light chains retain the stricter trust boundary. `LC:<HEX>` is a structural light-chain rearrangement identity, but LC identity alone is not treated as proof of clonality across cells. A multi-cell light-chain clone is represented as the paired state `HC:<HEX>+LC:<HEX>` (or the HC `~PNALT` family plus LC ID). The same LC recombination on different HC backgrounds is recurrence, not one clone.

### 3. Inspect the underlying mutation pattern

Lumrik already reconstructs both `naive_recombination` and an error-corrected `observed_receptor_sequence`. Valkyrn compares those sequences base by base.

`valkyrn_mutations.tsv` keeps the labels intentionally modest. A difference from a one-cell family is `singleton_observed_difference`. In multi-cell families, a difference found in every member is `family_shared_candidate`; one found in only part of the family is `branch_or_private_candidate`.

Valkyrn also counts whether the same aligned difference recurs in independent receptor families with the same chain and V call. Such events are labelled `recurrent_same_v_candidate` and the number of independent families is reported. This is useful evidence for a possible reference/germline-allele mismatch, but it is **not** called a germline allele yet: the current comparison is in reconstructed-receptor coordinates, not segment-aware germline coordinates.

These are descriptive categories, not final SHM calls. Once a common ancestral rearrangement is supported, changes inside the junction/CDR3 can be useful lineage information too; they should not be discarded merely because the original junction was created by V(D)J recombination.

### 4. Inspect isotype/class-switch state

For each receptor family Valkyrn retains the constant-region call and reports the observed isotype composition. Because the structural `HC:<HEX>` identity intentionally excludes the constant region, the same VDJ rearrangement can be followed across IgM/IgD/IgG/IgA/IgE states without manufacturing separate HC clones.

This is descriptive lineage evidence, not by itself a reconstructed class-switch tree.

### 5. Plot large HC families with ClonoMap

Large IGH families are handed to the workspace `clonomap` crate for an exploratory mutation-geometry view. Valkyrn aligns each observed HC back onto the naive-recombination coordinate system, then ClonoMap performs consensus-relative encoding, PCA and MST construction. The default threshold is 100 cells and can be changed with `--min-clonomap-family`.

The PCA/MST plots are deliberately not called phylogenies. They show sequence-state geometry inside an already-defined structural HC family and are intended to expose large branches and unusual mutation patterns for further analysis.

### 6. Prepare structural questions instead of blindly folding everything

Valkyrn does **not** make an antibody-structure package a core Rust dependency. Instead it creates a stable hand-off under `structure_candidates/`:

- `manifest.tsv` explains why each HC/LC pair was selected.
- `paired_receptors.fasta` contains the error-corrected observed heavy and light sequences.

The first prioritization rule selects HC/LC pairs from expanded IGH families that occur with more than one productive light-chain family. These are exactly the cases where structure can test a biological question: does a related heavy-chain state preserve a similar recognition surface while tolerating different light chains, or does each HC/LC combination create a different solution?

This hand-off is designed so Norn can later add an **optional** structure process using ABodyBuilder3 (or another antibody model) without coupling `valkyrn` to Python, CUDA, model weights, or a particular GPU environment. The structure runner should consume the manifest/FASTA contract and write model provenance, confidence and PDB/mmCIF outputs beside the candidate IDs.

## Usage

```bash
cargo run --release -p valkyrn -- \
    --vdj-dir results/sample/vdj \
    --out results/sample/valkyrn
```

Useful controls:

```text
--max-cdr3-distance 1     retained for CLI compatibility; not used for primary clone identity
--min-structure-family 3 minimum IGH family size for structural prioritization
--min-clonomap-family 100 minimum IGH family size for ClonoMap plots
--clonomap-k 30          PCA dimensions retained by ClonoMap
```

## Outputs

```text
valkyrn/
├── README.txt
├── valkyrn_cells.tsv
├── valkyrn_families.tsv
├── valkyrn_mutations.tsv
├── clonomap/
│   ├── clonomap_summary.tsv
│   └── <HC-family>/
│       ├── coords.tsv
│       ├── tree.tsv
│       ├── rows.tsv
│       ├── pca.png
│       └── mst.png
└── structure_candidates/
    ├── manifest.tsv
    └── paired_receptors.fasta
```

`README.txt` is a sample-level human-readable summary. The TSV files remain the implementation truth so that thresholds and biological interpretation can evolve without reducing the repertoire to hard PASS/FAIL labels.

## Why this is a separate Lumrik crate

Reconstruction and interpretation have different failure modes. `sc-vdj` must be conservative about what the reads prove. Valkyrn is allowed to compare cells, identify patterns, rank unusual families and formulate testable hypotheses. Keeping that boundary explicit prevents an attractive downstream biological story from changing the underlying receptor call.

The intended Lumrik flow is:

```text
reads
  -> nelrune / sc-vdj
  -> evidence-aware, error-corrected receptor reconstruction
  -> Valkyrn
       -> cell/receptor QC
       -> receptor families
       -> HC <-> LC relationships
       -> shared/private mutation patterns
       -> germline-allele hypotheses
       -> structurally interesting HC/LC pairs
  -> optional Norn structure process
       -> ABodyBuilder3 / future model
       -> structural comparison
       -> experimentally testable biological hypotheses
```

## Near-term development

The family builder deliberately reuses sc-vdj structural recombination identity rather than reclustering receptors from a lossy CDR3 summary. HC can define a lineage candidate; LC cannot define a multi-cell clone without a shared HC background. P/N-ambiguous HC calls have an explicit conservative fallback based on V/D/J and junction length information. The next important steps are to infer family ancestors within these structural families; detect probable unrepresented germline alleles across independent families; build mutation/lineage trees over the full reconstructed receptor, including established junction/CDR3 sequence; annotate HC lineage leaves with their observed LC partners; quantify LC-diverse versus LC-restricted HC families; and add structural comparison results back to the same family IDs.

Valkyrn should remain evidence-first: **hunt through immune repertoires for the interesting bastards, then say exactly why they are interesting.**
