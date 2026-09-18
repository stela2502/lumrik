# ClonoMap family model

The family layer owns biological family state. Plotting is downstream and must
not change membership, mutation measurements, or clone assignments.

## Ownership

```text
Analysis / FamilyCollection
│
├── Family (HC family)
│   ├── CellReceptor
│   │   ├── HC receptor
│   │   └── LC receptor observation(s)
│   ├── CellReceptor
│   │   ├── HC receptor
│   │   └── LC receptor observation(s)
│   └── LC clones
│       ├── LightClone
│       └── LightClone
├── Family (HC family)
└── unresolved CellReceptors
```

`CellReceptor` owns one cell's receptor evidence. HC and LC receptor records each
carry their own observed sequence and the Lumrik `naive_recombination` for that
cell's observed span. A compact recombination ID is clone/recombination identity;
it is not a globally unique row key and does not imply identical sequence span.

`Family` owns accepted HC-family members and the LC clone split inside that HC
family. It does **not** own rejected/ejected cells. `Family::align()` removes bad
candidates and returns the bare `CellReceptor`s to the outer analysis layer.
The outer layer decides whether another HC family can integrate them or whether
they remain unresolved.

## HC lifecycle

1. Construct `CellReceptor`s. Exactly one valid productive HC is required.
2. Build provisional HC families using structural information only: same IGH V,
   same J, and complete-link CDR3 edit distance <= configured maximum.
3. `Family::align()` measures each cell's observed HC against that cell's own
   naive reconstruction span. Missing terminal coverage is not mutation.
4. Alignment failures and mutation-depth outliers are returned to the caller.
5. The outer analysis tries each rejected `CellReceptor` against other HC
   families. `Family::try_integrate()` checks HC V/J compatibility and then the
   hard mutation/alignment gate; CDR3 distance is not a second rescue veto.
6. HC family membership is then final.

## LC lifecycle

LC splitting happens only inside a finalized HC `Family`; LC evidence can never
move a cell between HC families.

Initial LC clone collection uses locus (IGK/IGL), V, J, and strict CDR3 distance.
Mutation alignment is then measured from each LC observation against its own
naive reconstruction span. LC mutation outliers are removed from that LC clone,
try the other compatible LC clones in the same HC family, and if none accepts
them they become the root of a new LC clone. They are not discarded merely for
failing an existing LC clone.

## Import boundary

The standalone runner consumes a Lumrik VDJ output directory. AIRR supplies the
standardized cell/locus/V/D/J/C/CDR3/observed sequence fields; `vdj_calls.tsv`
supplies `naive_recombination`. These files are joined by
`(cell_id, lumrik_recombination_id) <-> (cell, recombination_id)`, never by
recombination ID alone.
