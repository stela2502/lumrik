# sc_primer AI Contract

## Role

`sc_primer` defines single-cell chemistry/read-structure grammars and detects molecule identity/provenance. Its output is upstream biological metadata used by mapping and quantification.

## Grammar provenance

`GrammarType` is semantically meaningful and MUST survive downstream processing:

- `Gex` serializes as `G`.
- `Vdj` serializes as `V`.
- `Other` serializes as `O`.

Do not collapse these classes because they share barcode/UMI machinery.

BD Rhapsody GEX and VDJ/custom-capture structures are separate chemistries/grammars. A run may intentionally accept multiple chemistries at once; detection should preserve which grammar matched each molecule.

## Detection order

When multiple grammars are enabled, cheap/highly selective rejection (for example strong `FIXED` anchors) may be used to improve speed, but optimization MUST NOT change grammar semantics or provenance.

## Barcode and UMI identity

Cell and UMI extraction/correction are molecule identity operations. Downstream crates should consume the resulting identity rather than re-infer chemistry from mapped sequence or gene names.
