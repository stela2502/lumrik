# scdata design contract

`scdata` owns sparse single-cell storage. Its public API should make ownership clear and should not force callers to reconstruct storage logic.

## The layers

### `CellHash<T>` is a public generic storage primitive

`CellHash` is intentionally public. It is the reusable sparse per-cell container for Lumrik crates that need to store their own cell-associated data. External crates may provide their own value type through the existing trait contract.

It is not an implementation detail of `Scdata` that we are trying to hide.

### `Scdata` is the feature/UMI matrix

`Scdata` owns one sparse cell-by-feature data set and the operations that naturally belong to that matrix: inserting observations, merging, filtering/finalizing, deriving counts, converting to analysis representations, and reading/writing matrix formats.

A standalone owner of an `Scdata` may use that API directly.

### `QuantData` owns a named collection of `Scdata`

When an `Scdata` is owned by `QuantData`, it does not escape from `QuantData`. There is no public `get()` or `get_mut()` returning the stored matrix.

Callers ask `QuantData` to perform an operation on a named data set and receive the result of that operation, not the underlying `Scdata`.

Examples:

- add an observation to `"exonic"`;
- ask for cell or feature counts from `"exonic"`;
- merge another `QuantData`;
- attach a completed externally produced `Scdata` under a unique name;
- apply one caller-provided allowed-cell set across all owned matrices;
- write all named matrices in one coordinated operation.

Dataset names are identities. Adding a second dataset with an existing name is a programming error and panics rather than replacing data silently.

## Run state is not quantification data

`MappingInfo`, timers, counters, progress state, and similar run bookkeeping do **not** belong to `QuantData`.

The runner/processor owns run state. If an `Scdata` operation needs a `MappingInfo` while inserting or merging, `QuantData` may accept it for that operation or return merge accounting to the caller. It must not retain the run report as part of the biological data object.

## Cell selection belongs to the caller; consistent application belongs to `QuantData`

A runner may choose which biological signal defines the canonical cells. `QuantData` may derive a cell set from a named matrix when asked, but the caller owns the policy decision.

Once the caller supplies the canonical cell set for export, `QuantData` applies it consistently to every matrix it owns. Runners should not fetch individual matrices and repeat the filtering/export lifecycle themselves.

Consistent application does **not** mean forcing every sparse matrix to have identical columns or rows. Every `Scdata` owns its own observed barcode and feature set. Raw output writes exactly that dataset's observed sparse data. Filtered output restricts each dataset to the caller's allowed cells, so absent cells remain absent rather than being manufactured as empty sparse columns.

## Feature/row order is supplied from outside

`QuantData` owns matrix values, not the biological reference that defines canonical feature order and feature serialization.

Writing therefore requires the caller to provide the row/index information for every named matrix being written. Missing row/index information is an error. The writing implementation belongs in one place so different runners cannot silently implement different matrix ordering or filtering rules.

## Analysis code gets analysis representations, not storage internals

If another crate needs data for analysis, prefer a deliberate representation or derived result such as counts, cell IDs, feature IDs, or a sparse matrix conversion (`sprs`) over exposing the internally owned `Scdata`/`CellData` structure.

Do not add getters merely because a field exists.

## Practical rule

**`CellHash` is reusable storage. `Scdata` owns one matrix. `QuantData` coordinates named matrices. The runner owns run policy and run state.**

Keep the normal call site short. Add API only for operations real callers need.
