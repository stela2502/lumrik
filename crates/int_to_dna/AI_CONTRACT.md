# int_to_dna AI Contract

## Role

Compact DNA encoding/conversion primitives used in high-throughput paths.

## Encoding stability

Encoded representation, base ordering, reverse-complement behavior, and handling of non-ACGT/IUPAC symbols are API semantics. Do not change them for convenience without checking every persisted/indexed consumer.

## Performance

This crate is used in hot paths. Avoid unnecessary allocation/string conversion in core encoding operations.
