# read-tag-table AI Contract

## Role

This crate stores/serializes per-read single-cell tags and deduplication identity used across preprocessing, mapping, and quantification.

## Round-trip identity

Serialization/deserialization must preserve cell sequence/quality, UMI sequence/quality, read identity, and grammar/provenance fields exactly according to the format contract.

QNAME encoding is a cross-crate interface. Changing it requires coordinated updates/tests in all producers and consumers; do not make an apparently local format change.

## Deduplication

Deduplication identity must not accidentally include/exclude fields that change biological molecule identity. Provenance is not disposable metadata.
