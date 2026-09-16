# sc-mapper

> **Lumrik crate.** This crate is part of [Lumrik](../../README.md) and is distributed as part of the Lumrik workspace. See the [Lumrik README](../../README.md#license) and root `LICENSE` for the workspace licensing terms.

## What this crate does

Streaming integration layer for external aligners such as STAR, minimap2 and BWA.

## Binaries

- `mapper-wrapper`

## Library use

Use `StreamingMapper`/`StreamingMapperCli` and the mapper implementations (`Star`, `Minimap2`, `Bwa`) to feed reads to an external mapper and consume mapping clusters without materializing an intermediate FASTQ.
