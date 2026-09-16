# lumrik-status

> **Lumrik crate.** This crate is part of [Lumrik](../../README.md) and is distributed as part of the Lumrik workspace. See the [Lumrik README](../../README.md#license) and root `LICENSE` for the workspace licensing terms.

## What this crate does

Reusable lightweight live HTTP status/dashboard server for long-running Lumrik commands.

## Binaries

None. This crate is library-only.

## Library use

Implement `ServerContent` for your run state and start a `StatusServer` with `spawn_status_server`.
