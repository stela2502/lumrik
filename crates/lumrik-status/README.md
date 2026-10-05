# lumrik-status

> **Lumrik crate.** This crate is part of [Lumrik](../../README.md) and is distributed as part of the Lumrik workspace. See the [Lumrik README](../../README.md#license) and root `LICENSE` for the workspace licensing terms.

## What this crate does

Reusable lightweight live HTTP status/dashboard server for long-running Lumrik commands.

## Binaries

None. This crate is library-only.

## Library use

Implement `ServerContent` for your run state and start a `StatusServer` with `spawn_status_server`.

The requested port is a starting port: if it is already in use, the server binds the first available port at or above it (up to 1000 ports higher). `StatusServer::addr()` reports the actual bound port. Port `0` retains the operating system's ephemeral-port behavior.
