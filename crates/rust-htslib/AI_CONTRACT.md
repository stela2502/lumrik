# rust-htslib AI Contract

## Role

This directory is a local/vendored `rust-htslib` dependency, not ordinary Lumrik application code.

## Modification policy

Do NOT modify this crate as part of a Lumrik feature/fix unless the task explicitly requires a fork-level HTSlib binding change and the reason has been established.

Prefer fixing Lumrik callers first. If a local rust-htslib modification is truly required, keep it minimal, document why upstream behavior is insufficient, and test compatibility with BAM/BCF consumers.
