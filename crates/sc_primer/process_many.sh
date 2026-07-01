#!/usr/bin/env bash

set -euo pipefail

while IFS= read -r seq; do
    target/release/identify_primers \
        --chemistry bd-v2-384 \
        --seq "$seq"
done
