#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <norn-output-root>"
    echo
    echo "Example:"
    echo "  $0 ~/sens05_home/NAS/NELRUNE/REAL_TEST_M39_v4_42cores"
    exit 1
fi

ROOT="$(realpath "$1")"


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONVERTER="$REPO_ROOT/target/release/convert-feature-observations"


SAMPLES=(
    ZD-4631-BcellsLaneF
    ZD-4631-PCLaneG
    ZD-4631-Undetermined
)

if [[ ! -x "$CONVERTER" ]]; then
    echo "ERROR: converter not found or not executable:"
    echo "  $CONVERTER"
    exit 1
fi

flatten_dir() {
    local parent="$1"
    local nested="$2"

    if [[ ! -d "$nested" ]]; then
        echo "  already flat / absent: $nested"
        return
    fi

    echo "  flattening: $nested -> $parent"

    # Refuse to overwrite anything already present.
    shopt -s dotglob nullglob
    for src in "$nested"/*; do
        dst="$parent/$(basename "$src")"

        if [[ -e "$dst" || -L "$dst" ]]; then
            echo "ERROR: refusing to overwrite existing:"
            echo "  $dst"
            exit 1
        fi

        mv -- "$src" "$dst"
    done
    shopt -u dotglob nullglob

    rmdir "$nested"
}

for sample in "${SAMPLES[@]}"; do
    SAMPLE="$ROOT/$sample"

    echo
    echo "================================================================"
    echo "$sample"
    echo "================================================================"

    [[ -d "$SAMPLE" ]] || {
        echo "ERROR: sample directory missing: $SAMPLE"
        exit 1
    }

    #
    # Flatten old Norn stage output directories.
    #
    flatten_dir "$SAMPLE/prepare" "$SAMPLE/prepare/prepare_out"
    flatten_dir "$SAMPLE/nelrune" "$SAMPLE/nelrune/nelrune_out"
    flatten_dir "$SAMPLE/vdj"     "$SAMPLE/vdj/vdj_out"

    #
    # Convert old feature_observations.bin exactly once.
    #
    OBS="$SAMPLE/prepare/feature_observations.bin"
    OLD="$SAMPLE/prepare/feature_observations_old.bin"

    if [[ -f "$OLD" ]]; then
        echo "  feature observations already converted:"
        echo "    legacy: $OLD"

        if [[ ! -f "$OBS" ]]; then
            echo "ERROR: legacy file exists but new feature_observations.bin is missing!"
            exit 1
        fi

        echo "    current: $OBS"
    else
        if [[ ! -f "$OBS" ]]; then
            echo "ERROR: feature observations missing:"
            echo "  $OBS"
            exit 1
        fi

        echo "  preserving legacy feature observations"
        mv -- "$OBS" "$OLD"

        echo "  converting -> new self-contained format"

        if ! "$CONVERTER" \
            --input "$OLD" \
            --out "$OBS"
        then
            echo "ERROR: conversion failed for $sample"
            echo "Legacy file is safe at:"
            echo "  $OLD"
            rm -f -- "$OBS"
            exit 1
        fi

        [[ -s "$OBS" ]] || {
            echo "ERROR: converter produced no usable output for $sample"
            exit 1
        }

        echo "  conversion complete"
    fi

    echo "  DONE: $sample"
done

echo
echo "================================================================"
echo "All samples converted and flattened."
echo "================================================================"
