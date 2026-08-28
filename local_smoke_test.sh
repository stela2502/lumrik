#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Nelrune local BD Rhapsody smoke test
#
# Processes at most 100,000 RAW read pairs PER R1/R2 input pair.
# With the two lanes below that means at most 200,000 raw pairs.
# ============================================================

ROOT="/data2/Elena"

R1=(
    "$ROOT/Rhapsody_229EB_20240625_S1_L001_R1_001.fastq.gz"
    "$ROOT/Rhapsody_229EB_20240625_S1_L002_R1_001.fastq.gz"
)

R2=(
    "$ROOT/Rhapsody_229EB_20240625_S1_L001_R2_001.fastq.gz"
    "$ROOT/Rhapsody_229EB_20240625_S1_L002_R2_001.fastq.gz"
)

# ------------------------------------------------------------
# Nelrune / mapper
# ------------------------------------------------------------

THREADS="${THREADS:-16}"
MAX_READS="${MAX_READS:-100000}"

# Prefer the freshly built local binary.
if [[ -x "./target/x86_64-unknown-linux-musl/release/nelrune" ]]; then
    NELRUNE="./target/x86_64-unknown-linux-musl/release/nelrune"
elif [[ -x "./target/release/nelrune" ]]; then
    NELRUNE="./target/release/nelrune"
else
    NELRUNE="$(command -v nelrune || true)"
fi

STAR="${STAR:-STAR}"

# ------------------------------------------------------------
# References
#
# Set these to your real local indices before running, e.g.
#
#   export STAR_INDEX=/data2/.../star_index
#   export SPLICE_INDEX=/data2/.../splice_index.idx
#
# ------------------------------------------------------------
STAR_INDEX=/data2/Elena/star_index
SPLICE_INDEX=/data2/Elena/splice_index.idx


# Optional, only needed if you actually want SNP/reference-aware work.
GENOME="${GENOME:-}"
VCF="${VCF:-}"

# ------------------------------------------------------------
# Additional short features
# ------------------------------------------------------------

ADDITIONAL_FEATURES=(
    bd_sample_mouse
)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT="${OUT:-/data2/Elena/nelrune_elena_100k}"

# ============================================================
# Sanity checks
# ============================================================

echo "HOST=$(hostname)"
echo "PWD=$PWD"
echo "PATH=$PATH"
echo

for exe in "$STAR" "$NELRUNE"; do
    printf '%-24s ' "$exe"
    if [[ "$exe" == */* ]]; then
        [[ -x "$exe" ]] && echo "$exe" || {
            echo "NOT FOUND / NOT EXECUTABLE"
            exit 1
        }
    else
        command -v "$exe" || {
            echo "NOT FOUND"
            exit 1
        }
    fi
done

if [[ -z "$STAR_INDEX" ]]; then
    echo
    echo "ERROR: STAR_INDEX is not set."
    echo "Example:"
    echo "  export STAR_INDEX=/data2/reference/star_index"
    exit 1
fi

if [[ -z "$SPLICE_INDEX" ]]; then
    echo
    echo "ERROR: SPLICE_INDEX is not set."
    echo "Example:"
    echo "  export SPLICE_INDEX=/data2/reference/splice_index.idx"
    exit 1
fi

if [[ ! -d "$STAR_INDEX" ]]; then
    echo "ERROR: STAR index directory not found: $STAR_INDEX" >&2
    exit 1
fi

if [[ ! -f "$SPLICE_INDEX" ]]; then
    echo "ERROR: splice index not found: $SPLICE_INDEX" >&2
    exit 1
fi

for f in "${R1[@]}" "${R2[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: FASTQ not found: $f" >&2
        exit 1
    fi
done

if (( ${#R1[@]} != ${#R2[@]} )); then
    echo "ERROR: R1/R2 count mismatch: ${#R1[@]} vs ${#R2[@]}" >&2
    exit 1
fi

echo
printf 'R1 files:\n'
printf '  %s\n' "${R1[@]}"

echo
printf 'R2 files:\n'
printf '  %s\n' "${R2[@]}"

echo
printf 'Input pairs: %d\n' "${#R1[@]}"
#printf 'Raw-read limit per pair: %s\n' "$MAX_READS"
printf 'Maximum raw pairs total: %d\n' "$(( MAX_READS * ${#R1[@]} ))"
printf 'Threads: %s\n' "$THREADS"
printf 'Output: %s\n' "$OUT"

# ============================================================
# Build command
# ============================================================

OPTIONAL_ARGS=()

if [[ -n "$VCF" ]]; then
    OPTIONAL_ARGS+=(--vcf "$VCF")

    if [[ -n "$GENOME" ]]; then
        OPTIONAL_ARGS+=(--genome "$GENOME")
    else
        echo "ERROR: VCF was supplied but GENOME is empty." >&2
        exit 1
    fi
fi

if (( ${#ADDITIONAL_FEATURES[@]} > 0 )); then
    OPTIONAL_ARGS+=(
        --additional-features
        "${ADDITIONAL_FEATURES[@]}"
    )
fi

NELRUNE_ARGS=(
    --r1 "${R1[@]}"
    --r2 "${R2[@]}"

    --chemistry bd-v2-384

    --mapper star
    --mapper-index "$STAR_INDEX"
    --mapper-threads "$THREADS"

    --threads "$THREADS"

    --index "$SPLICE_INDEX"
    --outpath "$OUT"
)

NELRUNE_ARGS+=("${OPTIONAL_ARGS[@]}")

mkdir -p "$OUT"

echo
echo "Nelrune command:"
printf ' %q' "$NELRUNE" "${NELRUNE_ARGS[@]}"
printf '\n\n'

"$NELRUNE" "${NELRUNE_ARGS[@]}"

