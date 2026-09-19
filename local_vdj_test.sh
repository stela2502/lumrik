#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Nelrune V(D)J local integration smoke test
#
# Uses the checked-in one-cell fixture, so this test is fast,
# deterministic, and does not depend on external FASTQ/reference data.
# It exercises the real nelrune-vdj binary and validates the four
# principal production outputs.
# ============================================================

FIXTURE="${FIXTURE:-crates/sc-vdj/tests/data/vdj-integration-cell}"
EXONIC="${EXONIC:-$FIXTURE/exonic}"
BAM="${BAM:-$FIXTURE/cell.bam}"
VDJ_INDEX="${VDJ_INDEX:-$FIXTURE/reference.vdjidx}"
OUT="${OUT:-target/sc-vdj-test-output/local-vdj-smoke}"
THREADS="${THREADS:-4}"

if [[ -n "${NELRUNE_VDJ:-}" ]]; then
    :
elif [[ -x "./target/release/nelrune-vdj" ]]; then
    NELRUNE_VDJ="./target/release/nelrune-vdj"
elif [[ -x "./target/release/nelrune-vdj" ]]; then
    NELRUNE_VDJ="./target/release/nelrune-vdj"
else
    NELRUNE_VDJ="$(command -v nelrune-vdj || true)"
fi

if [[ -z "$NELRUNE_VDJ" || ! -x "$NELRUNE_VDJ" ]]; then
    echo "ERROR: nelrune-vdj binary not found." >&2
    echo "Build it first, for example:" >&2
    echo "  cargo build --release --bin nelrune-vdj" >&2
    exit 1
fi

for f in "$BAM" "$VDJ_INDEX"; do
    [[ -f "$f" ]] || { echo "ERROR: fixture file not found: $f" >&2; exit 1; }
done
[[ -d "$EXONIC" ]] || { echo "ERROR: fixture exonic directory not found: $EXONIC" >&2; exit 1; }

rm -rf "$OUT"
mkdir -p "$OUT"

echo "HOST=$(hostname)"
echo "PWD=$PWD"
echo "nelrune-vdj=$NELRUNE_VDJ"
echo "fixture=$FIXTURE"
echo "output=$OUT"
echo

ARGS=(
    --exonic "$EXONIC"
    --bam "$BAM"
    --index "$VDJ_INDEX"
    --out "$OUT"
    --threads "$THREADS"
    --write-sequences
    --no-health-server
)

echo "nelrune-vdj command:"
printf ' %q' "$NELRUNE_VDJ" "${ARGS[@]}"
printf '\n\n'

"$NELRUNE_VDJ" "${ARGS[@]}"

required=(
    vdj_calls.tsv
    vdj_receptors.tsv
    airr_rearrangements.tsv
    vdj-mapping-info.txt
)

for name in "${required[@]}"; do
    path="$OUT/$name"
    [[ -s "$path" ]] || {
        echo "ERROR: expected non-empty output missing: $path" >&2
        exit 1
    }
done

# The established one-cell fixture contains exactly three receptor calls:
# one IGH, one IGK, and one IGL. Keep this assertion intentionally simple;
# the Rust integration test pins the detailed biological calls.
calls="$(( $(wc -l < "$OUT/vdj_calls.tsv") - 1 ))"
if (( calls != 3 )); then
    echo "ERROR: expected 3 V(D)J calls from one-cell fixture; got $calls" >&2
    exit 1
fi

for chain in IGH IGK IGL; do
    if ! awk -F '\t' -v chain="$chain" 'NR > 1 && $4 == chain { found=1 } END { exit !found }' "$OUT/vdj_calls.tsv"; then
        echo "ERROR: expected $chain call not found in $OUT/vdj_calls.tsv" >&2
        exit 1
    fi
done

echo
echo "V(D)J smoke test passed:"
echo "  calls: $calls (IGH + IGK + IGL)"
echo "  AIRR rows: $(( $(wc -l < "$OUT/airr_rearrangements.tsv") - 1 ))"
echo "  output: $OUT"
