#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/sens05_shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/DataDelivery_2026-05-29_13-28-30_snpseq01679/files/ZD-4631/20260522_LH00179_0469_B23K5C2LT3"
LOCAL="$HOME/NAS/primer-stress-data"
NREADS=1000000
NLINES=$((NREADS * 4))

mkdir -p "$LOCAL"

prepare_pair() {
    local name="$1"
    local src_r1="$2"
    local src_r2="$3"

    local dst_r1="$LOCAL/${name}_R1.fastq.gz"
    local dst_r2="$LOCAL/${name}_R2.fastq.gz"

    if [[ ! -s "$dst_r1" ]]; then
        echo "Creating local $name R1 fixture..."
        set +o pipefail
        zcat "$src_r1" | head -n "$NLINES" | gzip > "$dst_r1"
        set -o pipefail
    fi

    if [[ ! -s "$dst_r2" ]]; then
        echo "Creating local $name R2 fixture..."
        set +o pipefail
        zcat "$src_r2" | head -n "$NLINES" | gzip > "$dst_r2"
        set -o pipefail
    fi
}


run_test() {
    local name="$1"

    echo
    echo "================================================================"
    echo "  $name"
    echo "================================================================"

    target/release/primer-stress \
        --r1 "$LOCAL/${name}_R1.fastq.gz" \
        --r2 "$LOCAL/${name}_R2.fastq.gz" \
        --chemistry bd-v2-384 bd-v2-384vdj \
        --max-reads "$NREADS"
}

prepare_pair \
    "PCLaneG" \
    "$BASE/Sample_ZD-4631-PCLaneG/ZD-4631-PCLaneG_S84_L007_R1_001.fastq.gz" \
    "$BASE/Sample_ZD-4631-PCLaneG/ZD-4631-PCLaneG_S84_L007_R2_001.fastq.gz"

prepare_pair \
    "BcellsLaneF" \
    "$BASE/Sample_ZD-4631-BcellsLaneF/ZD-4631-BcellsLaneF_S83_L007_R1_001.fastq.gz" \
    "$BASE/Sample_ZD-4631-BcellsLaneF/ZD-4631-BcellsLaneF_S83_L007_R2_001.fastq.gz"

prepare_pair \
    "Undetermined" \
    "$BASE/Undetermined/Undetermined_S0_L007_R1_001.fastq.gz" \
    "$BASE/Undetermined/Undetermined_S0_L007_R2_001.fastq.gz"

run_test "PCLaneG"
run_test "BcellsLaneF"
run_test "Undetermined"
