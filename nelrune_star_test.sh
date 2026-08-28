#!/usr/bin/env bash


echo "HOST=$(hostname)"
echo "PATH=$PATH"

for exe in STAR gtf-splice-index nelrune; do
    printf '%-20s ' "$exe"
    command -v "$exe" || true
done


LOCAL="/home/stefanl/NAS/NELRUNE"

# ============================================================
# Nelrune test reference configuration
# ============================================================

# ------------------------------------------------------------
# Reference genome
# ------------------------------------------------------------

GENOME="/scale/gr01/shared/common/genome/genomes/mouse/GRCm38.p6/GRCm38.p6.genome.fa"
GTF="/scale/gr01/shared/common/genome/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf"


# ------------------------------------------------------------
# Generated indices
# ------------------------------------------------------------

# STAR genomeDir
STAR_INDEX="${LOCAL}/star_index"

# gtf-splice-index output used by Nelrune
SPLICE_INDEX="${LOCAL}/splice_index.idx"


# ------------------------------------------------------------
# Optional SNP calling
# ------------------------------------------------------------

#VCF="/path/to/variants.vcf.gz"



# ------------------------------------------------------------
# FASTQ input
#
# Bash expands these globs before passing them to Nelrune.
# Include BOTH the sample FASTQs and Undetermined FASTQs.
# ------------------------------------------------------------

#/home/stefanl/shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/DataDelivery_2026-05-29_13-28-30_snpseq01679/files/ZD-4631/20260522_LH00179_0469_B23K5C2LT3/Undetermined/Undetermined_S0_L007_R2_001.fastq.gz
#/home/stefanl/shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/DataDelivery_2026-05-29_13-28-30_snpseq01679/files/ZD-4631/20260522_LH00179_0469_B23K5C2LT3/Undetermined/Undetermined_S0_L007_R1_001.fastq.gz
#/home/stefanl/shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/DataDelivery_2026-05-29_13-28-30_snpseq01679/files/ZD-4631/20260522_LH00179_0469_B23K5C2LT3/Sample_ZD-4631-PCLaneG/ZD-4631-PCLaneG_S84_L007_R1_001.fastq.gz
#/home/stefanl/shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/DataDelivery_2026-05-29_13-28-30_snpseq01679/files/ZD-4631/20260522_LH00179_0469_B23K5C2LT3/Sample_ZD-4631-PCLaneG/ZD-4631-PCLaneG_S84_L007_R2_001.fastq.gz


FASTQ_ROOT="/home/stefanl/shared/jyuan/no_backup/giorgia_VDJ_single_cell_2026_06_01/DataDelivery_2026-05-29_13-28-30_snpseq01679/files/ZD-4631/20260522_LH00179_0469_B23K5C2LT3"

R1=(

    "$FASTQ_ROOT"/Sample_*/*_R1_*.fastq.gz
    "$FASTQ_ROOT"/Undetermined/*_R1_*.fastq.gz
)

R2=(
    "$FASTQ_ROOT"/Sample_*/*_R2_*.fastq.gz
    "$FASTQ_ROOT"/Undetermined/*_R2_*.fastq.gz
)


# ------------------------------------------------------------
# Fast features
# ------------------------------------------------------------

# Built-in BD sample tags (bd_sample_human, bd_sample_mouse)
ADDITIONAL_FEATURES=(
    bd_sample_mouse
)

# Later you can simply add FASTAs:
#
# ADDITIONAL_FEATURES+=(
#     "/path/to/hto.fa"
#     "/path/to/guides.fa"
# )

# Additional custom feature FASTAs can later be added here.
#
# e.g.
# GUIDE_FASTA="/path/to/guides.fa"
# HTO_FASTA="/path/to/hto.fa"


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT="./nelrune_test_out"


# ------------------------------------------------------------
# Programs
# ------------------------------------------------------------

NELRUNE="/scale/gr01/shared/common/software/Rustody/1.5/bin/nelrune"
STAR="STAR"

# ============================================================
# Nelrune server settings
# ============================================================

# Nelrune exposes a small live status dashboard while running.
#
# The server binds to 0.0.0.0 on the compute node, so this port
# must be unused on that node.
HEALTH_PORT=8080


# ------------------------------------------------------------
# Resources
# ------------------------------------------------------------

THREADS=16


# ------------------------------------------------------------
# STAR index
# ------------------------------------------------------------

if [[ ! -f "$STAR_INDEX/Genome" ]]; then
    echo "STAR index not found:"
    echo "  $STAR_INDEX"
    echo
    echo "Creating STAR index..."

    mkdir -p "$STAR_INDEX"

    "$STAR" \
        --runMode genomeGenerate \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME" \
        --sjdbGTFfile "$GTF"

    echo "STAR index created."
else
    echo "Using existing STAR index:"
    echo "  $STAR_INDEX"
fi


# ------------------------------------------------------------
# gtf-splice-index
# ------------------------------------------------------------

if [[ ! -f "$SPLICE_INDEX" ]]; then
    echo "Splice index not found:"
    echo "  $SPLICE_INDEX"
    echo
    echo "Creating splice index from:"
    echo "  $GTF"

    gtf-splice-index build \
        --annotation "$GTF" \
        --index "$SPLICE_INDEX"

    echo "Splice index created:"
    echo "  $SPLICE_INDEX"
else
    echo "Using existing splice index:"
    echo "  $SPLICE_INDEX"
fi

printf 'R1 files:\n'
printf '  %s\n' "${R1[@]}"

printf '\nR2 files:\n'
printf '  %s\n' "${R2[@]}"

printf '\nCounts: R1=%d R2=%d\n' \
    "${#R1[@]}" \
    "${#R2[@]}"

OPTIONAL_ARGS=()

if [[ -n "$VCF" ]]; then
    OPTIONAL_ARGS+=(
        --vcf "$VCF"
    )
fi
if [[ -n "$GENOME" ]]; then
    OPTIONAL_ARGS+=(
        --genome "$GENOME"
    )
fi
if (( ${#FAST_FEATURES[@]} > 0 )); then
    OPTIONAL_ARGS+=(
        --additional-features 
        "${ADDITIONAL_FEATURES[@]}"
    )
fi
if [[ -n "$HEALTH_PORT" ]]; then
    OPTIONAL_ARGS+=(
        --health-port "$HEALTH_PORT"
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

NELRUNE_ARGS+=(
    "${OPTIONAL_ARGS[@]}"
)

printf '\nNelrune command:\n'
printf ' %q' "$NELRUNE" "${NELRUNE_ARGS[@]}"
printf '\n\n'

"$NELRUNE" "${NELRUNE_ARGS[@]}"
