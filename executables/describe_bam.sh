#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
describe_bam.sh - Summarize a BAM file into a single, screen-friendly report.

Usage:
  describe_bam.sh -b <in.bam> [-o <report.txt>] [--no-index]
                  [--with-stats] [--head <N>]

Options:
  -b, --bam        Input BAM/CRAM file (required)
  -o, --out        Output report file (default: <bam>.report.txt)
  --no-index       Do not attempt to create an index if missing
  --with-stats     Also run 'samtools stats' (SLOW on large BAMs; default off)
  --head N         Number of 'samtools stats' lines to include (default: 60)
  -h, --help       Show help

Examples:
  describe_bam.sh -b sample.sorted.bam -o legacy/Bam_stats.txt
  describe_bam.sh -b sample.sorted.bam -o legacy/Bam_stats.txt --with-stats --head 120
EOF
}

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

# Print useful debug info on errors
trap 'ec=$?; echo "[ERROR] line $LINENO: $BASH_COMMAND (exit=$ec)" >&2; exit $ec' ERR

BAM=""
OUT=""
NO_INDEX=0
WITH_STATS=0
HEAD_N=60

while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--bam) BAM="${2:-}"; shift 2;;
    -o|--out) OUT="${2:-}"; shift 2;;
    --no-index) NO_INDEX=1; shift 1;;
    --with-stats) WITH_STATS=1; shift 1;;
    --head) HEAD_N="${2:-}"; shift 2;;
    -h|--help) usage; exit 0;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

[[ -n "${BAM}" ]] || { echo "ERROR: --bam is required" >&2; usage >&2; exit 2; }
[[ -f "${BAM}" ]] || { echo "ERROR: BAM not found: ${BAM}" >&2; exit 2; }

if [[ -z "${OUT}" ]]; then
  OUT="${BAM}.report.txt"
fi
OUT_DIR="$(dirname "$OUT")"
mkdir -p "$OUT_DIR"

need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing command: $1" >&2; exit 127; }; }
need_cmd samtools
need_cmd awk
need_cmd sed
need_cmd sort
need_cmd head
need_cmd date

# Create report early (so Ctrl-C still leaves something)
: > "${OUT}"
log "Writing report to: ${OUT}"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

hdr="${tmpdir}/header.txt"
flag="${tmpdir}/flagstat.txt"
idx="${tmpdir}/idxstats.txt"
stats="${tmpdir}/stats_head.txt"

# Index handling
index_exists() {
  [[ -f "${BAM}.bai" ]] || [[ -f "${BAM%.bam}.bai" ]] || [[ -f "${BAM}.crai" ]]
}

if [[ "${NO_INDEX}" -eq 0 ]]; then
  if ! index_exists; then
    log "Index not found; creating index..."
    samtools index -@ 4 "${BAM}"
  else
    log "Index found."
  fi
else
  log "Index creation disabled (--no-index)."
fi

# Write header section
{
  echo "BAM REPORT"
  echo "=========="
  echo "Generated: $(date -Is)"
  echo "Input:     ${BAM}"
  echo
} >> "${OUT}"

log "Collecting header..."
samtools view -H "${BAM}" > "${hdr}"

sort_order="$(grep -m1 '^@HD' "${hdr}" | sed -n 's/.*SO:\([^ \t]*\).*/\1/p' || true)"
pg_lines_count="$(grep -c '^@PG' "${hdr}" || true)"

{
  echo "Quick summary"
  echo "-------------"
  echo "Sort order:          ${sort_order:-unknown}"
  echo "Header @PG records:  ${pg_lines_count:-0}"
  echo
  echo "Header (first 80 lines)"
  echo "-----------------------"
  head -n 80 "${hdr}"
  echo
} >> "${OUT}"

log "Running flagstat..."
samtools flagstat "${BAM}" > "${flag}"

total_reads="$(awk '/in total/ {print $1; exit}' "${flag}" || true)"
mapped_reads="$(awk '/ mapped \(/ && !/primary/ {print $1; exit}' "${flag}" || true)"
mapped_pct="$(awk '/ mapped \(/ && !/primary/ {gsub(/[()%]/,"",$5); print $5; exit}' "${flag}" || true)"
paired="$(awk '/paired in sequencing/ {print $1; exit}' "${flag}" || true)"
properly_paired="$(awk '/properly paired/ {print $1; exit}' "${flag}" || true)"
dup="$(awk '/duplicates/ {print $1; exit}' "${flag}" || true)"

{
  echo "Read summary (from flagstat)"
  echo "----------------------------"
  echo "Total reads:         ${total_reads:-unknown}"
  echo "Mapped reads:        ${mapped_reads:-unknown} (${mapped_pct:-unknown}%)"
  echo "Paired in sequencing:${paired:-unknown}"
  echo "Properly paired:     ${properly_paired:-unknown}"
  echo "Duplicates:          ${dup:-unknown}"
  echo
  echo "flagstat"
  echo "--------"
  cat "${flag}"
  echo
} >> "${OUT}"

log "Running idxstats..."
if samtools idxstats "${BAM}" > "${idx}" 2> "${tmpdir}/idx.err"; then
  top_chr="${tmpdir}/top_chr.txt"
  awk 'BEGIN{OFS="\t"} $1!="*" {print $1,$3}' "${idx}" | sort -k2,2nr | head -n 15 > "${top_chr}"

  {
    echo "idxstats (top 15 chromosomes by mapped reads)"
    echo "--------------------------------------------"
    cat "${top_chr}"
    echo
  } >> "${OUT}"
else
  {
    echo "idxstats"
    echo "-------"
    echo "idxstats failed:"
    cat "${tmpdir}/idx.err"
    echo
  } >> "${OUT}"
fi

if [[ "${WITH_STATS}" -eq 1 ]]; then
  log "Running samtools stats (can be slow on large BAMs)..."
  samtools stats "${BAM}" | head -n "${HEAD_N}" > "${stats}"

  {
    echo "samtools stats (first ${HEAD_N} lines)"
    echo "------------------------------------"
    cat "${stats}"
    echo
  } >> "${OUT}"
else
  {
    echo "samtools stats"
    echo "-------------"
    echo "Skipped (run with --with-stats to include; can be slow on large BAMs)."
    echo
  } >> "${OUT}"
fi

log "DONE. Wrote report: ${OUT}"
echo "Wrote report: ${OUT}"

