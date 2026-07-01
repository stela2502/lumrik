#!/usr/bin/env bash
set -u  # (no -e on purpose)

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <bam> [bin_width]"
  exit 1
fi

BAM="$1"
BIN_WIDTH="${2:-50}"

bam_base="$(basename "$BAM")"
bam_stem="${bam_base%.bam}"
REPORT="/tmp/bw_compare_${bam_stem}_w${BIN_WIDTH}.txt"

FLAGS=(
  0        # none
  256      # secondary
  512      # QC-fail
  1024     # duplicate
  2048     # supplementary
  768      # secondary + QC-fail
  1280     # secondary + duplicate
  2304     # secondary + supplementary
  2816     # secondary + QC-fail + supplementary
  3840     # secondary + QC-fail + duplicate + supplementary
)

echo "FLAG    TOTAL"
echo "--------------------------------------------------------------"

for f in "${FLAGS[@]}"; do
  # Run compare; don't abort script if it fails
  OUT=$(./run_compare.sh "$BAM" "$BIN_WIDTH" --sam-flag-exclude $f  | grep "^TOTAL" )
  rc=$?

  if [[ $rc -ne 0 ]]; then
    printf "%-6s  ERROR (run_compare.sh exit %d)\n" "$f" "$rc"
    echo "./run_compare.sh $BAM $BIN_WIDTH $f"
    #echo "$OUT" | tail -n "$OUT"
    continue
  fi

  if [[ -z "$OUT" ]]; then
    printf "%-6s  ERROR (no TOTAL line found)\n" "$f"
  else
    printf "%-6s  %s\n" "$f" "$OUT"
  fi
done
