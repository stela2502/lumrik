#!/usr/bin/env bash
set -euo pipefail

# ---- user tools (override via env if needed) ----
: "${BAMCOVERAGE_PY:=bamCoverage}"    # deeptools bamCoverage
: "${BAMCOVERAGE_RS:=bam-coverage}"   # your rust tool
: "${BWCOMPARE:=bw-compare}"          # your comparator

# ---- args ----
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.bam> [bin_width] [rust_options]"
  exit 2
fi

BAM="$1"
BIN_WIDTH="${2:-50}"
if [[ $# -ge 2 ]]; then
  shift 2
else
  shift $#
fi
RUST_OPTS="$@"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi

# ---- derive names ----
bam_base="$(basename "$BAM")"
bam_stem="${bam_base%.bam}"
PY_BW="/tmp/python_${bam_stem}.bw"
RS_BW="/tmp/rust_${bam_stem}.bw"
TMP_REPORT="/tmp/bw_compare_${bam_stem}_w${BIN_WIDTH}.txt"

# ---- timing/memory logs ----
METRICS_DIR="/tmp/bw_metrics_${bam_stem}_w${BIN_WIDTH}"
mkdir -p "$METRICS_DIR"
PY_TIME_LOG="${METRICS_DIR}/01_python_time.txt"
RS_TIME_LOG="${METRICS_DIR}/02_rust_time.txt"
CMP_TIME_LOG="${METRICS_DIR}/03_compare_time.txt"

# ---- pick GNU time ----
TIME_BIN=""
if command -v /usr/bin/time >/dev/null 2>&1; then
  TIME_BIN="/usr/bin/time"
elif command -v gtime >/dev/null 2>&1; then
  TIME_BIN="gtime"
else
  echo "ERROR: GNU time not found (/usr/bin/time or gtime). Install 'time' (GNU time)." >&2
  exit 1
fi

run_timed() {
  # Usage: run_timed <logfile> <label> -- <command...>
  local logfile="$1"; shift
  local label="$1"; shift
  if [[ "${1:-}" != "--" ]]; then
    echo "BUG: run_timed expects -- before command" >&2
    exit 99
  fi
  shift

  echo "---- ${label} ----" | tee "$logfile"
  echo "CMD: $*" | tee -a "$logfile"
  echo | tee -a "$logfile"

  # GNU time -v writes stats to stderr, so we redirect stderr into the logfile.
  # Stdout passes through to terminal (and can be captured by caller if desired).
  set +e
  "$TIME_BIN" -v "$@" 2>>"$logfile"
  local rc=$?
  set -e

  echo | tee -a "$logfile"
  echo "EXIT: $rc" | tee -a "$logfile"
  echo "-------------------" | tee -a "$logfile"
  echo
  return $rc
}

extract_metrics() {
  # Usage: extract_metrics <logfile>
  local logfile="$1"

  # wall time
  local wall
  wall=$(grep "Elapsed (wall clock) time" "$logfile" \
        | awk -F': ' '{print $2}' \
        | awk '
          {
            if (index($0, ":") == 0) {
              print $0
            } else {
              n=split($0,a,":");
              if (n==2) print a[1]*60+a[2];
              else if (n==3) print a[1]*3600+a[2]*60+a[3];
            }
          }')

  # CPU percent
  local cpu
  cpu=$(grep "Percent of CPU" "$logfile" \
        | awk '{print $8}' | tr -d '%')

  # max RSS in MB
  local rss
  rss=$(grep "Maximum resident set size" "$logfile" \
        | awk '{print $6}')
  rss=$(awk "BEGIN {printf \"%.1f\", $rss/1024}")

  echo "$wall $cpu $rss"
}

echo "BAM:        $BAM"
echo "BIN_WIDTH:  $BIN_WIDTH"
echo "PY_BW:      $PY_BW"
echo "RS_BW:      $RS_BW"
echo "REPORT:     $TMP_REPORT"
echo "METRICS:    $METRICS_DIR"
echo

# ---- 1) create expected python bw (if missing) ----
if [[ -s "$PY_BW" ]]; then
  echo "[1/3] Python BW exists, skipping: $PY_BW"
  echo "---- Python bamCoverage (skipped; file exists) ----" > "$PY_TIME_LOG"
  echo "PY_BW: $PY_BW" >> "$PY_TIME_LOG"
else
  echo "[1/3] Creating python BW with deeptools bamCoverage..."
  run_timed "$PY_TIME_LOG" "Python bamCoverage" -- \
    "$BAMCOVERAGE_PY" -b "$BAM" -o "$PY_BW" -bs "$BIN_WIDTH" --minMappingQuality 0
fi

# ---- 2) create rust bw (always regenerate, safer while iterating) ----
echo "[2/3] Creating rust BW with $BAMCOVERAGE_RS -b $BAM -o $RS_BW -w $BIN_WIDTH $RUST_OPTS ..."
run_timed "$RS_TIME_LOG" "Rust bam-coverage" -- \
  "$BAMCOVERAGE_RS" -b "$BAM" -o "$RS_BW" -w "$BIN_WIDTH" $RUST_OPTS

# ---- 3) compare ----
echo "[3/3] Comparing..."
set +e
OUT="$(
  run_timed "$CMP_TIME_LOG" "bw-compare" -- \
    "$BWCOMPARE" --python-bw "$PY_BW" --rust-bw "$RS_BW" --bin-width "$BIN_WIDTH" --eps 1e-4 2>&1
)"
RC=$?
set -e

echo "$OUT" | tee "$TMP_REPORT"

if [[ $RC -ne 0 ]]; then
  echo "bw-compare exited with code $RC" >&2
  exit $RC
fi

echo
echo "Saved:"
echo "  Compare output:  $TMP_REPORT"
echo "  Metrics logs:"
echo "    $PY_TIME_LOG"
echo "    $RS_TIME_LOG"
echo "    $CMP_TIME_LOG"
echo
echo "================ Performance Summary ================"

printf "%-12s %10s %8s %12s\n" "Tool" "Time(s)" "CPU(%)" "MaxRSS(MB)"
printf "%-12s %10s %8s %12s\n" "------------" "--------" "------" "-----------"

if [[ -f "$PY_TIME_LOG" ]]; then
  read py_wall py_cpu py_rss < <(extract_metrics "$PY_TIME_LOG")
  printf "%-12s %10.1f %8s %12.1f\n" "deeptools" "$py_wall" "$py_cpu" "$py_rss"
fi

if [[ -f "$RS_TIME_LOG" ]]; then
  read rs_wall rs_cpu rs_rss < <(extract_metrics "$RS_TIME_LOG")
  printf "%-12s %10.1f %8s %12.1f\n" "rust" "$rs_wall" "$rs_cpu" "$rs_rss"
fi

if [[ -f "$CMP_TIME_LOG" ]]; then
  read cmp_wall cmp_cpu cmp_rss < <(extract_metrics "$CMP_TIME_LOG")
  printf "%-12s %10.1f %8s %12.1f\n" "compare" "$cmp_wall" "$cmp_cpu" "$cmp_rss"
fi

echo "====================================================="
