#!/usr/bin/env bash
set -euo pipefail

BAM="${1:?usage: bench.sh <sorted.bam>}"
OUTDIR="${2:-bench_out}"
BIN="${3:-50}"
THREADS="${4:-1}"

mkdir -p "$OUTDIR"

RUST_EXE="./target/release/bam-coverage"
PY_EXE="bamCoverage"

RUST_OUT="$OUTDIR/rust.bedgraph"
PY_OUT="$OUTDIR/python.bedgraph"

echo "== Inputs =="
echo "BAM: $BAM"
echo "BIN: $BIN"
echo "THREADS: $THREADS"
echo "OUTDIR: $OUTDIR"
echo

rm -f "$RUST_OUT" "$PY_OUT"

echo "== Rust (bam-coverage) =="
/usr/bin/time -v -o "$OUTDIR/time_rust.txt" \
  "$RUST_EXE" \
  -b "$BAM" \
  -o "$RUST_OUT" \
  -w "$BIN" \
  -n not \
  --min-mapping-quality 0 \
  2> "$OUTDIR/log_rust.txt"

ls -lh "$RUST_OUT" | tee "$OUTDIR/size_rust.txt"
echo

echo "== Python (deepTools bamCoverage) =="
/usr/bin/time -v -o "$OUTDIR/time_python.txt" \
  "$PY_EXE" \
  -b "$BAM" \
  -o "$PY_OUT" \
  --outFileFormat bedgraph \
  --binSize "$BIN" \
  --numberOfProcessors "$THREADS" \
  2> "$OUTDIR/log_python.txt"

ls -lh "$PY_OUT" | tee "$OUTDIR/size_python.txt"
echo

echo "== Extracted summary =="
echo "-- Rust time/mem --"
grep -E "Elapsed|User time|System time|Maximum resident set size" "$OUTDIR/time_rust.txt" || true
echo
echo "-- Python time/mem --"
grep -E "Elapsed|User time|System time|Maximum resident set size" "$OUTDIR/time_python.txt" || true

