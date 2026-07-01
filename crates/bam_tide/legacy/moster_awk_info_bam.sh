#!/usr/bin/env bash
set -euo pipefail

# bam_subset_diagnose.sh
# One-pass BAM/SAM text stream analysis using awk.
#
# Requirements:
#   - samtools
#   - awk with bitwise and() (gawk works; many awks do)
#
# Usage:
#   ./bam_subset_diagnose.sh <bam> [region] [topN]
#
# Examples:
#   ./bam_subset_diagnose.sh /tmp/pbmc_1k.bam
#   ./bam_subset_diagnose.sh /tmp/pbmc_1k.bam GL000220.1:109051-109100 50

BAM="${1:-}"
REGION="${2:-}"
TOPN="${3:-25}"

if [[ -z "${BAM}" || ! -f "${BAM}" ]]; then
  echo "Usage: $0 <bam> [region] [topN]" >&2
  exit 2
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in PATH" >&2
  exit 1
fi

AWK_BIN="${AWK_BIN:-awk}"

# Build samtools command
if [[ -n "${REGION}" ]]; then
  # Region must be passed to samtools view as a separate arg, not via grep.
  SAMCMD=(samtools view "${BAM}" "${REGION}")
else
  SAMCMD=(samtools view "${BAM}")
fi

# Stream once -> awk does everything
"${SAMCMD[@]}" | "${AWK_BIN}" -v TOPN="${TOPN}" '
function get_tag_i(field, prefix,   a) {
  # field like "NH:i:9"
  if (index(field, prefix) == 1) {
    split(field, a, ":");
    return a[3];
  }
  return "";
}

function cigar_sig(cig,   s) {
  # Coarse signature: replace run-lengths with op only, e.g. "56M1D35M" -> "MDM"
  # Also keep S/H/I/N/D/M/=X order intact.
  gsub(/[0-9]+/, "", cig);
  return cig;
}

function print_top_counts(title, arr, n,   k, i, line, cmd) {
  print "";
  print "== " title " (top " n ") ==";

  # We can’t sort assoc arrays in POSIX awk portably.
  # Use external sort by printing "count<TAB>key" to sort -nr.
  cmd = "sort -nr | head -n " n;
  for (k in arr) {
    printf("%d\t%s\n", arr[k], k) | cmd;
  }
  close(cmd);
}

BEGIN {
  FS = "\t";
  OFS = "\t";
  total = 0;
}

{
  total++;

  qname = $1;
  flag  = $2 + 0;
  rname = $3;
  pos1  = $4 + 0;   # SAM POS is 1-based in text
  mapq  = $5 + 0;
  cigar = $6;

  # ---- QNAME counts (dup detection) ----
  qcount[qname]++;

  # ---- FLAG-bit pattern summary ----
  sec  = and(flag, 256)  ? 1 : 0;
  supp = and(flag, 2048) ? 1 : 0;
  dup  = and(flag, 1024) ? 1 : 0;
  qcf  = and(flag, 512)  ? 1 : 0;
  unm  = and(flag, 4)    ? 1 : 0;

  fkey = "sec=" sec " supp=" supp " dup=" dup " qcf=" qcf " unm=" unm;
  flag_pat[fkey]++;

  # ---- MAPQ distribution ----
  mapq_cnt["MAPQ=" mapq]++;

  # ---- CIGAR signature distribution ----
  if (cigar != "*" && cigar != "") {
    sig = cigar_sig(cigar);
    cigar_cnt[sig]++;
  } else {
    cigar_cnt["*"]++;
  }

  # ---- NH/HI extraction ----
  nh = "NA";
  hi = "NA";
  for (i = 12; i <= NF; i++) {
    v = $i;
    t = get_tag_i(v, "NH:i:");
    if (t != "") nh = t;
    t = get_tag_i(v, "HI:i:");
    if (t != "") hi = t;
  }
  nhhi["NH=" nh " HI=" hi]++;

  # Also keep raw NH-only distribution if you want it
  nh_only["NH=" nh]++;

  # ---- Simple per-ref counter (optional) ----
  by_chr["CHR=" rname]++;
}

END {
  print "===== BAM subset diagnostic summary =====";
  print "Total SAM records:", total;

  print_top_counts("By chromosome (records)", by_chr, TOPN);
  print_top_counts("FLAG patterns", flag_pat, TOPN);
  print_top_counts("NH/HI pairs", nhhi, TOPN);
  print_top_counts("NH only", nh_only, TOPN);
  print_top_counts("MAPQ distribution", mapq_cnt, TOPN);
  print_top_counts("CIGAR signatures", cigar_cnt, TOPN);

  # QNAME duplicates: print only those with count>1, and show top TOPN by multiplicity
  print "";
  print "== QNAME duplicates (count > 1), top " TOPN " ==";
  cmd = "sort -nr | head -n " TOPN;
  dup_any = 0;
  for (k in qcount) {
    if (qcount[k] > 1) {
      dup_any = 1;
      printf("%d\t%s\n", qcount[k], k) | cmd;
    }
  }
  close(cmd);
  if (!dup_any) {
    print "(none)";
  }

  print "";
  print "Done.";
}
'

