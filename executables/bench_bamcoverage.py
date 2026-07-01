#!/usr/bin/env python3
import argparse
import os, re, shlex, subprocess, sys, time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple, List


@dataclass
class TimeMetrics:
    elapsed_seconds: Optional[float]
    max_rss_kb: Optional[int]
    raw_time_stderr: str


_ELAPSED_RE = re.compile(r"^Elapsed \(wall clock\) time .*:\s*(.+)\s*$", re.MULTILINE)
_MAXRSS_RE = re.compile(r"^Maximum resident set size \(kbytes\):\s*(\d+)\s*$", re.MULTILINE)


def parse_elapsed_to_seconds(s: str) -> Optional[float]:
    """
    GNU time -v elapsed formats:
      - "m:ss.xx"
      - "h:mm:ss.xx"
      - sometimes "d-hh:mm:ss" in some builds (rare)
    """
    s = s.strip()
    if not s:
        return None

    # Handle possible "d-hh:mm:ss.xx"
    if "-" in s:
        day_part, rest = s.split("-", 1)
        try:
            days = int(day_part)
        except ValueError:
            return None
    else:
        days = 0
        rest = s

    parts = rest.split(":")
    try:
        if len(parts) == 2:
            m = int(parts[0])
            sec = float(parts[1])
            return days * 86400 + m * 60 + sec
        if len(parts) == 3:
            h = int(parts[0])
            m = int(parts[1])
            sec = float(parts[2])
            return days * 86400 + h * 3600 + m * 60 + sec
    except ValueError:
        return None

    return None



@dataclass
class TimeMetrics:
    elapsed_seconds: Optional[float]
    max_rss_kb: Optional[int]
    raw_time_stderr: str  # keep name for compatibility; here we store stderr of the program


_RSS_RE = re.compile(r"^VmRSS:\s+(\d+)\s+kB$", re.MULTILINE)

def _read_rss_kb(pid: int) -> Optional[int]:
    """
    Linux-only. Reads current RSS (kB) from /proc/<pid>/status.
    Returns None if process already ended or /proc not available.
    """
    try:
        txt = Path(f"/proc/{pid}/status").read_text()
    except (FileNotFoundError, ProcessLookupError, PermissionError):
        return None

    m = _RSS_RE.search(txt)
    if not m:
        return None
    return int(m.group(1))

def run_with_time(cmd: List[str], log_prefix: Path, time_bin: Optional[str] = None,
                  poll_interval: float = 0.05) -> Tuple[int, str, TimeMetrics]:
    """
    Run cmd, measure wall time + peak RSS (kB) by polling /proc.
    No /usr/bin/time needed. Linux-only for RSS.
    """
    log_prefix.parent.mkdir(parents=True, exist_ok=True)

    # ---- ONE LINE: print command to stderr ----
    print("[RUN]", " ".join(shlex.quote(x) for x in cmd), file=sys.stderr, flush=True)

    # Start subprocess
    start = time.perf_counter()
    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert proc.stdout is not None and proc.stderr is not None

    peak_rss = None

    # Poll RSS while running
    while True:
        rc = proc.poll()
        rss = _read_rss_kb(proc.pid)
        if rss is not None:
            peak_rss = rss if peak_rss is None else max(peak_rss, rss)
        if rc is not None:
            break
        time.sleep(poll_interval)

    stdout, stderr = proc.communicate()
    elapsed = time.perf_counter() - start

    # Logs
    (log_prefix.with_suffix(".cmd.txt")).write_text(" ".join(shlex.quote(x) for x in cmd) + "\n")
    (log_prefix.with_suffix(".stdout.txt")).write_text(stdout or "")
    (log_prefix.with_suffix(".stderr.txt")).write_text(stderr or "")

    metrics = TimeMetrics(
        elapsed_seconds=elapsed,
        max_rss_kb=peak_rss,
        raw_time_stderr=stderr or ""
    )
    return proc.returncode, (stdout or ""), metrics


def run_plain(cmd: List[str], log_prefix: Path) -> Tuple[int, str, str]:
    log_prefix.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, text=True, capture_output=True)
    (log_prefix.with_suffix(".cmd.txt")).write_text(" ".join(shlex.quote(x) for x in cmd) + "\n")
    (log_prefix.with_suffix(".stdout.txt")).write_text(proc.stdout)
    (log_prefix.with_suffix(".stderr.txt")).write_text(proc.stderr)
    return proc.returncode, proc.stdout, proc.stderr


def extract_total_line(text: str) -> Optional[str]:
    for line in text.splitlines():
        if line.startswith("TOTAL"):
            return line
    return None


def main():
    ap = argparse.ArgumentParser(description="Benchmark bamCoverage vs bam-coverage and compare BigWigs.")
    ap.add_argument("--bam", required=True, help="Input BAM file")
    ap.add_argument("--flags", required=True, help="Comma-separated list of flags, e.g. 0,256,512")
    ap.add_argument("--bin-size", type=int, default=None, help="Optional: pass bin size to bamCoverage/bam-coverage (if you want)")
    ap.add_argument("--bamCoverage", default="bamCoverage", help="Path/name for deeptools bamCoverage")
    ap.add_argument("--bam-coverage", dest="bam_coverage", default="bam-coverage", help="Path/name for rust bam-coverage")
    ap.add_argument("--bw-compare", dest="bw_compare", default="bw-compare", help="Path/name for bw-compare")
    ap.add_argument("--out", default="results.tsv", help="Output TSV")
    ap.add_argument("--log-dir", default="bench_logs", help="Directory for per-run logs")
    args = ap.parse_args()

    bam = Path(args.bam)
    if not bam.exists():
        raise SystemExit(f"ERROR: BAM not found: {bam}")

    flags = []
    for x in args.flags.split(","):
        x = x.strip()
        if not x:
            continue
        try:
            flags.append(int(x))
        except ValueError:
            raise SystemExit(f"ERROR: invalid flag: {x}")

    core = bam.name
    if core.endswith(".bam"):
        core = core[:-4]

    log_dir = Path(args.log_dir)

    rows = []
    header = [
        "flag",
        "CHR", 
        "n_over_eps",
        "frac_n_over_eps",
        "mean_abs",
        "var_abs",
        "rmse",
        "max_abs", 
        "pearson_rho",
        "py_elapsed_s",
        "py_max_rss_kb",
        "rs_elapsed_s",
        "rs_max_rss_kb",
        "python_bw",
        "rust_bw",
    ]

    for flag in flags:
        python_bw = Path(f"python_{core}_flag{flag}.bw")
        rust_bw = Path(f"rust_{core}_flag{flag}.bw")

        # 1) Python deeptools bamCoverage
        py_cmd = [args.bamCoverage, "-b", str(bam), "--samFlagExclude", str(flag), "--outFileName", str(python_bw)]
        if args.bin_size is not None:
            # deeptools uses --binSize
            py_cmd += ["--binSize", str(args.bin_size)]
        rc, _out, py_metrics = run_with_time(py_cmd, log_dir / f"flag{flag}.python")
        if rc != 0:
            print(f"[flag {flag}] WARNING: bamCoverage exited with code {rc} (see logs)")

        # 2) Rust bam-coverage
        rs_cmd = [args.bam_coverage, "-b", str(bam), "--sam-flag-exclude", str(flag), "-o", str(rust_bw)]
        if args.bin_size is not None:
            # if your rust tool uses --bin-width or -bs, change here accordingly.
            rs_cmd += ["-w", str(args.bin_size)]
        rc, _out, rs_metrics = run_with_time(rs_cmd, log_dir / f"flag{flag}.rust")
        if rc != 0:
            print(f"[flag {flag}] WARNING: bam-coverage exited with code {rc} (see logs)")

        # 3) bw-compare (collect ^TOTAL line)
        cmp_cmd = [args.bw_compare, "--python-bw", str(python_bw), "--rust-bw", str(rust_bw)]
        rc, cmp_stdout, _cmp_stderr = run_plain(cmp_cmd, log_dir / f"flag{flag}.compare")
        if rc != 0:
            print(f"[flag {flag}] WARNING: bw-compare exited with code {rc} (see logs)")
        total_line = extract_total_line(cmp_stdout)

        rows.append([
            str(flag),
            total_line or "",
            "" if py_metrics.elapsed_seconds is None else f"{py_metrics.elapsed_seconds:.3f}",
            "" if py_metrics.max_rss_kb is None else str(py_metrics.max_rss_kb),
            "" if rs_metrics.elapsed_seconds is None else f"{rs_metrics.elapsed_seconds:.3f}",
            "" if rs_metrics.max_rss_kb is None else str(rs_metrics.max_rss_kb),
            str(python_bw),
            str(rust_bw),
        ])

    out_path = Path(args.out)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"Wrote {out_path} with {len(rows)} rows.")
    print(f"Logs in: {log_dir}/")


if __name__ == "__main__":
    main()

