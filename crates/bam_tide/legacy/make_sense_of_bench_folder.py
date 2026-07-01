#!/usr/bin/env python3
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Optional, Dict, List, Tuple


# ----------------------------
# Exceptions with full context
# ----------------------------

@dataclass
class ParseContext:
    path: Path
    lineno: int
    line: str
    key: Optional[str] = None
    value: Optional[str] = None
    rule: Optional[str] = None   # e.g. "rsplit(':',1)" or "elapsed split ':'"
    note: Optional[str] = None   # hint / explanation


class BenchParseError(Exception):
    def __init__(self, message: str, ctx: Optional[ParseContext] = None, cause: Optional[Exception] = None):
        super().__init__(message)
        self.ctx = ctx
        self.cause = cause


# ----------------------------
# Parsing helpers
# ----------------------------

def parse_elapsed(value: str, ctx: ParseContext) -> float:
    """
    GNU time formats commonly seen:
      - m:ss.xx   (minutes can be > 59)
      - h:mm:ss.xx
    We parse by splitting on ':'.
    """
    s = value.strip()
    parts = s.split(":")

    ctx = ParseContext(**{**ctx.__dict__, "rule": "elapsed split(':')", "value": value})

    try:
        if len(parts) == 2:
            # m:ss(.xx)
            m = int(parts[0])
            sec = float(parts[1])
            return m * 60.0 + sec
        if len(parts) == 3:
            # h:mm:ss(.xx)
            hours_str, minutes_str, sec_str = parts
            hours = int(hours_str)
            minutes = int(minutes_str)
            seconds = float(sec_str)
            return hours * 3600.0 + minutes * 60.0 + seconds

    except ValueError as e:
        raise BenchParseError(
            "Failed to parse elapsed time value as m:ss(.xx) or h:mm:ss(.xx).",
            ctx=ctx,
            cause=e,
        )

    raise BenchParseError(
        f"Unexpected elapsed time format (expected m:ss(.xx) or h:mm:ss(.xx)): '{s}'",
        ctx=ctx,
    )


def parse_time_v(path: Path) -> Dict[str, float]:
    """
    Parse GNU /usr/bin/time -v output file.

    We match keys robustly by prefix for the elapsed line; others are exact.
    On any problem, raise BenchParseError with full line context.
    """
    if not path.exists():
        raise BenchParseError(f"Missing timing file: {path}")

    lines = path.read_text(errors="replace").splitlines()
    out_raw: Dict[str, str] = {}

    for i, line in enumerate(lines, start=1):
        if ":" not in line:
            continue

        # Split at LAST colon because the elapsed key itself contains colons.
        try:
            if ": " not in line:
                continue
            k, v = line.split(": ", 1)
            k = k.strip()
            v = v.strip()
        except Exception as e:
            raise BenchParseError(
                "Failed to split line using rsplit(':', 1).",
                ctx=ParseContext(path=path, lineno=i, line=line, rule="rsplit(':',1)"),
                cause=e,
            )

        k = k.strip()
        v = v.strip()

        # Robust matching
        if k.startswith("Elapsed (wall clock) time"):
            out_raw["elapsed"] = v
        elif k == "User time (seconds)":
            out_raw["user_s"] = v
        elif k == "System time (seconds)":
            out_raw["sys_s"] = v
        elif k == "Maximum resident set size (kbytes)":
            out_raw["max_rss_kb"] = v

    required = ["elapsed", "user_s", "sys_s", "max_rss_kb"]
    missing = [k for k in required if k not in out_raw]
    if missing:
        # Provide a short preview and also tell how many lines were scanned
        preview = "\n".join(f"{n:>4}: {lines[n-1]}" for n in range(1, min(25, len(lines)) + 1))
        raise BenchParseError(
            f"Missing required keys {missing} in {path} (scanned {len(lines)} lines).",
            ctx=ParseContext(
                path=path,
                lineno=0,
                line="",
                rule="key collection",
                note=f"Preview (first {min(25, len(lines))} lines):\n{preview}",
            ),
        )

    # Convert to floats with per-field context
    def to_float(field: str) -> float:
        raw = out_raw[field]
        # We don't have the original line number here; still include field + value.
        ctx = ParseContext(path=path, lineno=0, line="", key=field, value=raw, rule="float()")
        try:
            return float(raw)
        except ValueError as e:
            raise BenchParseError(f"Failed to parse numeric field '{field}' as float.", ctx=ctx, cause=e)

    elapsed_s = parse_elapsed(
        out_raw["elapsed"],
        ParseContext(path=path, lineno=0, line="", key="elapsed", value=out_raw["elapsed"]),
    )

    user_s = to_float("user_s")
    sys_s = to_float("sys_s")

    # max_rss_kb is integer-like but float() is fine
    max_rss_kb = to_float("max_rss_kb")
    max_rss_mb = max_rss_kb / 1024.0

    return {
        "elapsed_s": elapsed_s,
        "user_s": user_s,
        "sys_s": sys_s,
        "max_rss_mb": max_rss_mb,
    }


def parse_size_mb(path: Path) -> float:
    if not path.exists():
        raise BenchParseError(f"Missing output file: {path}")
    return path.stat().st_size / (1024.0 * 1024.0)


def load_run(run_dir: Path) -> Dict[str, float]:
    if not run_dir.is_dir():
        raise BenchParseError(f"Not a directory: {run_dir}")

    rust_time = parse_time_v(run_dir / "time_rust.txt")
    py_time = parse_time_v(run_dir / "time_python.txt")

    rust_out = run_dir / "rust.bedgraph"
    py_out = run_dir / "python.bedgraph"

    rust_elapsed = rust_time["elapsed_s"]
    py_elapsed = py_time["elapsed_s"]
    if rust_elapsed <= 0:
        raise BenchParseError(f"Invalid rust elapsed time (<= 0): {rust_elapsed}")

    return {
        "name": run_dir.name,
        "rust_time_s": rust_elapsed,
        "python_time_s": py_elapsed,
        "speedup": py_elapsed / rust_elapsed,
        "rust_rss_mb": rust_time["max_rss_mb"],
        "python_rss_mb": py_time["max_rss_mb"],
        "rust_size_mb": parse_size_mb(rust_out),
        "python_size_mb": parse_size_mb(py_out),
    }


# ----------------------------
# Reporting
# ----------------------------

def format_error(e: BenchParseError) -> str:
    lines = []
    lines.append(str(e))

    if e.ctx is not None:
        ctx = e.ctx
        lines.append("")
        lines.append("Context:")
        lines.append(f"  file:   {ctx.path}")
        if ctx.lineno:
            lines.append(f"  line:   {ctx.lineno}")
        if ctx.rule:
            lines.append(f"  rule:   {ctx.rule}")
        if ctx.key is not None:
            lines.append(f"  key:    {ctx.key}")
        if ctx.value is not None:
            lines.append(f"  value:  {ctx.value}")
        if ctx.line:
            lines.append("  text:   " + ctx.line)
        if ctx.note:
            lines.append("")
            lines.append(ctx.note)

    if e.cause is not None:
        lines.append("")
        lines.append(f"Cause: {type(e.cause).__name__}: {e.cause}")

    return "\n".join(lines)


def print_markdown(rows: List[Dict[str, float]]) -> None:
    headers = [
        "run",
        "rust_time_s",
        "python_time_s",
        "speedup(py/rust)",
        "rust_rss_MB",
        "python_rss_MB",
        "rust_size_MB",
        "python_size_MB",
    ]
    print("| " + " | ".join(headers) + " |")
    print("|" + "|".join(["---"] * len(headers)) + "|")

    for r in rows:
        print(
            f"| {r['name']} "
            f"| {r['rust_time_s']:.2f} "
            f"| {r['python_time_s']:.2f} "
            f"| {r['speedup']:.2f} "
            f"| {r['rust_rss_mb']:.1f} "
            f"| {r['python_rss_mb']:.1f} "
            f"| {r['rust_size_mb']:.3f} "
            f"| {r['python_size_mb']:.3f} |"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description="Format bam-coverage vs deepTools benchmarks (strict + diagnostic)")
    ap.add_argument("dirs", nargs="+", help="benchmark output directories (e.g. bench_out)")
    args = ap.parse_args()

    rows = []
    try:
        for d in args.dirs:
            rows.append(load_run(Path(d).resolve()))
        print_markdown(rows)
    except BenchParseError as e:
        print("ERROR: benchmark parsing failed.\n", file=sys.stderr)
        print(format_error(e), file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

