#!/usr/bin/env python3
import sys
import math
from collections import defaultdict

def parse_bedgraph(path):
    """
    Yields (chr, start, end, value) in file order.
    bedGraph format: chr \t start \t end \t value
    """
    with open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                raise ValueError(f"{path}: bad line (expected 4 cols): {line!r}")
            chr_, s, e, v = parts
            yield chr_, int(s), int(e), float(v)

def chrom_order_key(chrname):
    # grouping only; not trying to impose chr9/chr10 semantics
    return chrname

def group_by_chr(records):
    d = defaultdict(list)
    for rec in records:
        d[rec[0]].append(rec)
    return d

def intervals_to_bins(intervals, binw, chrom_size=None):
    """
    Convert (start,end,value) intervals into dict bin_index->value.
    Missing bins imply 0.0.
    """
    bins = {}
    for (_chr, start, end, val) in intervals:
        if end <= start:
            continue
        b0 = start // binw
        b1 = (end + binw - 1) // binw  # ceil
        for b in range(b0, b1):
            bins[b] = val

    if chrom_size is not None:
        maxbin = (chrom_size + binw - 1) // binw
        for b in list(bins.keys()):
            if b >= maxbin:
                del bins[b]
    return bins

def main():
    if len(sys.argv) < 4:
        print("usage: compare_bedgraph_bins.py <rust.bedgraph> <python.bedgraph> <bin_width> [eps]", file=sys.stderr)
        sys.exit(2)

    rust_path, py_path = sys.argv[1], sys.argv[2]
    binw = int(sys.argv[3])
    eps = float(sys.argv[4]) if len(sys.argv) >= 5 else 1e-6

    rust_by_chr = group_by_chr(parse_bedgraph(rust_path))
    py_by_chr = group_by_chr(parse_bedgraph(py_path))
    all_chrs = sorted(set(rust_by_chr.keys()) | set(py_by_chr.keys()), key=chrom_order_key)

    total_bins = 0
    bins_with_signal = 0

    n_bad = 0
    sum_abs = 0.0
    max_abs = 0.0

    # relative error stats (only where max(x,y) > 0)
    sum_rel = 0.0
    max_rel = 0.0
    n_rel = 0

    # Pearson correlation accumulators (over all bins)
    n = 0
    sx = sy = sxx = syy = sxy = 0.0

    for chr_ in all_chrs:
        rints = rust_by_chr.get(chr_, [])
        pints = py_by_chr.get(chr_, [])

        rbins = intervals_to_bins(rints, binw)
        pbins = intervals_to_bins(pints, binw)

        chr_bins = set(rbins.keys()) | set(pbins.keys())
        for b in chr_bins:
            x = rbins.get(b, 0.0)
            y = pbins.get(b, 0.0)

            d = abs(x - y)
            total_bins += 1

            sum_abs += d
            if d > max_abs:
                max_abs = d
            if d > eps:
                n_bad += 1

            denom_rel = max(x, y)
            if denom_rel > 0:
                bins_with_signal += 1
                rel = d / denom_rel
                sum_rel += rel
                n_rel += 1
                if rel > max_rel:
                    max_rel = rel

            # correlation (all bins)
            n += 1
            sx += x
            sy += y
            sxx += x * x
            syy += y * y
            sxy += x * y

    mean_abs = sum_abs / total_bins if total_bins else 0.0
    frac_bad = n_bad / total_bins if total_bins else 0.0

    mean_rel = (sum_rel / n_rel) if n_rel else 0.0

    denom = math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
    corr = (n * sxy - sx * sy) / denom if denom != 0.0 else float("nan")

    print("== bedGraph binwise comparison ==")
    print(f"bin_width\t{binw}")
    print(f"eps\t\t{eps}")
    print(f"total_bins\t{total_bins}")
    print(f"bins_with_signal\t{bins_with_signal}")
    print(f"bad_bins\t{n_bad}")
    print(f"frac_bad\t{frac_bad:.6g}")
    print(f"mean_abs\t{mean_abs:.6g}")
    print(f"max_abs\t\t{max_abs:.6g}")
    print(f"mean_rel\t{mean_rel:.6g}")
    print(f"max_rel\t\t{max_rel:.6g}")
    print(f"corr\t\t{corr:.6g}")

if __name__ == "__main__":
    main()
