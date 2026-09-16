#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import pandas as pd


def read_barcodes(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return {
            line.strip().split("\t")[0].removesuffix("-1")
            for line in fh
            if line.strip()
        }


def find_sample(path: Path, known_samples: set[str]) -> str | None:
    # Prefer exact directory-component matches.
    for part in reversed(path.parts):
        if part in known_samples:
            return part

    # Valkyrn compare currently uses shorter sample labels.
    aliases = {
        "ZD-4631-BcellsLaneF": "Bcells",
        "ZD-4631-PCLaneG": "Plasma",
        "ZD-4631-Undetermined": "Undetermined",
    }

    text = str(path)
    for directory_name, sample in aliases.items():
        if directory_name in text and sample in known_samples:
            return sample

    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--results",
        required=True,
        help="Norn results directory",
    )
    ap.add_argument(
        "--annotations",
        required=True,
        help="Valkyrn compare cell_annotations.tsv",
    )
    ap.add_argument(
        "--out",
        default="barcode_vs_valkyrn",
    )
    args = ap.parse_args()

    results = Path(args.results)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    ann = pd.read_csv(args.annotations, sep="\t", dtype=str).fillna("")

    required = {"sample", "cell"}
    missing = required - set(ann.columns)
    if missing:
        raise SystemExit(
            f"Missing required cell_annotations columns: {sorted(missing)}"
        )

    # Normalise only the conventional 10x suffix.
    ann["cell_norm"] = ann["cell"].str.removesuffix("-1")

    samples = set(ann["sample"])

    barcode_files = sorted(
        list(results.rglob("barcodes.tsv.gz"))
        + list(results.rglob("barcodes.tsv"))
    )

    if not barcode_files:
        raise SystemExit(f"No barcodes.tsv[.gz] found under {results}")

    summary = []
    classifications = []

    print(f"Found {len(barcode_files)} barcode files\n")

    for barcode_file in barcode_files:
        sample = find_sample(barcode_file, samples)

        if sample is None:
            print(f"SKIP  {barcode_file}")
            print("      cannot associate with a Valkyrn sample")
            continue

        gex = read_barcodes(barcode_file)

        sample_ann = ann[ann["sample"] == sample]
        vdj = set(sample_ann["cell_norm"])

        both = gex & vdj
        gex_only = gex - vdj
        vdj_only = vdj - gex

        summary.append({
            "sample": sample,
            "barcode_file": str(barcode_file),
            "gex_barcodes": len(gex),
            "valkyrn_cells": len(vdj),
            "gex_and_vdj": len(both),
            "gex_only": len(gex_only),
            "vdj_only": len(vdj_only),
            "pct_gex_with_vdj":
                100.0 * len(both) / len(gex) if gex else 0,
            "pct_vdj_in_gex":
                100.0 * len(both) / len(vdj) if vdj else 0,
        })

        for cell in sorted(gex | vdj):
            if cell in both:
                status = "GEX+VDJ"
            elif cell in gex:
                status = "GEX_ONLY"
            else:
                status = "VDJ_ONLY"

            classifications.append({
                "sample": sample,
                "cell": cell,
                "status": status,
                "barcode_file": str(barcode_file),
            })

    summary_df = pd.DataFrame(summary)
    cells_df = pd.DataFrame(classifications)

    summary_df.to_csv(
        out / "barcode_vs_valkyrn_summary.tsv",
        sep="\t",
        index=False,
    )

    cells_df.to_csv(
        out / "barcode_vs_valkyrn_cells.tsv",
        sep="\t",
        index=False,
    )

    if summary_df.empty:
        raise SystemExit("No barcode files could be associated with samples")

    display = summary_df[
        [
            "sample",
            "gex_barcodes",
            "valkyrn_cells",
            "gex_and_vdj",
            "gex_only",
            "vdj_only",
            "pct_gex_with_vdj",
            "pct_vdj_in_gex",
        ]
    ].copy()

    display["pct_gex_with_vdj"] = display["pct_gex_with_vdj"].map(
        lambda x: f"{x:.1f}%"
    )
    display["pct_vdj_in_gex"] = display["pct_vdj_in_gex"].map(
        lambda x: f"{x:.1f}%"
    )

    print(display.to_string(index=False))
    print()
    print(f"Wrote {out / 'barcode_vs_valkyrn_summary.tsv'}")
    print(f"Wrote {out / 'barcode_vs_valkyrn_cells.tsv'}")


if __name__ == "__main__":
    main()
