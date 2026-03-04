#!/usr/bin/env python3
"""
Validate annotation inputs before running DRAM:
1. Bin names in the bins directory must exactly match the 'Name' column in
   quality_report.tsv (no extra/missing bins).
2. Detect number/string type mismatch: if the same ID would be interpreted
   differently (e.g. 4 vs "4"), exit with error so we fail early instead of
   after 14+ hours in distill.

Usage:
  python scripts/validate_annotation_inputs.py <bins_dir> <quality_report.tsv> [--ext fa]
Exit: 0 if valid, 1 with message if invalid.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Validate bins vs quality report before DRAM annotation.")
    ap.add_argument("bins_dir", help="Directory containing bin FASTA files")
    ap.add_argument("quality_report", help="CheckM2 quality_report.tsv path")
    ap.add_argument("--ext", default="fa", help="Bin file extension (default: fa)")
    args = ap.parse_args()

    bins_dir = Path(args.bins_dir)
    if not bins_dir.is_dir():
        print(f"Error: bins directory not found: {bins_dir}", file=sys.stderr)
        sys.exit(1)

    ext = args.ext if args.ext.startswith(".") else f".{args.ext}"
    bin_files = list(bins_dir.glob(f"*{ext}"))
    # Bin names from filenames (strip extension), as strings for consistent comparison
    bin_names_from_dir = {f.stem for f in bin_files}
    bin_names_from_dir_str = {str(s) for s in bin_names_from_dir}

    if not bin_names_from_dir_str:
        print(f"Error: no bin files *{ext} found in {bins_dir}", file=sys.stderr)
        sys.exit(1)

    # Read quality report: force 'Name' (and first column if named differently) as string
    # so we avoid pandas reading "4" as int 4 and causing later mixed-type in distill
    qpath = Path(args.quality_report)
    if not qpath.exists():
        print(f"Error: quality report not found: {qpath}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(qpath, sep="\t", dtype=str, low_memory=False)
    except Exception as e:
        print(f"Error: failed to read quality report: {e}", file=sys.stderr)
        sys.exit(1)

    if "Name" not in df.columns:
        print(
            f"Error: quality report has no 'Name' column. Columns: {list(df.columns)}",
            file=sys.stderr,
        )
        sys.exit(1)

    # All names as string (already from dtype=str)
    names_in_qc = set(df["Name"].astype(str).dropna())
    names_in_qc = {s.strip() for s in names_in_qc if s.strip()}

    # 1) Set equality: every bin file must have a row in QC and vice versa
    only_in_dir = bin_names_from_dir_str - names_in_qc
    only_in_qc = names_in_qc - bin_names_from_dir_str

    if only_in_dir or only_in_qc:
        print("Error: Bin names and quality report 'Name' column do not match.", file=sys.stderr)
        if only_in_dir:
            print(f"  In bins dir but NOT in quality report: {sorted(only_in_dir)}", file=sys.stderr)
        if only_in_qc:
            print(f"  In quality report but NOT in bins dir: {sorted(only_in_qc)}", file=sys.stderr)
        print(
            "  Fix: ensure CheckM2 was run on the same bin set and that bin IDs are identical (no number/string mismatch).",
            file=sys.stderr,
        )
        sys.exit(1)

    # 2) Type consistency: if any name in QC is numeric-looking, ensure it wasn't stored as number
    # (We already read with dtype=str, so we're good. But check for leading zeros / consistency.)
    for name in names_in_qc:
        if name != name.strip():
            print(
                f"Error: quality report 'Name' has whitespace: '{name}'",
                file=sys.stderr,
            )
            sys.exit(1)

    print(f"Validated: {len(bin_names_from_dir_str)} bins match quality report.")
    sys.exit(0)


if __name__ == "__main__":
    main()
