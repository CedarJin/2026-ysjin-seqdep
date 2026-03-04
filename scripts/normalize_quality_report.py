#!/usr/bin/env python3
"""
Write a copy of CheckM2 quality_report.tsv with the 'Name' column forced to string
so DRAM and downstream distill do not see mixed types (e.g. 4 vs "4") and crash.

Usage:
  python scripts/normalize_quality_report.py <quality_report.tsv> <output.tsv>
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Normalize quality report Name column to string for DRAM.")
    ap.add_argument("input_tsv", help="CheckM2 quality_report.tsv")
    ap.add_argument("output_tsv", help="Output path for normalized TSV")
    args = ap.parse_args()

    p = Path(args.input_tsv)
    if not p.exists():
        print(f"Error: input not found: {p}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(p, sep="\t", low_memory=False)
    if "Name" not in df.columns:
        print(f"Error: no 'Name' column in {p}", file=sys.stderr)
        sys.exit(1)

    # Force Name to string so downstream (e.g. Altair in distill) never sees int/str mix
    df["Name"] = df["Name"].astype(str)
    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False)
    print(f"Wrote normalized quality report to {args.output_tsv}")
    sys.exit(0)


if __name__ == "__main__":
    main()
