#!/usr/bin/env python3
"""
Normalize DRAM TSV files so numeric-only values in the 'fasta' column are
prefixed, preventing mixed int/str handling in downstream distill plotting.

Usage:
  python scripts/normalize_dram_fasta_ids.py <input.tsv> <output.tsv> [--prefix g_]
"""

import argparse
import csv
import re
import sys
from pathlib import Path


NUMERIC_ONLY = re.compile(r"^\d+$")


def normalize_fasta_value(value: str, prefix: str) -> str:
    value = value.strip()
    if NUMERIC_ONLY.fullmatch(value):
        return f"{prefix}{value}"
    return value


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Prefix numeric-only fasta IDs in a DRAM TSV file."
    )
    ap.add_argument("input_tsv", help="Input TSV path")
    ap.add_argument("output_tsv", help="Output TSV path")
    ap.add_argument(
        "--prefix",
        default="bin.",
        help="Prefix for numeric-only fasta IDs (default: bin.)",
    )
    args = ap.parse_args()

    in_path = Path(args.input_tsv)
    out_path = Path(args.output_tsv)
    if not in_path.exists():
        print(f"Error: input not found: {in_path}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)

    changed = 0
    with in_path.open("r", newline="") as fin, out_path.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            print(f"Error: empty TSV: {in_path}", file=sys.stderr)
            return 1
        if "fasta" not in reader.fieldnames:
            print(f"Error: no 'fasta' column in {in_path}", file=sys.stderr)
            return 1

        writer = csv.DictWriter(
            fout, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for row in reader:
            old = row.get("fasta", "")
            new = normalize_fasta_value(old, args.prefix)
            if old != new:
                changed += 1
            row["fasta"] = new
            writer.writerow(row)

    print(f"Wrote normalized DRAM TSV to {out_path} (updated {changed} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
