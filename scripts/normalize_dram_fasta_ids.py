#!/usr/bin/env python3
"""
Normalize DRAM TSV files for distill compatibility:
1) Prefix numeric-only values in the 'fasta' column to avoid mixed int/str IDs.
2) Strip trailing '.hmm' from 'cazy_best_hit' so DRAM distill function IDs match
   entries in function_heatmap_form (e.g., GH97.hmm -> GH97).

Usage:
  python scripts/normalize_dram_fasta_ids.py <input.tsv> <output.tsv> [--prefix g_]
"""

import argparse
import csv
import re
import sys
from pathlib import Path


NUMERIC_ONLY = re.compile(r"^\d+$")
HMM_SUFFIX = re.compile(r"\.hmm$", re.IGNORECASE)


def normalize_fasta_value(value: str, prefix: str) -> str:
    value = value.strip()
    if NUMERIC_ONLY.fullmatch(value):
        return f"{prefix}{value}"
    return value


def normalize_cazy_best_hit(value: str) -> str:
    value = value.strip()
    if not value:
        return value
    return HMM_SUFFIX.sub("", value)


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

    changed_fasta = 0
    changed_cazy = 0
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
                changed_fasta += 1
            row["fasta"] = new

            if "cazy_best_hit" in reader.fieldnames:
                old_cazy = row.get("cazy_best_hit", "")
                new_cazy = normalize_cazy_best_hit(old_cazy)
                if old_cazy != new_cazy:
                    changed_cazy += 1
                row["cazy_best_hit"] = new_cazy

            writer.writerow(row)

    print(
        "Wrote normalized DRAM TSV to "
        f"{out_path} (fasta rows updated: {changed_fasta}, "
        f"cazy_best_hit rows updated: {changed_cazy})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
