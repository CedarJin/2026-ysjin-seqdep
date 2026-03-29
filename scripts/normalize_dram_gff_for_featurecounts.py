#!/usr/bin/env python3
"""Normalize DRAM GFF seqids so featureCounts matches BAM reference names."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_gff", help="Original DRAM genes.gff")
    parser.add_argument("output_gff", help="Normalized GFF for featureCounts")
    return parser.parse_args()


def normalize_seqid(seqid: str) -> str:
    prefix = "final.contigs_"
    if seqid.startswith(prefix):
        return seqid[len(prefix) :]
    return seqid


def main() -> None:
    args = parse_args()
    in_path = Path(args.input_gff)
    out_path = Path(args.output_gff)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wrote_header = False
    records = 0

    with in_path.open() as inp, out_path.open("w") as out:
        for raw_line in inp:
            if raw_line.startswith("#"):
                if not wrote_header:
                    out.write("##gff-version 3\n")
                    wrote_header = True
                continue

            line = raw_line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected 9 GFF columns in {in_path}, got {len(fields)}: {line}")

            fields[0] = normalize_seqid(fields[0])
            out.write("\t".join(fields) + "\n")
            records += 1

    if records == 0:
        raise ValueError(f"No feature records were written from {in_path}")

    print(f"Wrote {records} normalized CDS records to {out_path}")


if __name__ == "__main__":
    main()
