#!/usr/bin/env python3
"""Merge featureCounts output with DRAM annotations for one min-contig threshold."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", required=True, help="DRAM annotations.tsv")
    parser.add_argument("--counts", required=True, help="featureCounts output table")
    parser.add_argument("--min-contig-len", required=True, type=int, help="Current min contig length")
    parser.add_argument("--output", required=True, help="Merged abundance table")
    return parser.parse_args()


def read_annotations(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header in {path}")

        ann_fields = ["gene_id" if name == "" else name for name in reader.fieldnames]
        rows: dict[str, dict[str, str]] = {}

        for row in reader:
            normalized = {}
            for original_name, new_name in zip(reader.fieldnames, ann_fields):
                normalized[new_name] = row.get(original_name, "")

            gene_id = normalized["gene_id"]
            if gene_id:
                rows[gene_id] = normalized

    return ann_fields, rows


def read_featurecounts(path: Path) -> tuple[str, dict[str, dict[str, str]]]:
    header = None
    rows: dict[str, dict[str, str]] = {}

    with path.open() as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue

            line = raw_line.rstrip("\n")
            if not line:
                continue

            if header is None:
                header = line.split("\t")
                continue

            values = line.split("\t")
            row = dict(zip(header, values))
            gene_id = row["Geneid"]
            rows[gene_id] = row

    if header is None:
        raise ValueError(f"Could not find featureCounts header in {path}")

    sample_column = header[-1]
    return sample_column, rows


def main() -> None:
    args = parse_args()
    ann_fields, annotations = read_annotations(Path(args.annotations))
    sample_column, counts = read_featurecounts(Path(args.counts))

    output_fields = (
        ["gene_id", "min_contig_len"]
        + [field for field in ann_fields if field != "gene_id"]
        + ["fc_chr", "fc_start", "fc_end", "fc_strand", "fc_length_bp", "assigned_fragments"]
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for gene_id in sorted(counts):
            row = {"gene_id": gene_id, "min_contig_len": str(args.min_contig_len)}
            row.update(annotations.get(gene_id, {}))

            fc_row = counts[gene_id]
            row.update(
                {
                    "fc_chr": fc_row.get("Chr", ""),
                    "fc_start": fc_row.get("Start", ""),
                    "fc_end": fc_row.get("End", ""),
                    "fc_strand": fc_row.get("Strand", ""),
                    "fc_length_bp": fc_row.get("Length", ""),
                    "assigned_fragments": fc_row.get(sample_column, "0"),
                }
            )
            writer.writerow(row)

    print(f"Wrote merged abundance table to {out_path}")


if __name__ == "__main__":
    main()
