#!/usr/bin/env python3
"""
Filter FragGeneScanRs predicted proteins by metadata score.

This script streams metadata and FAA in lockstep to avoid loading huge files into memory.

Usage:
  python scripts/filter_fraggenescanrs_by_score.py \
    --faa-in genes_raw.faa \
    --meta-in genes.meta \
    --faa-out genes_filtered.faa \
    --stats-out filter_stats.tsv \
    --min-score 1.30
"""

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter FragGeneScanRs FAA by metadata prediction score."
    )
    parser.add_argument("--faa-in", required=True, help="Input FAA from FragGeneScanRs (-a)")
    parser.add_argument("--meta-in", required=True, help="Input metadata from FragGeneScanRs (-m)")
    parser.add_argument("--faa-out", required=True, help="Output FAA after score filtering")
    parser.add_argument("--stats-out", required=True, help="Output TSV with filtering stats")
    parser.add_argument("--min-score", type=float, required=True, help="Keep predictions with score >= min-score")
    return parser.parse_args()


def iter_meta(meta_path: Path):
    current_read = None
    with meta_path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_read = line[1:].strip()
                continue
            if current_read is None:
                raise ValueError(f"Malformed metadata line before any read header: {line}")
            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(f"Malformed metadata row: {line}")
            start = parts[0].strip()
            end = parts[1].strip()
            strand = parts[2].strip()
            score = float(parts[4].strip())
            gene_id = f"{current_read}_{start}_{end}_{strand}"
            yield gene_id, score


def iter_faa_records(faa_path: Path):
    header = None
    seq_lines = []
    with faa_path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line[1:].strip()
                seq_lines = []
            else:
                if header is None:
                    raise ValueError("Malformed FAA: sequence line found before header")
                seq_lines.append(line)
        if header is not None:
            yield header, seq_lines


def main():
    args = parse_args()
    faa_in = Path(args.faa_in)
    meta_in = Path(args.meta_in)
    faa_out = Path(args.faa_out)
    stats_out = Path(args.stats_out)
    min_score = float(args.min_score)

    faa_out.parent.mkdir(parents=True, exist_ok=True)
    stats_out.parent.mkdir(parents=True, exist_ok=True)

    meta_iter = iter_meta(meta_in)
    total = 0
    kept = 0

    with faa_out.open("w") as out_handle:
        for faa_header, seq_lines in iter_faa_records(faa_in):
            total += 1
            try:
                meta_id, score = next(meta_iter)
            except StopIteration as exc:
                raise ValueError(
                    "Metadata ended before FAA; input files are not aligned from the same FragGeneScanRs run."
                ) from exc

            if faa_header != meta_id:
                raise ValueError(
                    f"Header mismatch between FAA and metadata: faa='{faa_header}' meta='{meta_id}'"
                )

            if score >= min_score:
                kept += 1
                out_handle.write(f">{faa_header}\n")
                out_handle.write("\n".join(seq_lines))
                out_handle.write("\n")

        try:
            extra_meta = next(meta_iter)
            raise ValueError(
                f"Metadata has extra entries beyond FAA (first extra id: {extra_meta[0]})."
            )
        except StopIteration:
            pass

    removed = total - kept
    keep_fraction = (kept / total) if total else 0.0

    with stats_out.open("w") as stats:
        stats.write("metric\tvalue\n")
        stats.write(f"min_score\t{min_score}\n")
        stats.write(f"total_predictions\t{total}\n")
        stats.write(f"kept_predictions\t{kept}\n")
        stats.write(f"removed_predictions\t{removed}\n")
        stats.write(f"keep_fraction\t{keep_fraction:.6f}\n")


if __name__ == "__main__":
    main()
