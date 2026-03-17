#!/usr/bin/env python3
"""
Scan FragGeneScanRs metadata across multiple score thresholds.

This script reads only the metadata file (-m output), so threshold scanning is
fast and does not duplicate huge FAA filtering work.

It also reports score distribution summary statistics:
  - min / max / mean / std
  - q25 / median / q75 (estimated from a bounded reservoir sample)

Usage:
  python scripts/scan_fraggenescanrs_scores.py \
    --meta-in genes_predicted.meta.tsv \
    --summary-out genes_predicted.score_scan.tsv \
    --thresholds 1.20,1.25,1.30,1.35
"""

import argparse
import math
import random
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute FragGeneScanRs keep/remove stats for multiple score thresholds."
    )
    parser.add_argument("--meta-in", required=True, help="FragGeneScanRs metadata file (-m)")
    parser.add_argument("--summary-out", required=True, help="Output TSV summary path")
    parser.add_argument(
        "--thresholds",
        required=True,
        help="Comma-separated score thresholds, e.g. 1.20,1.25,1.30,1.35",
    )
    parser.add_argument(
        "--quantile-sample-size",
        type=int,
        default=2_000_000,
        help="Maximum sample size for quantile estimation (default: 2000000)",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=11,
        help="Random seed for quantile sampling (default: 11)",
    )
    return parser.parse_args()


def iter_scores(meta_path: Path):
    with meta_path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                raise ValueError(f"Malformed metadata row: {line}")
            yield float(parts[4].strip())


def percentile_from_sorted(values, p):
    if not values:
        return float("nan")
    if len(values) == 1:
        return values[0]
    pos = (len(values) - 1) * p
    low = int(math.floor(pos))
    high = int(math.ceil(pos))
    if low == high:
        return values[low]
    weight = pos - low
    return values[low] * (1 - weight) + values[high] * weight


def main():
    args = parse_args()
    meta_in = Path(args.meta_in)
    summary_out = Path(args.summary_out)
    sample_limit = max(1, int(args.quantile_sample_size))
    random.seed(int(args.random_seed))

    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]
    if not thresholds:
        raise ValueError("No thresholds parsed from --thresholds")

    # Keep original order as provided by user/config.
    kept_counts = [0 for _ in thresholds]
    total = 0
    score_min = None
    score_max = None
    mean = 0.0
    m2 = 0.0
    sample = []

    for score in iter_scores(meta_in):
        total += 1
        if score_min is None or score < score_min:
            score_min = score
        if score_max is None or score > score_max:
            score_max = score

        # Welford one-pass mean/variance
        delta = score - mean
        mean += delta / total
        delta2 = score - mean
        m2 += delta * delta2

        # Reservoir sampling for scalable quantile estimation
        if len(sample) < sample_limit:
            sample.append(score)
        else:
            j = random.randint(1, total)
            if j <= sample_limit:
                sample[j - 1] = score

        for i, t in enumerate(thresholds):
            if score >= t:
                kept_counts[i] += 1

    if total > 1:
        variance = m2 / (total - 1)
    else:
        variance = 0.0
    std = math.sqrt(max(variance, 0.0))

    sample.sort()
    q25 = percentile_from_sorted(sample, 0.25)
    q50 = percentile_from_sorted(sample, 0.50)
    q75 = percentile_from_sorted(sample, 0.75)

    summary_out.parent.mkdir(parents=True, exist_ok=True)
    with summary_out.open("w") as out:
        out.write(
            "threshold\ttotal_predictions\tkept_predictions\tremoved_predictions\tkeep_fraction\t"
            "score_min\tscore_mean\tscore_std\tscore_max\tq25_est\tmedian_est\tq75_est\tquantile_sample_size\n"
        )
        for t, kept in zip(thresholds, kept_counts):
            removed = total - kept
            keep_fraction = (kept / total) if total else 0.0
            out.write(
                f"{t}\t{total}\t{kept}\t{removed}\t{keep_fraction:.6f}\t"
                f"{score_min:.6f}\t{mean:.6f}\t{std:.6f}\t{score_max:.6f}\t"
                f"{q25:.6f}\t{q50:.6f}\t{q75:.6f}\t{len(sample)}\n"
            )


if __name__ == "__main__":
    main()
