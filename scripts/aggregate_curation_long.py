#!/usr/bin/env python3
"""
Concatenate per-subsample curated *_per_module.tsv tables into one long table
for the SAP sequencing-depth benchmarking analysis.

Each input path is expected to look like:
  <...>/{omic}/{sample}/{prefix}/custom_modules/{id_variant}/{name}_per_module.tsv
where {prefix} = {sample}_{depth}_seed{seed}, e.g. MG0001_50M_seed11.

The output adds identifying columns parsed from the path so the table can be
grouped by depth / seed / omic / sample / id_variant for the depth-benchmarking
models:
  omic  sample  subsample  depth  depth_reads  seed  id_variant  <original module columns...>

`depth_reads` is the numeric read-pair count (e.g. 50M -> 50000000) for easy
sorting / plotting on a continuous axis.

Usage (from Snakemake, passing the per-module files explicitly):
  python scripts/aggregate_curation_long.py --out <long.tsv> <file1> <file2> ...

Or discover them under a root:
  python scripts/aggregate_curation_long.py --out <long.tsv> \
    --root annotation/dram_contigs_T2T --subdir custom_modules
"""
import argparse
import re
from pathlib import Path

import pandas as pd

PREFIX_RE = re.compile(r"^(?P<sample>.+)_(?P<depth>\d+M)_seed(?P<seed>\d+)$")
DEPTH_UNIT = {"M": 1_000_000, "K": 1_000, "G": 1_000_000_000}


def depth_to_reads(depth: str) -> int:
    m = re.fullmatch(r"(\d+)([KMG])", depth)
    if not m:
        return -1
    return int(m.group(1)) * DEPTH_UNIT[m.group(2)]


def parse_path(path: Path):
    """Return (omic, sample, subsample, name, id_variant) from a per_module path."""
    parts = path.parts
    # .../{omic}/{sample}/{prefix}/{subdir}/{id_variant}/{name}_per_module.tsv
    name = path.name.replace("_per_module.tsv", "")
    if parts[-2] in ("id_only", "id_ec"):
        id_variant = parts[-2]
        subsample = parts[-4]
        sample = parts[-5]
        omic = parts[-6]
    else:
        id_variant = ""
        subsample = path.parent.parent.name
        sample = path.parent.parent.parent.name
        omic = path.parent.parent.parent.parent.name
    return omic, sample, subsample, name, id_variant


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("inputs", nargs="*", help="per_module.tsv files")
    ap.add_argument("--root", default=None, help="Root to glob for per_module files")
    ap.add_argument("--subdir", default="custom_modules", help="Per-subsample subfolder name")
    ap.add_argument("--out", required=True, help="Output long TSV")
    args = ap.parse_args()

    paths = [Path(p) for p in args.inputs]
    if args.root:
        paths += sorted(
            Path(args.root).glob(f"*/*/*/{args.subdir}/*_per_module.tsv")
        )
    # de-dup while preserving order
    seen, uniq = set(), []
    for p in paths:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            uniq.append(p)
    paths = uniq
    if not paths:
        raise SystemExit("ERROR: no per_module.tsv inputs found")

    frames = []
    for p in paths:
        if not p.is_file():
            print(f"WARN: skipping missing {p}")
            continue
        df = pd.read_csv(p, sep="\t")
        if df.empty:
            continue
        omic, sample, subsample, name, id_variant = parse_path(p)
        if not id_variant:
            id_variant = "id_only"
        m = PREFIX_RE.match(subsample)
        depth = m.group("depth") if m else ""
        seed = m.group("seed") if m else ""
        # drop any pre-existing 'sample' label column to avoid clashes
        if "sample" in df.columns:
            df = df.drop(columns=["sample"])
        meta = {
            "omic": omic,
            "sample": sample,
            "subsample": subsample,
            "depth": depth,
            "depth_reads": depth_to_reads(depth) if depth else -1,
            "seed": seed,
            "id_variant": id_variant,
        }
        for k, v in reversed(list(meta.items())):
            df.insert(0, k, v)
        frames.append(df)

    if not frames:
        raise SystemExit("ERROR: all inputs were empty")

    long_df = pd.concat(frames, ignore_index=True)
    long_df = long_df.sort_values(
        ["feature_class", "module", "id_variant", "omic", "sample", "depth_reads", "seed"],
        kind="stable",
    ).reset_index(drop=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(out, sep="\t", index=False)
    print(
        f"Wrote {len(long_df):,} rows from {len(frames)} subsample tables -> {out}\n"
        f"  feature_classes: {sorted(long_df['feature_class'].unique())}\n"
        f"  id_variants: {sorted(long_df['id_variant'].unique())}\n"
        f"  depths: {sorted(long_df['depth'].unique())}  seeds: {sorted(long_df['seed'].unique())}\n"
        f"  omics: {sorted(long_df['omic'].unique())}  n_subsamples: {long_df['subsample'].nunique()}"
    )


if __name__ == "__main__":
    main()
