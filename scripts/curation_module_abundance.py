#!/usr/bin/env python3
"""
Standalone curated-module abundance for one contig-annotation subsample.

This is *separate* from aggregate_ko_module_abundance.py (the DRAM-native
module pipeline). It computes ONLY the SAP gene-curation modules
(carotenoids / glycan-LBP / LPS / polyphenol / bile acid / SCFA-48), so the
DRAM-native SCFA / polyphenol(CAMPER) / CAZy tables produced by the original
pipeline are left completely untouched -- you keep both versions of each.

It reuses the exact TPM + ID-extraction logic of the original pipeline by
importing the helper functions from aggregate_ko_module_abundance.py (that
module is imported read-only; nothing in it is modified).

Pipeline:
  featureCounts gene_counts.tsv + DRAM distill_input_fixed/annotations.tsv
      -> gene TPM (length + library-size normalised)
      -> ID TPM (KO/EC/Pfam/CAZy, fractional allocation)
      -> for each curated reference: per-feature + per-module abundance

Outputs (in --out-dir, one pair per curated reference):
  <name>_per_feature.tsv   one row per curated gene/enzyme + its TPM/detection
  <name>_per_module.tsv    one row per module: detection fraction, TPM rollups
  id_tpm.tsv               (optional shared intermediate, --write-id-tpm)

Reference TSVs are built in two ID-matching variants (see build_curation_references.py):
  id_only/   KO + CAZy + CAMPER D-id
  id_ec/     KO + EC + CAZy + CAMPER D-id
Pass the variant subfolder via --curation-refs-dir (e.g. scripts/curation_refs/id_only).

Usage:
  /path/to/dram/python scripts/curation_module_abundance.py \
    --counts       abundance/.../{prefix}/gene_counts.tsv \
    --annotations  annotation/.../{prefix}/distill_input_fixed/annotations.tsv \
    --curation-refs-dir scripts/curation_refs \
    --out-dir      annotation/.../{prefix}/custom_modules \
    --sample-label {prefix}

Run with the DRAM conda env python (mag_annotator must be importable).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Reuse the original pipeline's TPM / ID logic without modifying it.
from aggregate_ko_module_abundance import (
    compute_tpm,
    extract_ids_by_gene,
    gene_tpm_to_id_tpm,
    parse_featurecounts,
    split_ids_string,
)


def summarise_custom_module_abundance(ref_df: pd.DataFrame, id_tpm_map: dict):
    """
    Generic curated-module aggregation for the SAP gene-curation reference
    tables built by scripts/build_curation_references.py.

    Each reference row lists a feature's alternative IDs in an `ids` column
    (comma-separated, already normalised to id_tpm format: bare KO/CAZy/Pfam,
    'EC:'-prefixed EC). A gene's TPM was fractionally split across all of its
    IDs upstream, so summing a row's alternative IDs recovers that reaction's
    TPM without double counting.

    Returns (per_feature_df, per_module_df).
    """
    feat_rows = []
    for _, r in ref_df.iterrows():
        ids = sorted(split_ids_string(r.get("ids", "")))
        per_id = {i: id_tpm_map.get(i, 0.0) for i in ids}
        feature_tpm = float(sum(per_id.values()))
        hits = [i for i, t in per_id.items() if t > 0]
        feat_rows.append({
            "feature_class": r.get("feature_class", ""),
            "module": r.get("module", ""),
            "feature": r.get("feature", ""),
            "gene_symbol": r.get("gene_symbol", ""),
            "n_ids": len(ids),
            "n_ids_hit": len(hits),
            "hit_ids": ",".join(hits),
            "feature_tpm": feature_tpm,
            "feature_detected": feature_tpm > 0,
            "confidence": r.get("confidence", ""),
            "description": r.get("description", ""),
        })
    per_feature = pd.DataFrame(feat_rows)

    mod_rows = []
    if not per_feature.empty:
        for (fclass, module), g in per_feature.groupby(
            ["feature_class", "module"], sort=False
        ):
            n_feat = len(g)
            n_det = int(g["feature_detected"].sum())
            tpms = g["feature_tpm"]
            completeness = (n_det / n_feat) if n_feat else float("nan")
            mod_rows.append({
                "feature_class": fclass,
                "module": module,
                "n_features": n_feat,
                "n_features_detected": n_det,
                "detection_fraction": completeness,
                "module_sum_tpm": float(tpms.sum()),
                "module_mean_feature_tpm": float(tpms.mean()) if n_feat else float("nan"),
                "module_min_feature_tpm": float(tpms.min()) if n_feat else float("nan"),
                "module_weighted_tpm": float(tpms.sum() * completeness)
                if not (isinstance(completeness, float) and np.isnan(completeness))
                else float("nan"),
            })
    per_module = (
        pd.DataFrame(mod_rows)
        .sort_values(
            ["feature_class", "module_weighted_tpm", "module"],
            ascending=[True, False, True],
        )
        .reset_index(drop=True)
        if mod_rows
        else pd.DataFrame(mod_rows)
    )
    return per_feature, per_module


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--counts", required=True, help="featureCounts output TSV")
    ap.add_argument(
        "--annotations",
        required=True,
        help="DRAM distill_input_fixed/annotations.tsv (the normalised one)",
    )
    ap.add_argument(
        "--curation-ref",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Curated reference TSV. Repeatable, e.g. "
        "--curation-ref lbp=scripts/curation_refs/lbp.tsv.",
    )
    ap.add_argument(
        "--curation-refs-dir",
        default=None,
        help="Directory whose *.tsv are each loaded with NAME = filename stem.",
    )
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-label", default=None, help="Value for a 'sample' column")
    ap.add_argument(
        "--write-id-tpm",
        action="store_true",
        help="Also write the shared id_tpm.tsv intermediate into --out-dir.",
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    curation_specs = {}
    for entry in args.curation_ref:
        if "=" not in entry:
            raise SystemExit(f"ERROR: --curation-ref must be NAME=PATH, got: {entry}")
        name, path = entry.split("=", 1)
        curation_specs[name.strip()] = Path(path.strip())
    if args.curation_refs_dir:
        for p in sorted(Path(args.curation_refs_dir).glob("*.tsv")):
            curation_specs.setdefault(p.stem, p)
    if not curation_specs:
        raise SystemExit("ERROR: no curated references given (--curation-ref / --curation-refs-dir)")

    print("=== curation_module_abundance ===")
    print(f"counts:       {args.counts}")
    print(f"annotations:  {args.annotations}")
    print(f"out-dir:      {out_dir}")
    print(f"references:   {', '.join(sorted(curation_specs))}")

    # gene counts -> TPM
    gene_df = parse_featurecounts(Path(args.counts).resolve())
    gene_tpm = compute_tpm(gene_df)
    print(
        f"[TPM] {len(gene_tpm):,} genes; sum TPM = {gene_tpm['TPM'].sum():,.1f} (~1e6)"
    )

    # gene -> ID TPM
    annotations = pd.read_csv(args.annotations, sep="\t", index_col=0)
    gene_to_ids = extract_ids_by_gene(annotations)
    id_tpm_df = gene_tpm_to_id_tpm(gene_tpm, gene_to_ids)
    id_tpm_map = dict(zip(id_tpm_df["id"], id_tpm_df["tpm"]))
    print(f"[ID TPM] {len(id_tpm_df):,} unique IDs")
    if args.write_id_tpm:
        out = id_tpm_df.copy()
        if args.sample_label:
            out.insert(0, "sample", args.sample_label)
        out.to_csv(out_dir / "id_tpm.tsv", sep="\t", index=False)

    for name, ref_path in curation_specs.items():
        ref_df = pd.read_csv(ref_path, sep="\t", dtype=str).fillna("")
        per_feature, per_module = summarise_custom_module_abundance(ref_df, id_tpm_map)
        for df in (per_feature, per_module):
            if args.sample_label and not df.empty:
                df.insert(0, "sample", args.sample_label)
        per_feature.to_csv(out_dir / f"{name}_per_feature.tsv", sep="\t", index=False)
        per_module.to_csv(out_dir / f"{name}_per_module.tsv", sep="\t", index=False)
        n_det = int(per_feature["feature_detected"].sum()) if not per_feature.empty else 0
        print(
            f"[{name}] {len(per_feature)} features ({n_det} detected) "
            f"across {len(per_module)} modules"
        )

    print(f"\nWrote curated-module outputs under {out_dir}")


if __name__ == "__main__":
    main()
