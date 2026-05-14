#!/usr/bin/env python3
"""
Extra SCFA-production annotation for a DRAM contig-annotation sample.

DRAM's distill product.tsv only checks a sparse subset of SCFA marker
genes. In particular, it misses (a) the gut-dominant
butyryl-CoA:acetate CoA-transferase ("But") pathway for butyrate, and
(b) the succinate/methylmalonyl-CoA mutase pathway for propionate that
is the major route in Bacteroidetes. This script uses a curated SCFA
reference table (scripts/scfa_reference.tsv) to re-annotate a sample's
SCFA-production potential using DRAM's own ID extraction logic, then
prints a diff against DRAM's product.tsv SCFA section.

Usage:
  /path/to/dram/python scripts/extra_scfa_check.py \\
    --sample-dir annotation/dram_contigs_T2T/metaG/MG0001/MG0001_2M_seed11 \\
    --reference scripts/scfa_reference.tsv \\
    --out-dir   annotation/dram_contigs_T2T/metaG/MG0001/MG0001_2M_seed11/scfa_extra

Run this with the DRAM conda env's python (mag_annotator must be importable).
"""

import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd


# DRAM's SCFA function names in product.tsv (the columns we want to diff
# against). We map our pathway-level summary to these where applicable.
DRAM_SCFA_COLUMNS = [
    "SCFA and alcohol conversions: pyruvate => acetyl CoA v1",
    "SCFA and alcohol conversions: pyruvate => acetyl CoA v2",
    "SCFA and alcohol conversions: pyruvate => acetylCoA f+ formate v3",
    "SCFA and alcohol conversions: acetate, pt 1",
    "SCFA and alcohol conversions: acetate, pt 2",
    "SCFA and alcohol conversions: acetate, pt 3",
    "SCFA and alcohol conversions: lactate L",
    "SCFA and alcohol conversions: lactate D",
    "SCFA and alcohol conversions: Butyrate, pt 1",
    "SCFA and alcohol conversions: Butyrate, pt 2",
    "SCFA and alcohol conversions: Propionate, pt 1",
    "SCFA and alcohol conversions: Propionate, pt 2",
    "SCFA and alcohol conversions: Alcohol production",
]


def split_ids(value: str):
    """Split a comma/space-separated ID string into a set of stripped IDs."""
    if value is None:
        return set()
    parts = re.split(r"[,;\s]+", str(value).strip())
    return {p for p in parts if p}


def build_observed_id_counter(annotations: pd.DataFrame) -> Counter:
    """
    Build a Counter mapping ID -> gene-count using DRAM's own extractor so
    that we match KO/EC/Pfam/CAZy IDs exactly the way DRAM does.
    """
    from mag_annotator.summarize_genomes import get_ids_from_annotations_all
    return get_ids_from_annotations_all(annotations)


def annotate_reference(
    reference: pd.DataFrame,
    observed_counter: Counter,
) -> pd.DataFrame:
    """
    For every reference row, decide which of its IDs are present in the
    sample, count hit genes, and record the leftover IDs.
    """
    observed_set = set(observed_counter.keys())
    out_rows = []
    for _, row in reference.iterrows():
        ids = sorted(split_ids(row["ids"]))
        hits = sorted(set(ids) & observed_set)
        miss = sorted(set(ids) - observed_set)
        hit_count = sum(observed_counter[k] for k in hits)
        out_rows.append({
            **row.to_dict(),
            "n_ids": len(ids),
            "n_ids_hit": len(hits),
            "hit_gene_count": hit_count,
            "hit_ids": ",".join(hits) if hits else "",
            "missing_ids": ",".join(miss) if miss else "",
            "step_present": bool(hits),
        })
    return pd.DataFrame(out_rows)


def summarise_pathways(per_step: pd.DataFrame) -> pd.DataFrame:
    """
    For each (scfa, pathway), compute essential-step completeness:
      complete = all essential steps present
      essential_steps_present / essential_steps_total
    Non-essential rows (essential != 'yes') are counted as auxiliary.
    """
    rows = []
    for (scfa, pathway), grp in per_step.groupby(["scfa", "pathway"], sort=False):
        essential = grp[grp["essential"].str.lower() == "yes"]
        n_ess = len(essential)
        n_ess_present = int(essential["step_present"].sum())
        aux = grp[grp["essential"].str.lower() != "yes"]
        n_aux = len(aux)
        n_aux_present = int(aux["step_present"].sum())
        rows.append({
            "scfa": scfa,
            "pathway": pathway,
            "essential_steps_total": n_ess,
            "essential_steps_present": n_ess_present,
            "essential_completeness": (n_ess_present / n_ess) if n_ess else float("nan"),
            "pathway_complete": (n_ess_present == n_ess and n_ess > 0),
            "aux_steps_total": n_aux,
            "aux_steps_present": n_aux_present,
        })
    return pd.DataFrame(rows)


def summarise_scfa(pathway_summary: pd.DataFrame) -> pd.DataFrame:
    """
    For each SCFA: does at least one pathway complete the essential steps?
    Also report 'partial' (>=1 essential step present in some pathway).
    """
    rows = []
    for scfa, grp in pathway_summary.groupby("scfa", sort=False):
        rows.append({
            "scfa": scfa,
            "n_pathways_checked": len(grp),
            "n_pathways_complete": int(grp["pathway_complete"].sum()),
            "n_pathways_partial": int((grp["essential_steps_present"] > 0).sum()),
            "max_essential_completeness": float(grp["essential_completeness"].max()),
            "scfa_producible": bool(grp["pathway_complete"].any()),
            "scfa_partial": bool((grp["essential_steps_present"] > 0).any()),
        })
    return pd.DataFrame(rows)


def compare_to_dram(
    scfa_summary: pd.DataFrame,
    pathway_summary: pd.DataFrame,
    product_tsv: Path,
) -> pd.DataFrame:
    """Build a row-per-DRAM-column diff vs our scfa annotation."""
    product = pd.read_csv(product_tsv, sep="\t", index_col=0)
    # we expect exactly one genome row for contigs distill
    genome = product.index[0]

    # Manual mapping from DRAM product columns -> our scfa+pathway concept
    mapping = {
        "SCFA and alcohol conversions: pyruvate => acetyl CoA v1": ("pyruvate->acetyl-CoA", "PDH (v1)"),
        "SCFA and alcohol conversions: pyruvate => acetyl CoA v2": ("pyruvate->acetyl-CoA", "POR/OFOR (v2)"),
        "SCFA and alcohol conversions: pyruvate => acetylCoA f+ formate v3": ("pyruvate->acetyl-CoA", "PFL (v3)"),
        "SCFA and alcohol conversions: acetate, pt 1": ("acetate", "PTA-ACK"),
        "SCFA and alcohol conversions: acetate, pt 2": ("acetate", "Acd (ADP-forming)"),
        "SCFA and alcohol conversions: acetate, pt 3": ("acetate", "Acetyl-CoA hydrolase"),
        "SCFA and alcohol conversions: lactate L": ("lactate", "L-lactate from pyruvate"),
        "SCFA and alcohol conversions: lactate D": ("lactate", "D-lactate from pyruvate"),
        "SCFA and alcohol conversions: Butyrate, pt 1": ("butyrate", "Buk path (terminal) [ptb-only marker]"),
        "SCFA and alcohol conversions: Butyrate, pt 2": ("butyrate", "Buk path (terminal) [buk marker]"),
        "SCFA and alcohol conversions: Propionate, pt 1": ("propionate", "propionate kinase only"),
        "SCFA and alcohol conversions: Propionate, pt 2": ("propionate", "propionate CoA-transferase only"),
        "SCFA and alcohol conversions: Alcohol production": ("ethanol", "any adh marker"),
    }

    pathway_lookup = pathway_summary.set_index(["scfa", "pathway"])

    rows = []
    for dram_col in DRAM_SCFA_COLUMNS:
        if dram_col not in product.columns:
            continue
        dram_val = product.at[genome, dram_col]
        scfa_key, pw_key = mapping.get(dram_col, (None, None))
        ours_call = None
        ours_evidence = ""
        if scfa_key is not None:
            scfa_row = scfa_summary[scfa_summary["scfa"] == scfa_key]
            if not scfa_row.empty:
                ours_call = bool(scfa_row["scfa_producible"].iloc[0])
                ours_evidence = (
                    f"max_essential_completeness="
                    f"{float(scfa_row['max_essential_completeness'].iloc[0]):.2f}; "
                    f"complete pathways="
                    f"{int(scfa_row['n_pathways_complete'].iloc[0])}/"
                    f"{int(scfa_row['n_pathways_checked'].iloc[0])}; "
                    f"partial pathways="
                    f"{int(scfa_row['n_pathways_partial'].iloc[0])}"
                )
        rows.append({
            "dram_product_column": dram_col,
            "dram_value": dram_val,
            "scfa_mapped": scfa_key,
            "pathway_mapped": pw_key,
            "our_scfa_producible": ours_call,
            "evidence": ours_evidence,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-dir", required=True, help="DRAM contig-annotation sample directory")
    ap.add_argument(
        "--reference",
        default=str(Path(__file__).parent / "scfa_reference.tsv"),
        help="Path to the curated SCFA reference TSV",
    )
    ap.add_argument(
        "--out-dir",
        default=None,
        help="Where to write annotated TSVs (default: <sample-dir>/scfa_extra/)",
    )
    args = ap.parse_args()

    sample_dir = Path(args.sample_dir).resolve()
    annotations_path = sample_dir / "distill_input_fixed" / "annotations.tsv"
    product_tsv = sample_dir / "distillate" / "product.tsv"
    if not annotations_path.is_file():
        raise SystemExit(f"ERROR: missing {annotations_path}")
    if not product_tsv.is_file():
        raise SystemExit(f"ERROR: missing {product_tsv}")
    reference_path = Path(args.reference).resolve()
    if not reference_path.is_file():
        raise SystemExit(f"ERROR: missing reference table {reference_path}")
    out_dir = Path(args.out_dir).resolve() if args.out_dir else (sample_dir / "scfa_extra")
    out_dir.mkdir(parents=True, exist_ok=True)

    annotations = pd.read_csv(annotations_path, sep="\t", index_col=0)
    reference = pd.read_csv(reference_path, sep="\t", dtype=str).fillna("")

    observed = build_observed_id_counter(annotations)

    per_step = annotate_reference(reference, observed)
    pathway_summary = summarise_pathways(per_step)
    scfa_summary = summarise_scfa(pathway_summary)

    per_step.to_csv(out_dir / "scfa_per_step.tsv", sep="\t", index=False)
    pathway_summary.to_csv(out_dir / "scfa_per_pathway.tsv", sep="\t", index=False)
    scfa_summary.to_csv(out_dir / "scfa_per_scfa.tsv", sep="\t", index=False)

    diff = compare_to_dram(scfa_summary, pathway_summary, product_tsv)
    diff.to_csv(out_dir / "scfa_vs_dram_product.tsv", sep="\t", index=False)

    print(f"Sample: {sample_dir}")
    print(f"Wrote per-step    -> {out_dir / 'scfa_per_step.tsv'}")
    print(f"Wrote per-pathway -> {out_dir / 'scfa_per_pathway.tsv'}")
    print(f"Wrote per-scfa    -> {out_dir / 'scfa_per_scfa.tsv'}")
    print(f"Wrote diff vs DRAM-> {out_dir / 'scfa_vs_dram_product.tsv'}")
    print()
    print("=" * 78)
    print("Per-SCFA summary (essential-step completeness):")
    print("=" * 78)
    print(scfa_summary.to_string(index=False))
    print()
    print("=" * 78)
    print("Per-pathway summary:")
    print("=" * 78)
    print(pathway_summary.to_string(index=False))
    print()
    print("=" * 78)
    print("Diff vs DRAM product.tsv SCFA columns:")
    print("=" * 78)
    print(diff.to_string(index=False))


if __name__ == "__main__":
    main()
