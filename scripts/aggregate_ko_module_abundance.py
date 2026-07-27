#!/usr/bin/env python3
"""
Aggregate gene-level read counts (from featureCounts) into per-sample
ID-level abundance (KO/EC/Pfam/CAZy) and module-level abundance, with
TPM normalisation. Covers everything DRAM puts into product.tsv plus
the curated SCFA reference table.

Pipeline:
  featureCounts gene_counts.tsv  +  DRAM annotations.tsv
        |
        v
  Step 1: gene -> TPM (length + library-size normalised)
  Step 2: gene -> ID (KO/EC/Pfam/CAZy/peptidase) via fractional allocation
            (a gene's TPM is split evenly across the IDs it maps to)
  Step 3a: ID TPM -> SCFA module abundance (scripts/scfa_reference.tsv)
  Step 3b: ID TPM -> KEGG module step abundance (module_step_form.tsv)
  Step 3c: ID TPM -> ETC complex abundance (etc_module_database.tsv)
  Step 3d: ID TPM -> functional True/False abundance
            (function_heatmap_form.tsv: CAZy, Nitrogen, Sulfur,
            Other Reductases, Photosynthesis, Methanogenesis,
            SCFA/Alcohol pyruvate routes, CAMPER)

Outputs (in --out-dir):
  gene_tpm.tsv
  id_tpm.tsv                          (KO/EC/Pfam/CAZy/peptidase combined)
  scfa_per_step_abundance.tsv
  scfa_per_pathway_abundance.tsv
  scfa_per_scfa_abundance.tsv
  kegg_module_abundance.tsv           (only if --module-step-form given)
  etc_complex_abundance.tsv           (only if --etc-module-database given)
  functional_abundance.tsv            (only if --function-heatmap-form given)
"""
import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def parse_featurecounts(counts_path: Path) -> pd.DataFrame:
    """
    Read the default featureCounts TSV (leading '# Program: ...' line is
    handled by `comment='#'`). Returns DataFrame with columns
    gene_id, length, count.
    """
    df = pd.read_csv(counts_path, sep="\t", comment="#")
    meta_cols = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    sample_cols = [c for c in df.columns if c not in meta_cols]
    if not sample_cols:
        raise SystemExit(f"ERROR: no sample count columns in {counts_path}")
    if len(sample_cols) > 1:
        print(
            f"WARN: {counts_path} has {len(sample_cols)} sample columns; "
            f"using only first: {sample_cols[0]}"
        )
    df = df.rename(
        columns={"Geneid": "gene_id", "Length": "length", sample_cols[0]: "count"}
    )
    return df[["gene_id", "length", "count"]].copy()


def compute_tpm(gene_df: pd.DataFrame) -> pd.DataFrame:
    """
    TPM = (count/length_kb) / sum(count/length_kb) * 1e6.
    Result has columns gene_id, length, count, RPK, TPM.
    """
    gene_df = gene_df.copy()
    length_kb = gene_df["length"].astype(float) / 1000.0
    gene_df["RPK"] = gene_df["count"].astype(float) / length_kb
    total = gene_df["RPK"].sum()
    gene_df["TPM"] = (gene_df["RPK"] / total * 1e6) if total > 0 else 0.0
    return gene_df


def extract_ids_by_gene(annotations: pd.DataFrame):
    """Use DRAM's own row-wise ID extractor so we match what distill uses."""
    from mag_annotator.summarize_genomes import get_ids_from_annotations_by_row
    by_row = get_ids_from_annotations_by_row(annotations)
    return {g: set(ids) for g, ids in by_row.items()}


def gene_tpm_to_id_tpm(
    gene_tpm: pd.DataFrame,
    gene_to_ids: dict,
) -> pd.DataFrame:
    """
    Fractional allocation: a gene's TPM is split evenly across its IDs.
    Returns (id, tpm, n_genes) sorted by tpm desc.
    """
    id_tpm = defaultdict(float)
    id_genecount = defaultdict(int)
    for _, row in gene_tpm.iterrows():
        gid = row["gene_id"]
        tpm = float(row["TPM"])
        ids = gene_to_ids.get(gid, set())
        if not ids:
            continue
        share = tpm / len(ids)
        for i in ids:
            id_tpm[i] += share
            id_genecount[i] += 1
    out = pd.DataFrame({"id": sorted(id_tpm.keys())})
    out["tpm"] = out["id"].map(id_tpm)
    out["n_genes"] = out["id"].map(id_genecount)
    # Stable: tie-break on `id` so output is bytewise reproducible.
    return out.sort_values(
        ["tpm", "id"], ascending=[False, True]
    ).reset_index(drop=True)


def split_ids_string(s) -> set:
    if not s or (isinstance(s, float) and pd.isna(s)):
        return set()
    parts = re.split(r"[,;\s]+", str(s).strip())
    return {p for p in parts if p}


def summarise_scfa_abundance(scfa_ref: pd.DataFrame, id_tpm_map: dict):
    """
    Per-step / per-pathway / per-SCFA TPM aggregation following
    scripts/scfa_reference.tsv.

    Step rule:
        step_tpm = sum of TPM over the IDs listed in 'ids' (these are
                   alternative IDs for that step; we sum because each
                   ID's TPM was already fractionally allocated upstream).

    Pathway rule (only essential steps count for completeness):
        completeness         = #essential steps with tpm>0 / #essential steps
        pathway_sum_tpm      = sum over essential step_tpm
        pathway_mean_tpm     = mean over essential step_tpm
        pathway_min_tpm      = min over essential step_tpm (bottleneck)
        pathway_weighted_tpm = pathway_sum_tpm * completeness
        pathway_complete     = (completeness == 1.0 and n_essential > 0)

    SCFA rule:
        complete_pathway_max_mean_tpm = max pathway_mean_tpm among
                                        complete pathways
        scfa_total_sum_tpm            = sum over pathways of pathway_sum_tpm
        scfa_max_completeness         = max essential_completeness over pathways
    """
    per_step = []
    for _, row in scfa_ref.iterrows():
        # `split_ids_string` returns a set, which has hash-randomised
        # iteration order. Sort before summing so floating-point round-off
        # is identical across runs (was the source of the 1.42e-14 epsilon
        # drift we previously saw in scfa_per_step_abundance.tsv).
        ids = sorted(split_ids_string(row.get("ids", "")))
        step_tpm = float(sum(id_tpm_map.get(i, 0.0) for i in ids))
        per_step.append({
            **row.to_dict(),
            "step_tpm": step_tpm,
            "step_present": step_tpm > 0,
        })
    per_step_df = pd.DataFrame(per_step)

    pathway_rows = []
    for (scfa, pathway), g in per_step_df.groupby(["scfa", "pathway"], sort=False):
        ess = g[g["essential"].astype(str).str.lower() == "yes"]
        n_ess = len(ess)
        present = int((ess["step_tpm"] > 0).sum())
        completeness = (present / n_ess) if n_ess else float("nan")
        complete = bool(n_ess > 0 and present == n_ess)
        ess_tpm = ess["step_tpm"]
        pathway_rows.append({
            "scfa": scfa,
            "pathway": pathway,
            "n_essential_steps": n_ess,
            "n_essential_steps_present": present,
            "essential_completeness": completeness,
            "pathway_complete": complete,
            "pathway_sum_tpm": float(ess_tpm.sum()) if n_ess else 0.0,
            "pathway_mean_tpm": float(ess_tpm.mean()) if n_ess else float("nan"),
            "pathway_min_essential_tpm": float(ess_tpm.min()) if n_ess else float("nan"),
            "pathway_weighted_tpm": float(ess_tpm.sum() * completeness)
            if not (isinstance(completeness, float) and np.isnan(completeness))
            else float("nan"),
        })
    pathway_df = pd.DataFrame(pathway_rows)

    scfa_rows = []
    for scfa, g in pathway_df.groupby("scfa", sort=False):
        complete = g[g["pathway_complete"]]
        scfa_rows.append({
            "scfa": scfa,
            "n_pathways_checked": len(g),
            "n_pathways_complete": len(complete),
            "scfa_max_pathway_completeness": float(g["essential_completeness"].max()),
            "scfa_complete_pathway_max_mean_tpm": (
                float(complete["pathway_mean_tpm"].max()) if len(complete) else 0.0
            ),
            "scfa_total_sum_tpm": float(g["pathway_sum_tpm"].sum()),
            "scfa_total_weighted_tpm": float(g["pathway_weighted_tpm"].fillna(0).sum()),
        })
    scfa_df = pd.DataFrame(scfa_rows)

    return per_step_df, pathway_df, scfa_df


def _heatmap_modules() -> set:
    try:
        from mag_annotator.summarize_genomes import HEATMAP_MODULES
        return set(HEATMAP_MODULES)
    except Exception:
        return set()


def summarise_kegg_module_abundance(
    module_step_form: pd.DataFrame, id_tpm_map: dict
) -> pd.DataFrame:
    """
    Re-aggregate KEGG modules at the step level using TPM.

    A step is defined by the first integer in `path` (e.g. "0,0" -> step 0).
    Within a step the multiple `ko` rows are OR alternatives, so we sum
    their TPM as the step abundance.

    Returns one row per module with:
      n_steps, n_steps_present, completeness,
      module_sum_tpm   = sum over steps,
      module_mean_step_tpm,
      module_min_step_tpm,
      module_weighted_tpm = module_sum_tpm * completeness.
    """
    msf = module_step_form.copy()
    msf["step"] = msf["path"].astype(str).str.split(",").str[0].astype(int)
    heatmap_modules = _heatmap_modules()

    rows = []
    for module, g in msf.groupby("module", sort=False):
        module_name = g["module_name"].iloc[0]
        step_tpms = []
        for step, sg in g.groupby("step"):
            kos = sg["ko"].astype(str).str.strip().tolist()
            step_tpm = float(sum(id_tpm_map.get(k, 0.0) for k in kos))
            step_tpms.append(step_tpm)
        n_steps = len(step_tpms)
        if n_steps == 0:
            continue
        n_present = sum(1 for t in step_tpms if t > 0)
        completeness = n_present / n_steps
        rows.append({
            "module": module,
            "module_name": module_name,
            "in_product_heatmap": module in heatmap_modules,
            "n_steps": n_steps,
            "n_steps_present": n_present,
            "completeness": completeness,
            "module_sum_tpm": float(sum(step_tpms)),
            "module_mean_step_tpm": float(np.mean(step_tpms)),
            "module_min_step_tpm": float(min(step_tpms)),
            "module_weighted_tpm": float(sum(step_tpms) * completeness),
        })
    return pd.DataFrame(rows).sort_values(
        ["module_weighted_tpm", "module"], ascending=[False, True]
    ).reset_index(drop=True)


def summarise_etc_complex_abundance(
    etc_module_df: pd.DataFrame, id_tpm_map: dict
) -> pd.DataFrame:
    """
    For each ETC complex in etc_module_database, parse `definition` into a
    DAG (same as DRAM's make_etc_coverage_df), enumerate all simple paths
    from start to end, score each path by summing TPM over its KOs, and
    take the best path. Optional subunits (-Knnnnn) are stripped first,
    consistent with DRAM.

    Output columns:
      module_id, module_name, complex, complex_module_name,
      path_length, n_kos_present, percent_coverage,
      best_path_sum_tpm   = sum TPM on the path with max coverage
      best_path_mean_tpm  = mean TPM across KOs on that path
      best_path_min_tpm   = bottleneck KO TPM on that path
      module_weighted_tpm = best_path_sum_tpm * percent_coverage
    """
    import networkx as nx
    from mag_annotator.summarize_genomes import make_module_network

    rows = []
    for _, module_row in etc_module_df.iterrows():
        definition = re.sub(r"-K\d\d\d\d\d", "", module_row["definition"])
        module_net, _ = make_module_network(definition)
        no_out = [n for n in module_net.nodes() if module_net.out_degree(n) == 0]
        for n in no_out:
            module_net.add_edge(n, "end")

        observed = {ko for ko in id_tpm_map.keys() if ko in module_net.nodes}
        best = None
        for path in nx.all_simple_paths(module_net, source="start", target="end"):
            kos = [n for n in path if n not in ("start", "end")]
            if not kos:
                continue
            present = [k for k in kos if k in observed]
            coverage = len(present) / len(kos)
            path_tpms = [id_tpm_map.get(k, 0.0) for k in kos]
            path_sum = float(sum(path_tpms))
            cand = {
                "kos": kos,
                "present": present,
                "coverage": coverage,
                "sum_tpm": path_sum,
                "mean_tpm": float(np.mean(path_tpms)),
                "min_tpm": float(min(path_tpms)) if path_tpms else 0.0,
            }
            if best is None or (
                cand["coverage"] > best["coverage"]
                or (cand["coverage"] == best["coverage"] and cand["sum_tpm"] > best["sum_tpm"])
            ):
                best = cand
        if best is None:
            continue
        complex_module_name = "Complex %s: %s" % (
            str(module_row["complex"]).replace("Complex ", ""),
            module_row["module_name"],
        )
        rows.append({
            "module_id": module_row["module_id"],
            "module_name": module_row["module_name"],
            "complex": module_row["complex"],
            "complex_module_name": complex_module_name,
            "path_length": len(best["kos"]),
            "n_kos_present": len(best["present"]),
            "percent_coverage": best["coverage"],
            "best_path_sum_tpm": best["sum_tpm"],
            "best_path_mean_tpm": best["mean_tpm"],
            "best_path_min_tpm": best["min_tpm"],
            "module_weighted_tpm": best["sum_tpm"] * best["coverage"],
            "best_path_kos": ",".join(best["kos"]),
            "best_path_present_kos": ",".join(sorted(best["present"])),
        })
    return pd.DataFrame(rows).sort_values(
        ["module_weighted_tpm", "complex_module_name"], ascending=[False, True]
    ).reset_index(drop=True)


def summarise_functional_abundance(
    function_heatmap_form: pd.DataFrame, id_tpm_map: dict
) -> pd.DataFrame:
    """
    For each (category, function_name) in function_heatmap_form, score
    each rule row by summing TPM over its `function_ids` (OR alternatives),
    then aggregate rules (AND relationship in DRAM's True/False logic).

    DRAM's `.hmm` normalisation: function_ids may be bare CAZy families
    (e.g., 'GH13'); after our upstream `normalize_dram_fasta_ids.py` step
    the annotation IDs are also bare (no `.hmm` suffix), so an exact
    match works. To be safe, we also try the `.hmm`-suffixed variant.

    Output columns:
      category, subcategory, function_name, product_column,
      n_rules, n_rules_passed, completeness,
      function_complete   = (completeness == 1.0)
      function_sum_tpm    = sum over rules of rule_sum_tpm
      function_mean_rule_tpm
      function_min_rule_tpm  (bottleneck across rules)
      function_weighted_tpm  = function_sum_tpm * completeness
      per_rule_detail
    """
    fhf = function_heatmap_form.copy()
    # strip whitespace consistently with DRAM's make_functional_df
    fhf = fhf.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    def lookup_tpm(token: str) -> float:
        if token in id_tpm_map:
            return id_tpm_map[token]
        if token + ".hmm" in id_tpm_map:
            return id_tpm_map[token + ".hmm"]
        return 0.0

    rows = []
    for (cat, subcat, fname), g in fhf.groupby(
        ["category", "subcategory", "function_name"], sort=False
    ):
        rule_tpms = []
        rule_details = []
        for _, rule in g.iterrows():
            ids = sorted(split_ids_string(rule.get("function_ids", "")))
            per_id = {i: lookup_tpm(i) for i in ids}
            rule_sum = float(sum(per_id.values()))
            hits = [i for i, t in per_id.items() if t > 0]
            rule_tpms.append(rule_sum)
            rule_details.append(
                f"[{rule.get('long_function_name', fname)}] "
                f"sum_tpm={rule_sum:.2f}; "
                f"hits={len(hits)}/{len(ids)}; "
                f"hit_ids={','.join(hits) if hits else '-'}"
            )

        n_rules = len(rule_tpms)
        n_passed = sum(1 for t in rule_tpms if t > 0)
        completeness = (n_passed / n_rules) if n_rules else float("nan")
        rows.append({
            "category": cat,
            "subcategory": subcat,
            "function_name": fname,
            "product_column": f"{cat}: {fname}",
            "n_rules": n_rules,
            "n_rules_passed": n_passed,
            "completeness": completeness,
            "function_complete": bool(n_rules > 0 and n_passed == n_rules),
            "function_sum_tpm": float(sum(rule_tpms)),
            "function_mean_rule_tpm": float(np.mean(rule_tpms)) if n_rules else float("nan"),
            "function_min_rule_tpm": float(min(rule_tpms)) if n_rules else float("nan"),
            "function_weighted_tpm": float(sum(rule_tpms) * completeness)
            if not (isinstance(completeness, float) and np.isnan(completeness))
            else float("nan"),
            "per_rule_detail": " || ".join(rule_details),
        })
    return pd.DataFrame(rows).sort_values(
        ["category", "function_weighted_tpm", "function_name"],
        ascending=[True, False, True],
    ).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--counts", required=True, help="featureCounts output TSV")
    ap.add_argument("--annotations", required=True, help="DRAM annotations.tsv")
    ap.add_argument(
        "--scfa-reference",
        default=str(Path(__file__).parent / "scfa_reference.tsv"),
        help="Curated SCFA reference TSV (default: scripts/scfa_reference.tsv)",
    )
    ap.add_argument(
        "--module-step-form",
        default=None,
        help="DRAM module_step_form.tsv path (KEGG modules in product.tsv).",
    )
    ap.add_argument(
        "--etc-module-database",
        default=None,
        help="DRAM etc_module_database.tsv path (ETC complexes in product.tsv).",
    )
    ap.add_argument(
        "--function-heatmap-form",
        default=None,
        help="DRAM function_heatmap_form.tsv path (True/False columns in product.tsv).",
    )
    ap.add_argument("--out-dir", required=True)
    ap.add_argument(
        "--sample-label",
        default=None,
        help="Optional label embedded as a 'sample' column in each output",
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== aggregate_ko_module_abundance ===")
    print(f"counts:              {args.counts}")
    print(f"annotations:         {args.annotations}")
    print(f"scfa-reference:      {args.scfa_reference}")
    print(f"module-step-form:    {args.module_step_form}")
    print(f"etc-module-database: {args.etc_module_database}")
    print(f"function-heatmap:    {args.function_heatmap_form}")
    print(f"out-dir:             {out_dir}")
    if args.sample_label:
        print(f"sample-label:     {args.sample_label}")

    # Step 1: gene counts -> TPM
    gene_df = parse_featurecounts(Path(args.counts).resolve())
    gene_tpm = compute_tpm(gene_df)
    if args.sample_label:
        gene_tpm.insert(0, "sample", args.sample_label)
    gene_tpm.to_csv(out_dir / "gene_tpm.tsv", sep="\t", index=False)
    print(
        f"\n[Step 1] gene -> TPM: {len(gene_tpm):,} genes; "
        f"sum TPM = {gene_tpm['TPM'].sum():,.1f} (should be ~1e6)"
    )

    # Step 2: gene -> ID TPM
    annotations = pd.read_csv(args.annotations, sep="\t", index_col=0)
    gene_to_ids = extract_ids_by_gene(annotations)
    n_with_ids = sum(1 for v in gene_to_ids.values() if v)
    id_tpm_df = gene_tpm_to_id_tpm(gene_tpm, gene_to_ids)
    if args.sample_label:
        id_tpm_df.insert(0, "sample", args.sample_label)
    id_tpm_df.to_csv(out_dir / "id_tpm.tsv", sep="\t", index=False)
    print(
        f"[Step 2] gene -> ID TPM: {len(id_tpm_df):,} unique IDs; "
        f"{n_with_ids:,}/{len(gene_tpm):,} genes had >=1 ID"
    )

    id_tpm_map = dict(zip(id_tpm_df["id"], id_tpm_df["tpm"]))

    # Step 3a: SCFA module abundance
    scfa_ref = pd.read_csv(args.scfa_reference, sep="\t", dtype=str).fillna("")
    per_step, pathway_df, scfa_df = summarise_scfa_abundance(scfa_ref, id_tpm_map)
    for df in (per_step, pathway_df, scfa_df):
        if args.sample_label:
            df.insert(0, "sample", args.sample_label)
    per_step.to_csv(out_dir / "scfa_per_step_abundance.tsv", sep="\t", index=False)
    pathway_df.to_csv(out_dir / "scfa_per_pathway_abundance.tsv", sep="\t", index=False)
    scfa_df.to_csv(out_dir / "scfa_per_scfa_abundance.tsv", sep="\t", index=False)
    print(
        f"[Step 3a] SCFA: {len(scfa_df)} SCFAs, {len(pathway_df)} pathways, "
        f"{len(per_step)} steps"
    )

    # Step 3b: KEGG module abundance (optional)
    if args.module_step_form:
        msf = pd.read_csv(args.module_step_form, sep="\t")
        kegg_df = summarise_kegg_module_abundance(msf, id_tpm_map)
        if args.sample_label:
            kegg_df.insert(0, "sample", args.sample_label)
        kegg_df.to_csv(out_dir / "kegg_module_abundance.tsv", sep="\t", index=False)
        print(f"[Step 3b] KEGG modules: {len(kegg_df)} modules")

    # Step 3c: ETC complex abundance (optional)
    if args.etc_module_database:
        etc_df = pd.read_csv(args.etc_module_database, sep="\t")
        etc_out = summarise_etc_complex_abundance(etc_df, id_tpm_map)
        if args.sample_label:
            etc_out.insert(0, "sample", args.sample_label)
        etc_out.to_csv(out_dir / "etc_complex_abundance.tsv", sep="\t", index=False)
        print(f"[Step 3c] ETC complexes: {len(etc_out)} entries")

    # Step 3d: Functional True/False abundance (optional)
    if args.function_heatmap_form:
        fhf = pd.read_csv(args.function_heatmap_form, sep="\t")
        func_out = summarise_functional_abundance(fhf, id_tpm_map)
        if args.sample_label:
            func_out.insert(0, "sample", args.sample_label)
        func_out.to_csv(out_dir / "functional_abundance.tsv", sep="\t", index=False)
        cats = func_out["category"].value_counts().to_dict()
        print(
            f"[Step 3d] Functional True/False: {len(func_out)} functions "
            f"across {len(cats)} categories ({cats})"
        )

    print(f"\nWrote outputs under {out_dir}")


if __name__ == "__main__":
    main()
