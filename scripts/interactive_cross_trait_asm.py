#!/usr/bin/env python3
"""
Interactive & CLI Cross-Trait DNA ASM Recurrence, FDR Filtering, and Functional Concordance Tool.

Loads per-trait `dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv` files, collapses multi-row
CpG/TSS/GTEx entries per unique QUERY_SNP within each trait, applies Benjamini-Hochberg (BH)
or Storey's q-value multiple-testing correction, filters by recurrence across >= k traits,
and outputs a wide cross-trait summary matrix (`results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv`).
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

try:
    from compare_overlap_snps import storey_qvalues as spline_storey_qvalues
except ImportError:
    spline_storey_qvalues = None

DEFAULT_INPUT_DIR = "results/methylation"
DEFAULT_OUTPUT_FILE = "results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv"


def compute_storey_qvalues(p_values: np.ndarray) -> tuple:
    """
    Calculate Storey's q-values:
    Uses spline-based pi0 estimation when available, falling back to the
    Storey (2002) lambda=0.5 estimator: pi0 = min(1, sum(p > 0.5) / (0.5 * m)).
    """
    p = np.asarray(p_values, dtype=float)
    valid_mask = ~np.isnan(p)
    p_valid = p[valid_mask]
    m = len(p_valid)
    if m == 0:
        return np.full_like(p, np.nan), 1.0

    if spline_storey_qvalues is not None and m >= 50:
        return spline_storey_qvalues(p)

    est_pi0 = float(min(1.0, np.sum(p_valid > 0.5) / (0.5 * m)))
    est_pi0 = max(est_pi0, 1.0 / m)

    sort_idx = np.argsort(p_valid)
    p_sorted = p_valid[sort_idx]
    ranks = np.arange(1, m + 1, dtype=float)
    q_sorted = (est_pi0 * m * p_sorted) / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    q_valid = np.empty_like(p_valid)
    q_valid[sort_idx] = q_sorted

    q_full = np.full_like(p, np.nan, dtype=float)
    q_full[valid_mask] = q_valid
    return q_full, est_pi0


def discover_trait_comparison_files(input_dir: str = DEFAULT_INPUT_DIR) -> dict:
    """
    Discover per-trait `dna_asm_*_tss_vs_gtex_full_comparison.tsv` files in either:
      - results/methylation/dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv
      - results/methylation/{trait}/dna_asm_tss_vs_gtex_full_comparison.tsv
    Returns dict {trait_name: filepath}.
    """
    trait_files = {}

    # Pattern 1: Flat layout
    flat_pattern = os.path.join(input_dir, "dna_asm_*_tss_vs_gtex_full_comparison.tsv")
    for path in sorted(glob.glob(flat_pattern)):
        base = os.path.basename(path)
        m = re.match(r"^dna_asm_(.+)_tss_vs_gtex_full_comparison\.tsv$", base)
        if m:
            trait_files[m.group(1)] = path

    # Pattern 2: Nested layout
    nested_pattern = os.path.join(input_dir, "*", "dna_asm_tss_vs_gtex_full_comparison.tsv")
    for path in sorted(glob.glob(nested_pattern)):
        trait = os.path.basename(os.path.dirname(path))
        if trait not in trait_files:
            trait_files[trait] = path

    return trait_files


def join_unique_clean(values, exclude: set = None) -> str:
    exclude_upper = {e.upper() for e in (exclude or {"", "NA", "NAN", "NONE", "<NA>", "NONE_WITHIN_5KB"})}
    items = []
    seen = set()
    for val in values:
        if val is None or (isinstance(val, float) and np.isnan(val)):
            continue
        for token in str(val).split(","):
            t = token.strip()
            if t and t.upper() not in exclude_upper and t.upper() not in seen:
                seen.add(t.upper())
                items.append(t)
    return ", ".join(sorted(items)) if items else "None"


def aggregate_snp_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fast dictionary-based aggregation of CpG, TSS, and GTEx metadata per QUERY_SNP.
    Avoids slow Pandas groupby(Series) object instantiation across tens of thousands of groups.
    """
    if df.empty:
        return pd.DataFrame(
            columns=[
                "QUERY_SNP",
                "CPG_IDS",
                "TAG_SNPS",
                "RELATIONSHIP",
                "TSS_GENES_5KB",
                "MIN_DISTANCE_TO_TSS",
                "GTEX_BRAIN_GENES",
                "CONCORDANT_GENES",
                "GTEX_TISSUES",
            ]
        )

    min_dist = (
        df.groupby("QUERY_SNP", sort=False)["DISTANCE_TO_TSS"]
        .min()
        .reset_index()
        .rename(columns={"DISTANCE_TO_TSS": "MIN_DISTANCE_TO_TSS"})
    )

    cols = [
        ("CPG_IDS", "CPG_ID"),
        ("TAG_SNPS", "TAG_SNP"),
        ("RELATIONSHIP", "RELATIONSHIP"),
        ("TSS_GENES_5KB", "TSS_GENE_SYMBOL"),
        ("GTEX_BRAIN_GENES", "GTEX_GENE_SYMBOL"),
        ("CONCORDANT_GENES", "_CONC_GENE"),
        ("GTEX_TISSUES", "GTEX_TISSUES"),
    ]

    snp_arr = df["QUERY_SNP"].values
    agg_data = {out_col: {} for out_col, _ in cols}

    for out_col, src_col in cols:
        bucket = {}
        vals = df[src_col].values if src_col in df.columns else np.repeat("", len(df))
        for snp, v in zip(snp_arr, vals):
            if snp not in bucket:
                bucket[snp] = []
            bucket[snp].append(v)
        agg_data[out_col] = {snp: join_unique_clean(vlist) for snp, vlist in bucket.items()}

    unique_snps = min_dist["QUERY_SNP"].tolist()
    meta_df = pd.DataFrame({"QUERY_SNP": unique_snps})
    for out_col, _ in cols:
        col_map = agg_data[out_col]
        meta_df[out_col] = [col_map.get(s, "None") for s in unique_snps]

    return pd.merge(meta_df, min_dist, on="QUERY_SNP", how="left")


def deduplicate_trait_table(df: pd.DataFrame, trait_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Collapse multi-row CpG/TSS/GTEx entries in a single trait's comparison file
    to 1 row per unique QUERY_SNP so each variant is tested once for FDR.
    Returns (base_snp_df, cleaned_raw_df).
    """
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()

    df = df.copy()
    df["QUERY_SNP"] = df["QUERY_SNP"].astype(str).str.strip()
    df["GWAS_PVALUE"] = pd.to_numeric(df["GWAS_PVALUE"], errors="coerce")
    df["GWAS_BETA"] = pd.to_numeric(df["GWAS_BETA"], errors="coerce")
    df["DISTANCE_TO_TSS"] = pd.to_numeric(df["DISTANCE_TO_TSS"], errors="coerce")
    df = df.dropna(subset=["QUERY_SNP", "GWAS_PVALUE"])

    # Identify concordant gene symbol where IS_CONCORDANT is True
    is_conc = df["IS_CONCORDANT"].astype(str).str.lower().isin(["true", "1"])
    df["_CONC_GENE"] = np.where(is_conc, df["TSS_GENE_SYMBOL"].astype(str).str.strip(), "")
    conc_snps = set(df.loc[is_conc & (df["_CONC_GENE"] != ""), "QUERY_SNP"].unique())

    # 1 row per unique QUERY_SNP for multiple-testing correction
    base_snp = (
        df.sort_values("DISTANCE_TO_TSS", ascending=True, na_position="last")
        .drop_duplicates(subset=["QUERY_SNP"], keep="first")[
            ["QUERY_SNP", "GWAS_PVALUE", "GWAS_BETA", "DISTANCE_TO_TSS"]
        ]
        .rename(columns={"DISTANCE_TO_TSS": "MIN_DISTANCE_TO_TSS"})
        .copy()
    )
    base_snp["HAS_CONCORDANT"] = base_snp["QUERY_SNP"].isin(conc_snps)
    base_snp["TRAIT"] = trait_name
    return base_snp, df


def run_cross_trait_asm_analysis(
    input_dir: str = DEFAULT_INPUT_DIR,
    method: str = "storey",
    threshold: float = 0.05,
    min_traits: int = 1,
    output_path: str = DEFAULT_OUTPUT_FILE,
    selected_traits: list = None,
    save_output: bool = True,
) -> tuple:
    """
    Core engine for both CLI and Streamlit Dashboard:
    Returns (wide_recurrent_df, all_snps_recurrence_counts_df, long_effects_df, trait_stats_list).
    """
    trait_files = discover_trait_comparison_files(input_dir)
    if selected_traits:
        trait_files = {t: p for t, p in trait_files.items() if t in selected_traits}

    if not trait_files:
        raise FileNotFoundError(
            f"No 'dna_asm_*_tss_vs_gtex_full_comparison.tsv' files found in {input_dir}"
        )

    traits = sorted(trait_files.keys())
    per_trait_dedup = {}
    raw_frames = []
    trait_stats = []

    for trait in traits:
        raw_df = pd.read_csv(trait_files[trait], sep="\t", low_memory=False)
        dedup_df, cleaned_raw = deduplicate_trait_table(raw_df, trait)
        if dedup_df.empty:
            continue

        pvals = dedup_df["GWAS_PVALUE"].values
        if method.lower() == "storey":
            adj_p, est_pi0 = compute_storey_qvalues(pvals)
        else:
            _, adj_p, _, _ = multipletests(pvals, alpha=threshold, method="fdr_bh")
            est_pi0 = 1.0

        dedup_df["ADJ_P"] = adj_p
        dedup_df["IS_SIG"] = dedup_df["ADJ_P"] <= threshold
        per_trait_dedup[trait] = dedup_df
        raw_frames.append(cleaned_raw)

        n_sig = int(dedup_df["IS_SIG"].sum())
        n_conc_sig = int((dedup_df["IS_SIG"] & dedup_df["HAS_CONCORDANT"]).sum())
        trait_stats.append(
            {
                "Trait": trait,
                "Candidate_SNPs": len(dedup_df),
                "Pi0": round(est_pi0, 4),
                "Significant_SNPs": n_sig,
                "Significant_Concordant_SNPs": n_conc_sig,
            }
        )
        print(
            f"[cross_trait_asm] {trait}: {len(dedup_df):,} unique SNPs | "
            f"{n_sig:,} significant at {method.upper()} <= {threshold} "
            f"(pi0={est_pi0:.4f}, {n_conc_sig:,} with concordant TSS=GTEx gene)."
        )

    if not per_trait_dedup:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), []

    all_long = pd.concat(list(per_trait_dedup.values()), ignore_index=True)

    # Filter to trait-significant rows (ADJ_P <= threshold) for recurrence counting
    sig_long = all_long[all_long["IS_SIG"]].copy()
    if sig_long.empty:
        recurrence_df = pd.DataFrame(columns=["QUERY_SNP", "N_TRAITS", "TRAITS"])
    else:
        recurrence_df = (
            sig_long.groupby("QUERY_SNP", as_index=False)
            .agg(
                N_TRAITS=("TRAIT", "nunique"),
                TRAITS=("TRAIT", lambda s: ", ".join(sorted(set(s)))),
            )
        )

    # Aggregate metadata for SNPs that pass significance in at least 1 trait (or all if min_traits <= 0)
    all_raw_concat = pd.concat(raw_frames, ignore_index=True)
    if min_traits >= 1:
        sig_snps_set = set(recurrence_df["QUERY_SNP"].unique())
        raw_subset = all_raw_concat[all_raw_concat["QUERY_SNP"].isin(sig_snps_set)]
    else:
        raw_subset = all_raw_concat

    meta_agg = aggregate_snp_metadata(raw_subset)
    if meta_agg.empty:
        full_matrix = pd.merge(meta_agg, recurrence_df, on="QUERY_SNP", how="left")
    else:
        full_matrix = pd.merge(meta_agg, recurrence_df, on="QUERY_SNP", how="left")
    full_matrix["N_TRAITS"] = full_matrix["N_TRAITS"].fillna(0).astype(int)
    full_matrix["TRAITS"] = full_matrix["TRAITS"].fillna("None")

    # Attach dynamic per-trait columns: BETA_{T} and ADJ_P_{T} (significant hits)
    dynamic_cols = []
    adj_p_cols = []
    for trait in traits:
        t_df = per_trait_dedup[trait]
        t_sig = t_df[t_df["IS_SIG"]][["QUERY_SNP", "GWAS_BETA", "ADJ_P"]].rename(
            columns={
                "GWAS_BETA": f"BETA_{trait}",
                "ADJ_P": f"ADJ_P_{trait}",
            }
        )
        full_matrix = pd.merge(full_matrix, t_sig, on="QUERY_SNP", how="left")
        dynamic_cols.extend([f"BETA_{trait}", f"ADJ_P_{trait}"])
        adj_p_cols.append(f"ADJ_P_{trait}")

    full_matrix = full_matrix.rename(columns={"QUERY_SNP": "SNP"})
    full_matrix["MIN_DISTANCE_TO_TSS"] = pd.to_numeric(
        full_matrix["MIN_DISTANCE_TO_TSS"], errors="coerce"
    ).astype("Int64")

    # Distribution of N_TRAITS across all SNPs with N_TRAITS >= 1
    sig_any = full_matrix[full_matrix["N_TRAITS"] >= 1].copy()

    # Filter for N_TRAITS >= min_traits
    filtered = full_matrix[full_matrix["N_TRAITS"] >= min_traits].copy()
    if not filtered.empty and adj_p_cols:
        filtered["_MEDIAN_ADJ_P"] = filtered[adj_p_cols].median(axis=1, skipna=True)
        filtered["_HAS_CONC"] = filtered["CONCORDANT_GENES"] != "None"
        filtered = (
            filtered.sort_values(
                by=["N_TRAITS", "_HAS_CONC", "_MEDIAN_ADJ_P"],
                ascending=[False, False, True],
                na_position="last",
            )
            .drop(columns=["_MEDIAN_ADJ_P", "_HAS_CONC"])
            .reset_index(drop=True)
        )

    lead_cols = [
        "SNP",
        "N_TRAITS",
        "TRAITS",
        "CONCORDANT_GENES",
        "CPG_IDS",
        "TAG_SNPS",
        "RELATIONSHIP",
        "TSS_GENES_5KB",
        "MIN_DISTANCE_TO_TSS",
        "GTEX_BRAIN_GENES",
        "GTEX_TISSUES",
    ]
    ordered_cols = lead_cols + dynamic_cols
    filtered = filtered[ordered_cols]

    if save_output and output_path:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        filtered.to_csv(output_path, sep="\t", index=False, na_rep="NA")
        print(
            f"[cross_trait_asm] Saved {len(filtered):,} recurrent SNPs (N_TRAITS >= {min_traits}) to {output_path}"
        )

    return filtered, sig_any, all_long, trait_stats


def interactive_cli(input_dir: str, default_threshold: float, default_output: str) -> None:
    trait_files = discover_trait_comparison_files(input_dir)
    if not trait_files:
        print(f"Error: No 'dna_asm_*_tss_vs_gtex_full_comparison.tsv' files found in '{input_dir}'.")
        sys.exit(1)

    traits = sorted(trait_files.keys())
    n_total = len(traits)

    print("\n" + "=" * 64)
    print("   DNA ASM Cross-Trait FDR & Functional Convergence Explorer")
    print("=" * 64)
    print(f"Discovered {n_total} GWAS traits: {', '.join(traits)}")

    print("\nMultiple-Testing Correction Method:")
    print("  [1] Benjamini-Hochberg FDR (fdr_bh)")
    print("  [2] Storey's q-value (pi0-adjusted FDR) [Default]")
    m_in = input("Select method [1 or 2, default: 2]: ").strip()
    method = "bh" if m_in == "1" else "storey"

    t_in = input(f"\nEnter significance threshold [default: {default_threshold}]: ").strip()
    try:
        threshold = float(t_in) if t_in else default_threshold
    except ValueError:
        threshold = default_threshold

    print(f"\nTrait Recurrence Threshold k (1 to {n_total}, or 'all'):")
    default_k = min(2, n_total)
    k_in = input(f"Enter minimum traits k [default: {default_k}]: ").strip().lower()
    if k_in == "all":
        min_k = n_total
    elif k_in:
        try:
            min_k = max(1, min(n_total, int(k_in)))
        except ValueError:
            min_k = default_k
    else:
        min_k = default_k

    out_in = input(f"\nOutput TSV path [default: {default_output}]: ").strip()
    output_path = out_in if out_in else default_output

    print()
    run_cross_trait_asm_analysis(
        input_dir=input_dir,
        method=method,
        threshold=threshold,
        min_traits=min_k,
        output_path=output_path,
        save_output=True,
    )


def main() -> None:
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        output_tsv = smk.output.tsv if hasattr(smk.output, "tsv") else smk.output[0]
        method = getattr(smk.params, "method", "storey")
        min_k = int(getattr(smk.params, "k", 1))
        threshold = float(getattr(smk.params, "threshold", 0.05))
        input_dir = getattr(smk.params, "input_dir", DEFAULT_INPUT_DIR)
        run_cross_trait_asm_analysis(
            input_dir=input_dir,
            method=method,
            threshold=threshold,
            min_traits=min_k,
            output_path=output_tsv,
            save_output=True,
        )
        return

    parser = argparse.ArgumentParser(
        description="Interactive & CLI cross-trait DNA ASM FDR filtering and concordance summary."
    )
    parser.add_argument(
        "--input-dir",
        default=DEFAULT_INPUT_DIR,
        help="Directory containing per-trait DNA ASM comparison files (default: results/methylation)",
    )
    parser.add_argument(
        "--method",
        choices=["bh", "storey"],
        default=None,
        help="Multiple-testing correction method: 'bh' or 'storey'",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="FDR / q-value significance threshold (default: 0.05)",
    )
    parser.add_argument(
        "--min-traits",
        default=None,
        help="Minimum trait recurrence k (integer 1..N or 'all')",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=DEFAULT_OUTPUT_FILE,
        help="Output summary TSV path (default: results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv)",
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="Force interactive prompt mode",
    )
    args = parser.parse_args()

    if args.interactive or (args.method is None and args.min_traits is None):
        interactive_cli(args.input_dir, args.threshold, args.output)
    else:
        trait_files = discover_trait_comparison_files(args.input_dir)
        n_total = max(1, len(trait_files))
        method = args.method or "storey"
        if args.min_traits is None:
            min_k = 1
        elif str(args.min_traits).lower() == "all":
            min_k = n_total
        else:
            min_k = int(args.min_traits)

        run_cross_trait_asm_analysis(
            input_dir=args.input_dir,
            method=method,
            threshold=args.threshold,
            min_traits=min_k,
            output_path=args.output,
            save_output=True,
        )


if __name__ == "__main__":
    main()

