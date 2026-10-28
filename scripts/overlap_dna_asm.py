#!/usr/bin/env python3
"""
Step A1: Tag-SNP GWAS Overlap with BH FDR (5%), then LD-Proxy Expansion of Significant Tags.

Reads data/methylation/DNA_ASM_mapping_table.tsv (CPG_ID, TAG_SNP, QUERY_SNP, RELATIONSHIP[, LD_R2])
and runs a two-stage overlap against one trait's GWAS summary statistics:

  Stage 1: one test per unique TAG_SNP (its own GWAS P-value). With --tag-fallback best_ld_proxy,
           a tag absent from the GWAS is represented by its highest-LD_R2 proxy present in the
           GWAS (chosen on LD, never on P). BH FDR is applied across the unique tested SNPs.
  Stage 2: every mapping row (tag + all LD proxies, with CpG and LD_R2) of each tested TAG_SNP
           inherits that tag's test, so proxies absent from the GWAS are retained.

GWAS_P, GWAS_BETA, GWAS_FDR_BH and BH_SIGNIFICANT describe the block's tested SNP (TESTED_SNP);
QUERY_IN_GWAS, QUERY_GWAS_P and QUERY_GWAS_BETA describe the row's own QUERY_SNP.

Outputs:
    --output-all : every mapping row of each tested TAG_SNP
    --output-sig : the subset with BH_SIGNIFICANT == True (input to Steps A2-A4)

Usage:
    python scripts/overlap_dna_asm.py \
        --mapping-table data/methylation/DNA_ASM_mapping_table.tsv \
        --gwas data/raw_data/SUD2023_gwas.txt --trait SUD2023 --alpha 0.05 \
        --tag-fallback none \
        --output-all results/dna/SUD2023_overlap_all.tsv \
        --output-sig results/dna/SUD2023_overlap_sig.tsv
"""

import argparse
import os

import pandas as pd

from common_utils import bh_fdr, gwas_lookup

MAPPING_COLUMNS = ["CPG_ID", "TAG_SNP", "QUERY_SNP", "RELATIONSHIP"]
TAG_FALLBACKS = ("none", "best_ld_proxy")
OUTPUT_COLUMNS = [
    "TRAIT",
    "CPG_ID",
    "TAG_SNP",
    "QUERY_SNP",
    "RELATIONSHIP",
    "GWAS_P",
    "GWAS_BETA",
    "GWAS_FDR_BH",
    "BH_SIGNIFICANT",
    "TESTED_SNP",
    "LD_R2",
    "QUERY_IN_GWAS",
    "QUERY_GWAS_P",
    "QUERY_GWAS_BETA",
]


def load_mapping_table(mapping_table_path: str) -> pd.DataFrame:
    """Read the CpG/tag/proxy mapping table; LD_R2 is optional (NA when absent, 1.0 for tags)."""
    mapping_df = pd.read_csv(mapping_table_path, sep="\t", dtype=str)
    missing = [c for c in MAPPING_COLUMNS if c not in mapping_df.columns]
    if missing:
        raise ValueError(f"Mapping table {mapping_table_path} is missing required columns: {missing}")

    for col in MAPPING_COLUMNS:
        mapping_df[col] = mapping_df[col].astype(str).str.strip()
    mapping_df["LD_R2"] = pd.to_numeric(mapping_df.get("LD_R2"), errors="coerce")
    is_tag = mapping_df["QUERY_SNP"] == mapping_df["TAG_SNP"]
    mapping_df.loc[is_tag, "LD_R2"] = 1.0
    return mapping_df[MAPPING_COLUMNS + ["LD_R2"]].drop_duplicates(subset=MAPPING_COLUMNS)


def select_tested_snps(
    mapping_df: pd.DataFrame, gwas_df: pd.DataFrame, tag_fallback: str
) -> pd.DataFrame:
    """Stage 1 test assignment: one (TAG_SNP, TESTED_SNP) pair per tag with a valid GWAS P-value."""
    valid_snps = set(gwas_df.loc[gwas_df["GWAS_P"].between(0.0, 1.0), "SNP"])
    candidates = mapping_df.loc[
        mapping_df["QUERY_SNP"].isin(valid_snps), ["TAG_SNP", "QUERY_SNP", "LD_R2"]
    ].drop_duplicates(subset=["TAG_SNP", "QUERY_SNP"])
    candidates["_IS_TAG"] = candidates["QUERY_SNP"] == candidates["TAG_SNP"]
    if tag_fallback == "none":
        candidates = candidates[candidates["_IS_TAG"]]
    elif candidates.loc[~candidates["_IS_TAG"], "LD_R2"].isna().any():
        raise ValueError(
            "--tag-fallback best_ld_proxy needs LD_R2 for every proxy; regenerate the mapping "
            "table with scripts/ldsnps_filtering.py."
        )

    # Prefer the tag itself, then the highest-LD proxy (ties broken by rsID for determinism)
    tested = (
        candidates.sort_values(
            by=["TAG_SNP", "_IS_TAG", "LD_R2", "QUERY_SNP"], ascending=[True, False, False, True]
        )
        .drop_duplicates(subset=["TAG_SNP"], keep="first")
        .rename(columns={"QUERY_SNP": "TESTED_SNP"})
    )
    return tested[["TAG_SNP", "TESTED_SNP"]].reset_index(drop=True)


def run_overlap_dna_asm(
    mapping_table_path: str,
    gwas_path: str,
    trait: str,
    output_all: str,
    output_sig: str,
    alpha: float = 0.05,
    tag_fallback: str = "none",
) -> pd.DataFrame:
    if tag_fallback not in TAG_FALLBACKS:
        raise ValueError(f"tag_fallback must be one of {TAG_FALLBACKS}, got '{tag_fallback}'.")

    mapping_df = load_mapping_table(mapping_table_path)
    n_tags = mapping_df["TAG_SNP"].nunique()
    print(
        f"[overlap_dna_asm] {trait}: {len(mapping_df):,} mapping rows "
        f"({n_tags:,} unique TAG_SNPs, {mapping_df['QUERY_SNP'].nunique():,} unique QUERY_SNPs)."
    )

    # One GWAS pass over tags + proxies: tag P-values feed Stage 1, proxy values are annotation only
    gwas_df = gwas_lookup(gwas_path, set(mapping_df["QUERY_SNP"]))

    # Stage 1: BH FDR across the unique tested SNPs (one per tag)
    tested = select_tested_snps(mapping_df, gwas_df, tag_fallback)
    family = gwas_df[gwas_df["SNP"].isin(tested["TESTED_SNP"])].copy()
    family["GWAS_FDR_BH"] = bh_fdr(family["GWAS_P"].to_numpy())
    family["BH_SIGNIFICANT"] = family["GWAS_FDR_BH"] <= alpha
    tag_tests = tested.merge(family.rename(columns={"SNP": "TESTED_SNP"}), on="TESTED_SNP")

    n_surrogate = int((tag_tests["TESTED_SNP"] != tag_tests["TAG_SNP"]).sum())
    n_sig_tags = int(tag_tests["BH_SIGNIFICANT"].sum())
    print(
        f"[overlap_dna_asm] {trait} Stage 1: {len(tag_tests):,}/{n_tags:,} TAG_SNPs tested "
        f"({len(family):,} unique SNPs in the BH family; {n_surrogate:,} tags via an LD-proxy "
        f"surrogate) | {n_sig_tags:,} BH-significant at FDR <= {alpha}."
    )

    # Stage 2: expand each tested tag to all its mapping rows (proxies kept even if absent from GWAS)
    overlap_df = mapping_df.merge(tag_tests, on="TAG_SNP", how="inner")
    query_gwas = gwas_df.rename(
        columns={"SNP": "QUERY_SNP", "GWAS_P": "QUERY_GWAS_P", "GWAS_BETA": "QUERY_GWAS_BETA"}
    )
    overlap_df = overlap_df.merge(query_gwas, on="QUERY_SNP", how="left")
    overlap_df["QUERY_IN_GWAS"] = overlap_df["QUERY_SNP"].isin(gwas_df["SNP"])
    overlap_df["TRAIT"] = trait
    overlap_df = overlap_df[OUTPUT_COLUMNS].sort_values(
        by=["CPG_ID", "TAG_SNP", "RELATIONSHIP", "QUERY_SNP"], ascending=[True, True, False, True]
    )
    sig_df = overlap_df[overlap_df["BH_SIGNIFICANT"]]

    for path in (output_all, output_sig):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    overlap_df.to_csv(output_all, sep="\t", index=False, na_rep="NA")
    sig_df.to_csv(output_sig, sep="\t", index=False, na_rep="NA")

    n_absent = int((~sig_df["QUERY_IN_GWAS"]).sum()) if not sig_df.empty else 0
    print(
        f"[overlap_dna_asm] {trait} Stage 2: {len(overlap_df):,} rows of tested tags -> {output_all}; "
        f"{len(sig_df):,} rows of BH-significant tags ({sig_df['QUERY_SNP'].nunique():,} SNPs, "
        f"{n_absent:,} rows whose SNP is absent from the GWAS, "
        f"{sig_df['CPG_ID'].nunique():,} CpGs) -> {output_sig}"
    )
    return sig_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Overlap DNA ASM tag SNPs with one GWAS (BH FDR), then expand significant tags to their LD proxies."
    )
    parser.add_argument("--mapping-table", required=True, help="DNA_ASM_mapping_table.tsv")
    parser.add_argument("--gwas", required=True, help="Raw GWAS summary statistics file")
    parser.add_argument("--trait", required=True, help="Trait identifier")
    parser.add_argument("--alpha", type=float, default=0.05, help="BH FDR level (default: 0.05)")
    parser.add_argument(
        "--tag-fallback",
        choices=TAG_FALLBACKS,
        default="none",
        help="Test for tags absent from the GWAS: skip them (none) or use their highest-LD proxy",
    )
    parser.add_argument("--output-all", required=True, help="Output TSV of all rows of tested tags")
    parser.add_argument("--output-sig", required=True, help="Output TSV of rows of BH-significant tags")
    args = parser.parse_args()
    run_overlap_dna_asm(
        args.mapping_table,
        args.gwas,
        args.trait,
        args.output_all,
        args.output_sig,
        args.alpha,
        args.tag_fallback,
    )


if __name__ == "__main__":
    main()
