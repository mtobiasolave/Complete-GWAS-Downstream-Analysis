#!/usr/bin/env python3
"""
Step 5a: Intersect GWAS Summary Statistics with DNA ASM Tag SNPs and LD Proxies.

Reads data/methylation/DNA_ASM_mapping_table.tsv (CPG_ID, TAG_SNP, QUERY_SNP, RELATIONSHIP)
and performs an inner join (GWAS_SNP == QUERY_SNP) against raw GWAS summary statistics
without any p-value filtering, retaining GWAS P-value and Effect Size (BETA / ln(OR))
columns alongside CPG_ID, TAG_SNP, QUERY_SNP, and RELATIONSHIP.

Usage:
    python scripts/overlap_dna_asm.py \
        --mapping-table data/methylation/DNA_ASM_mapping_table.tsv \
        --gwas data/raw_data/BipolarDisorder2024_gwas.txt ... \
        --output results/methylation/dna_asm_gwas_overlap.tsv
"""

import argparse
import os
from functools import reduce

import numpy as np
import pandas as pd

from common_utils import detect_column, extract_matching_gwas_df


def parse_trait_name(gwas_path: str) -> str:
    base = os.path.basename(gwas_path)
    for suffix in (
        "_gwas.txt.gz",
        "_gwas.tsv.gz",
        "_gwas.txt",
        "_gwas.tsv",
        ".txt.gz",
        ".tsv.gz",
        ".txt",
        ".tsv",
    ):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


def run_overlap_dna_asm(mapping_table_path: str, gwas_files: list, output_path: str) -> pd.DataFrame:
    if isinstance(gwas_files, str):
        gwas_files = [gwas_files]

    mapping_df = pd.read_csv(mapping_table_path, sep="\t", dtype=str)
    required_cols = ["CPG_ID", "TAG_SNP", "QUERY_SNP", "RELATIONSHIP"]
    missing = [c for c in required_cols if c not in mapping_df.columns]
    if missing:
        raise ValueError(f"Mapping table {mapping_table_path} is missing required columns: {missing}")

    for col in required_cols:
        mapping_df[col] = mapping_df[col].astype(str).str.strip()

    mapping_df = mapping_df[required_cols].drop_duplicates()
    cand_snps = set(mapping_df["QUERY_SNP"].dropna().unique())
    print(
        f"[overlap_dna_asm] Loaded {len(mapping_df):,} mapping rows "
        f"({len(cand_snps):,} unique QUERY_SNPs)."
    )

    single_gwas_mode = len(gwas_files) == 1
    trait_dfs = []
    matched_snp_union = set()

    for gwas_path in sorted(gwas_files):
        trait_name = parse_trait_name(gwas_path)
        df_raw = extract_matching_gwas_df(gwas_path, cand_snps)

        snp_col = detect_column(df_raw.columns, "SNP")
        p_col = detect_column(df_raw.columns, "P")
        beta_col = detect_column(df_raw.columns, "BETA")
        or_col = detect_column(df_raw.columns, "OR")

        if not snp_col or not p_col:
            raise ValueError(
                f"Could not detect SNP ({snp_col}) or P ({p_col}) column in {gwas_path}. "
                f"Columns found: {list(df_raw.columns)}"
            )

        df_raw[snp_col] = df_raw[snp_col].astype(str).str.strip()
        p_vals = pd.to_numeric(df_raw[p_col], errors="coerce")

        if beta_col:
            beta_vals = pd.to_numeric(df_raw[beta_col], errors="coerce")
        elif or_col:
            or_vals = pd.to_numeric(df_raw[or_col], errors="coerce")
            beta_vals = np.where(or_vals > 0, np.log(or_vals), np.nan)
        else:
            beta_vals = pd.Series(np.nan, index=df_raw.index)

        p_out_col = "GWAS_P" if single_gwas_mode else f"P_{trait_name}"
        beta_out_col = "GWAS_BETA" if single_gwas_mode else f"BETA_{trait_name}"

        trait_sub = (
            pd.DataFrame({
                "QUERY_SNP": df_raw[snp_col],
                p_out_col: p_vals,
                beta_out_col: beta_vals,
            })
            .dropna(subset=["QUERY_SNP"])
            .drop_duplicates(subset=["QUERY_SNP"])
        )

        matched_snp_union.update(trait_sub["QUERY_SNP"].unique())
        trait_dfs.append(trait_sub)
        print(
            f"[overlap_dna_asm] {trait_name}: matched {len(trait_sub):,} unique QUERY_SNPs "
            f"(P col='{p_col}', Effect col='{beta_col or or_col}', "
            f"converted OR->ln(OR)={bool(not beta_col and or_col)})."
        )

    if not trait_dfs or not matched_snp_union:
        gwas_cols = []
        for t_df in trait_dfs:
            gwas_cols.extend([c for c in t_df.columns if c != "QUERY_SNP"])
        overlap_df = pd.DataFrame(columns=required_cols + gwas_cols)
    elif single_gwas_mode:
        overlap_df = pd.merge(mapping_df, trait_dfs[0], on="QUERY_SNP", how="inner")
    else:
        combined_gwas = reduce(
            lambda left, right: pd.merge(left, right, on="QUERY_SNP", how="outer"),
            trait_dfs,
        )
        overlap_df = pd.merge(mapping_df, combined_gwas, on="QUERY_SNP", how="inner")

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    overlap_df.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    print(f"[overlap_dna_asm] Saved {len(overlap_df):,} overlapping rows to {output_path}")
    return overlap_df


def main() -> None:
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        run_overlap_dna_asm(smk.input.mapping_table, smk.input.gwas, smk.output[0])
        return

    parser = argparse.ArgumentParser(
        description="Intersect GWAS summary statistics with DNA ASM Tag SNPs and LD proxies."
    )
    parser.add_argument(
        "--mapping-table",
        required=True,
        help="Path to data/methylation/DNA_ASM_mapping_table.tsv",
    )
    parser.add_argument(
        "--gwas",
        nargs="+",
        required=True,
        help="One or more raw GWAS summary statistics files",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output TSV path (results/methylation/dna_asm_gwas_overlap.tsv)",
    )
    args = parser.parse_args()
    run_overlap_dna_asm(args.mapping_table, args.gwas, args.output)


if __name__ == "__main__":
    main()
