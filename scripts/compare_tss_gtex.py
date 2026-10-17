#!/usr/bin/env python3
"""
Step 5d: Compare Proximal CpG TSS Genes (< 5 kb) with GTEx v10 Brain eQTL Target Genes.

Logic:
1. Nearest TSS Filtering:
   Deduplicates `dna_asm_{trait}_cpg_gene_tss_matches.tsv` by (CPG_ID, TAG_SNP, QUERY_SNP, GENE_SYMBOL)
   after sorting on DISTANCE_TO_TSS ascending, retaining only the single closest TSS per gene.
2. Outer Join & Concordance Classification:
   Merges deduplicated TSS matches with aggregated GTEx brain eQTL hits on QUERY_SNP and assigns:
   - IS_CONCORDANT: True when a valid TSS_GENE_SYMBOL (< 5 kb) matches GTEX_GENE_SYMBOL (case-insensitive).
   - MATCH_CATEGORY:
     * "FULL_CONCORDANCE": TSS_GENE_SYMBOL == GTEX_GENE_SYMBOL
     * "DISCORDANT_GENES": Both proximal TSS gene (< 5 kb) and brain eQTL gene exist, but differ
     * "TSS_ONLY": Proximal TSS gene (< 5 kb) exists, but no GTEx brain eQTL hit
     * "GTEX_ONLY": GTEx brain eQTL hit exists, but no gene TSS within 5 kb of CpG
     * "NONE_WITHIN_5KB": Neither proximal TSS gene nor GTEx brain eQTL hit
3. Sorting & Export:
   Exports full comparison table (`dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv`) and
   filtered high-priority concordant summary (`dna_asm_{trait}_concordant_genes_summary.tsv`).
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

OUTPUT_COLUMNS = [
    "TRAIT",
    "CPG_ID",
    "TAG_SNP",
    "QUERY_SNP",
    "RELATIONSHIP",
    "TSS_GENE_SYMBOL",
    "DISTANCE_TO_TSS",
    "GTEX_GENE_SYMBOL",
    "IS_CONCORDANT",
    "MATCH_CATEGORY",
    "GTEX_TISSUES",
    "N_BRAIN_TISSUES",
    "MIN_GTEX_PVALUE",
    "MEAN_GTEX_SLOPE",
    "GWAS_PVALUE",
    "GWAS_BETA",
]


def infer_trait_name(filepath: str) -> str:
    base = os.path.basename(filepath)
    m = re.match(r"^dna_asm_(.+)_cpg_gene_tss_matches\.tsv$", base)
    if m:
        return m.group(1)
    parts = os.path.normpath(filepath).split(os.sep)
    if len(parts) >= 2 and parts[-2] not in ("methylation", "results"):
        return parts[-2]
    return "UNKNOWN_TRAIT"


def run_compare_tss_gtex(
    tss_path: str,
    gtex_path: str,
    output_full_path: str,
    output_concordant_path: str,
    trait_name: str = None,
) -> tuple:
    trait = trait_name or infer_trait_name(tss_path)

    tss_df = pd.read_csv(tss_path, sep="\t", low_memory=False)
    gtex_df = pd.read_csv(gtex_path, sep="\t", low_memory=False)

    # Standardize column names from TSS table
    if "GENE_SYMBOL" in tss_df.columns and "TSS_GENE_SYMBOL" not in tss_df.columns:
        tss_df = tss_df.rename(columns={"GENE_SYMBOL": "TSS_GENE_SYMBOL"})
    if "GWAS_P" in tss_df.columns and "GWAS_PVALUE" not in tss_df.columns:
        tss_df = tss_df.rename(columns={"GWAS_P": "GWAS_PVALUE"})

    tss_df["DISTANCE_TO_TSS"] = pd.to_numeric(
        tss_df["DISTANCE_TO_TSS"], errors="coerce"
    ).astype("Int64")
    tss_df["GWAS_PVALUE"] = pd.to_numeric(tss_df["GWAS_PVALUE"], errors="coerce")
    tss_df["GWAS_BETA"] = pd.to_numeric(tss_df["GWAS_BETA"], errors="coerce")

    for col in ("CPG_ID", "TAG_SNP", "QUERY_SNP", "RELATIONSHIP", "TSS_GENE_SYMBOL"):
        if col in tss_df.columns:
            tss_df[col] = tss_df[col].astype(str).str.strip()

    # 1. Nearest TSS Filtering: sort by DISTANCE_TO_TSS ascending and keep head(1) per (CPG_ID, TAG_SNP, QUERY_SNP, upper(TSS_GENE_SYMBOL))
    tss_df["_GENE_UPPER"] = tss_df["TSS_GENE_SYMBOL"].str.upper()
    tss_dedup = (
        tss_df.sort_values(by=["DISTANCE_TO_TSS"], ascending=True, na_position="last")
        .drop_duplicates(
            subset=["CPG_ID", "TAG_SNP", "QUERY_SNP", "RELATIONSHIP", "_GENE_UPPER"],
            keep="first",
        )
        .drop(columns=["_GENE_UPPER"])
        .copy()
    )
    print(
        f"[compare_tss_gtex] {trait}: Deduplicated TSS rows from {len(tss_df):,} -> {len(tss_dedup):,} "
        f"(nearest TSS per CPG_ID x QUERY_SNP x GENE_SYMBOL)."
    )

    # 2. Clean GTEx table
    gtex_cols = [
        "QUERY_SNP",
        "GTEX_GENE_SYMBOL",
        "GTEX_TISSUES",
        "N_BRAIN_TISSUES",
        "MIN_GTEX_PVALUE",
        "MEAN_GTEX_SLOPE",
    ]
    if gtex_df.empty:
        gtex_clean = pd.DataFrame(columns=gtex_cols)
    else:
        gtex_clean = gtex_df.copy()
        gtex_clean["QUERY_SNP"] = gtex_clean["QUERY_SNP"].astype(str).str.strip()
        gtex_clean["GTEX_GENE_SYMBOL"] = (
            gtex_clean["GTEX_GENE_SYMBOL"].astype(str).str.strip()
        )
        gtex_clean["N_BRAIN_TISSUES"] = pd.to_numeric(
            gtex_clean["N_BRAIN_TISSUES"], errors="coerce"
        ).astype("Int64")
        gtex_clean["MIN_GTEX_PVALUE"] = pd.to_numeric(
            gtex_clean["MIN_GTEX_PVALUE"], errors="coerce"
        )
        gtex_clean["MEAN_GTEX_SLOPE"] = pd.to_numeric(
            gtex_clean["MEAN_GTEX_SLOPE"], errors="coerce"
        )
        gtex_clean = gtex_clean[gtex_cols]

    # 3. Outer join on QUERY_SNP
    merged = pd.merge(tss_dedup, gtex_clean, on="QUERY_SNP", how="outer")
    merged["TRAIT"] = trait

    # 4. Vectorized concordance and category classification
    tss_sym = merged["TSS_GENE_SYMBOL"].fillna("").astype(str).str.strip()
    gtex_sym = merged["GTEX_GENE_SYMBOL"].fillna("").astype(str).str.strip()

    tss_upper = tss_sym.str.upper()
    gtex_upper = gtex_sym.str.upper()

    invalid_tokens = {"", "NA", "NAN", "NONE", "<NA>", "NONE_WITHIN_5KB"}
    has_tss = ~tss_upper.isin(invalid_tokens)
    has_gtex = ~gtex_upper.isin(invalid_tokens)

    is_concordant = has_tss & has_gtex & (tss_upper == gtex_upper)
    merged["IS_CONCORDANT"] = is_concordant

    conditions = [
        is_concordant,
        has_tss & has_gtex & (~is_concordant),
        has_tss & (~has_gtex),
        (~has_tss) & has_gtex,
    ]
    choices = [
        "FULL_CONCORDANCE",
        "DISCORDANT_GENES",
        "TSS_ONLY",
        "GTEX_ONLY",
    ]
    merged["MATCH_CATEGORY"] = np.select(
        conditions, choices, default="NONE_WITHIN_5KB"
    )

    # Restore clean NA representations where GTEx or TSS was absent
    merged["GTEX_GENE_SYMBOL"] = np.where(has_gtex, gtex_sym, pd.NA)
    merged["TSS_GENE_SYMBOL"] = np.where(
        tss_upper.isin({"", "NA", "NAN", "NONE", "<NA>"}),
        "NONE_WITHIN_5KB",
        tss_sym,
    )
    merged["DISTANCE_TO_TSS"] = pd.to_numeric(
        merged["DISTANCE_TO_TSS"], errors="coerce"
    ).astype("Int64")
    merged["N_BRAIN_TISSUES"] = pd.to_numeric(
        merged["N_BRAIN_TISSUES"], errors="coerce"
    ).astype("Int64")

    # 5. Sort by IS_CONCORDANT descending, MATCH_CATEGORY rank, lowest GWAS_PVALUE, shortest DISTANCE_TO_TSS
    cat_rank_map = {
        "FULL_CONCORDANCE": 0,
        "DISCORDANT_GENES": 1,
        "GTEX_ONLY": 2,
        "TSS_ONLY": 3,
        "NONE_WITHIN_5KB": 4,
    }
    merged["_cat_rank"] = merged["MATCH_CATEGORY"].map(cat_rank_map).fillna(5)
    merged = (
        merged.sort_values(
            by=["IS_CONCORDANT", "_cat_rank", "GWAS_PVALUE", "DISTANCE_TO_TSS"],
            ascending=[False, True, True, True],
            na_position="last",
        )
        .drop(columns=["_cat_rank"])
        .reset_index(drop=True)
    )

    full_df = merged[OUTPUT_COLUMNS]
    concordant_df = full_df[full_df["IS_CONCORDANT"] == True].copy()

    os.makedirs(os.path.dirname(output_full_path) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(output_concordant_path) or ".", exist_ok=True)

    full_df.to_csv(output_full_path, sep="\t", index=False, na_rep="NA")
    concordant_df.to_csv(output_concordant_path, sep="\t", index=False, na_rep="NA")

    cat_counts = full_df["MATCH_CATEGORY"].value_counts().to_dict()
    n_conc_genes = concordant_df["TSS_GENE_SYMBOL"].nunique() if not concordant_df.empty else 0
    n_conc_cpgs = concordant_df["CPG_ID"].nunique() if not concordant_df.empty else 0
    print(
        f"[compare_tss_gtex] {trait} Summary: "
        f"FULL_CONCORDANCE={cat_counts.get('FULL_CONCORDANCE', 0):,} rows "
        f"({n_conc_genes:,} unique concordant genes across {n_conc_cpgs:,} CpGs) | "
        f"DISCORDANT_GENES={cat_counts.get('DISCORDANT_GENES', 0):,} | "
        f"GTEX_ONLY={cat_counts.get('GTEX_ONLY', 0):,} | "
        f"TSS_ONLY={cat_counts.get('TSS_ONLY', 0):,}"
    )
    print(f"[compare_tss_gtex] Full comparison saved to: {output_full_path}")
    print(f"[compare_tss_gtex] Concordant summary saved to: {output_concordant_path}")
    return full_df, concordant_df


def main() -> None:
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        tss_path = smk.input.tss
        gtex_path = smk.input.gtex
        output_full = smk.output.full
        output_concordant = smk.output.concordant
        trait = getattr(smk.wildcards, "trait", None)
        run_compare_tss_gtex(tss_path, gtex_path, output_full, output_concordant, trait)
        return

    parser = argparse.ArgumentParser(
        description="Compare CpG proximal TSS genes (< 5 kb) with GTEx brain eQTL target genes."
    )
    parser.add_argument("--tss", required=True, help="Input dna_asm_{trait}_cpg_gene_tss_matches.tsv")
    parser.add_argument("--gtex", required=True, help="Input dna_asm_{trait}_gtex_brain_eqtls.tsv")
    parser.add_argument(
        "--output-full",
        required=True,
        help="Output dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv",
    )
    parser.add_argument(
        "--output-concordant",
        required=True,
        help="Output dna_asm_{trait}_concordant_genes_summary.tsv",
    )
    parser.add_argument("--trait", default=None, help="Trait identifier")
    args = parser.parse_args()

    run_compare_tss_gtex(
        args.tss,
        args.gtex,
        args.output_full,
        args.output_concordant,
        trait_name=args.trait,
    )


if __name__ == "__main__":
    main()

