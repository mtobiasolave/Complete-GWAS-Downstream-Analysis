#!/usr/bin/env python3
"""
Steps A5 & B3: Cross-Trait Pleiotropy Summary of BH-Significant Methylation SNPs.

A SNP is "present" in a trait when it is BH-significant (FDR <= 5%) in that trait's
GWAS overlap (Step A1 / B1). No recurrence filter is applied here: the table keeps
every SNP significant in >= 1 trait, and N_TRAITS is filtered interactively in app.py.

--mode dna (inputs: results/dna/{trait}_concordance.tsv, already BH-significant only)
    SNP, N_TRAITS, TRAITS_PRESENT, CONCORDANT, CONCORDANT_GENES, CPG_GENES_5KB,
    GTEX_BRAIN_GENES, CPG_IDS, GTEX_TISSUES, MIN_FDR_BH

--mode rna (inputs: results/rna/{trait}_overlap_bh.tsv, filtered on BH_SIGNIFICANT;
            --genes: Step B2 results/rna/rna_sig_genes.tsv)
    SNP, GENES, N_TRAITS, TRAITS_PRESENT, MIN_FDR_BH

Usage:
    python scripts/pleiotropy.py --mode dna --inputs results/dna/*_concordance.tsv \
        --output results/dna/pleiotropy_dna.tsv
"""

import argparse
import os

import numpy as np
import pandas as pd

from common_utils import join_unique_clean

DNA_COLUMNS = [
    "SNP",
    "N_TRAITS",
    "TRAITS_PRESENT",
    "CONCORDANT",
    "CONCORDANT_GENES",
    "CPG_GENES_5KB",
    "GTEX_BRAIN_GENES",
    "CPG_IDS",
    "GTEX_TISSUES",
    "MIN_FDR_BH",
]
RNA_COLUMNS = ["SNP", "GENES", "N_TRAITS", "TRAITS_PRESENT", "MIN_FDR_BH"]


def is_true(series: pd.Series) -> pd.Series:
    """Parse boolean columns read back from TSV ('True'/'False'/NA)."""
    return series.astype(str).str.strip().str.lower().isin(["true", "1"])


def summarize_dna(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=DNA_COLUMNS)

    df = df.copy()
    df["_CONC_GENE"] = np.where(is_true(df["IS_CONCORDANT"]), df["TSS_GENE_SYMBOL"], pd.NA)
    df["GWAS_FDR_BH"] = pd.to_numeric(df["GWAS_FDR_BH"], errors="coerce")

    out = (
        df.groupby("QUERY_SNP", sort=False)
        .agg(
            N_TRAITS=("TRAIT", "nunique"),
            TRAITS_PRESENT=("TRAIT", join_unique_clean),
            CONCORDANT_GENES=("_CONC_GENE", join_unique_clean),
            CPG_GENES_5KB=("TSS_GENE_SYMBOL", join_unique_clean),
            GTEX_BRAIN_GENES=("GTEX_GENE_SYMBOL", join_unique_clean),
            CPG_IDS=("CPG_ID", join_unique_clean),
            GTEX_TISSUES=("GTEX_TISSUES", join_unique_clean),
            MIN_FDR_BH=("GWAS_FDR_BH", "min"),
        )
        .reset_index()
        .rename(columns={"QUERY_SNP": "SNP"})
    )
    out["CONCORDANT"] = out["CONCORDANT_GENES"] != "NA"
    out = out.sort_values(
        by=["N_TRAITS", "CONCORDANT", "MIN_FDR_BH"], ascending=[False, False, True]
    )
    return out[DNA_COLUMNS].reset_index(drop=True)


def summarize_rna(df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    if not df.empty:
        df = df[is_true(df["BH_SIGNIFICANT"])].copy()
    if df.empty:
        return pd.DataFrame(columns=RNA_COLUMNS)

    df["GWAS_FDR_BH"] = pd.to_numeric(df["GWAS_FDR_BH"], errors="coerce")
    out = (
        df.groupby("SNP", sort=False)
        .agg(
            N_TRAITS=("TRAIT", "nunique"),
            TRAITS_PRESENT=("TRAIT", join_unique_clean),
            MIN_FDR_BH=("GWAS_FDR_BH", "min"),
        )
        .reset_index()
        .sort_values(by=["N_TRAITS", "MIN_FDR_BH"], ascending=[False, True])
    )
    out["GENES"] = out["SNP"].map(dict(zip(genes_df["SNP"], genes_df["GENES"]))).fillna("NA")
    return out[RNA_COLUMNS].reset_index(drop=True)


def run_pleiotropy(mode: str, inputs: list, output_path: str, genes_path: str = None) -> pd.DataFrame:
    frames = [pd.read_csv(path, sep="\t", low_memory=False) for path in inputs]
    frames = [f for f in frames if not f.empty]
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    if mode == "dna":
        summary = summarize_dna(combined)
    else:
        genes_df = pd.read_csv(genes_path, sep="\t", usecols=["SNP", "GENES"], dtype=str, keep_default_na=False)
        summary = summarize_rna(combined, genes_df)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    summary.to_csv(output_path, sep="\t", index=False, na_rep="NA")

    n_multi = int((summary["N_TRAITS"] >= 2).sum()) if not summary.empty else 0
    print(
        f"[pleiotropy:{mode}] {len(inputs)} traits | {len(summary):,} BH-significant SNPs "
        f"({n_multi:,} in >= 2 traits) -> {output_path}"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cross-trait pleiotropy summary of BH-significant DNA ASM or RNA m6A SNPs."
    )
    parser.add_argument("--mode", choices=["dna", "rna"], required=True, help="Methylation track")
    parser.add_argument("--inputs", nargs="+", required=True, help="Per-trait result TSVs")
    parser.add_argument("--output", "-o", required=True, help="Output pleiotropy TSV")
    parser.add_argument("--genes", help="Step B2 rna_sig_genes.tsv (required for --mode rna)")
    args = parser.parse_args()
    if args.mode == "rna" and not args.genes:
        parser.error("--genes is required for --mode rna")
    run_pleiotropy(args.mode, args.inputs, args.output, args.genes)


if __name__ == "__main__":
    main()
