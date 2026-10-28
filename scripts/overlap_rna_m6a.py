#!/usr/bin/env python3
"""
Step B1: Overlap the RNA m6A SNP List with One GWAS and Apply BH FDR (5%).

Streams the GWAS summary statistics against the RNA m6A SNP list (first column of
each line) and applies Benjamini-Hochberg FDR across the unique overlapping SNPs.

Output columns:
    TRAIT, SNP, GWAS_P, GWAS_BETA, GWAS_FDR_BH, BH_SIGNIFICANT

Usage:
    python scripts/overlap_rna_m6a.py --snp-list data/methylation/RNA_m6A_SNPlist.txt \
        --gwas data/raw_data/SUD2023_gwas.txt --trait SUD2023 --alpha 0.05 \
        --output results/rna/SUD2023_overlap_bh.tsv
"""

import argparse
import os

from common_utils import gwas_overlap_bh


def run_overlap_rna_m6a(snp_list_path: str, gwas_path: str, trait: str, output_path: str, alpha: float = 0.05):
    cand_snps = set()
    with open(snp_list_path, "r", encoding="utf-8", errors="replace") as f_snp:
        for line in f_snp:
            parts = line.strip().split()
            if parts:
                cand_snps.add(parts[0])
    print(f"[overlap_rna_m6a] {trait}: {len(cand_snps):,} RNA m6A candidate SNPs loaded.")

    overlap_df = gwas_overlap_bh(gwas_path, cand_snps, alpha=alpha)
    overlap_df.insert(0, "TRAIT", trait)

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    overlap_df.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    print(f"[overlap_rna_m6a] {trait}: {len(overlap_df):,} overlapping SNPs written to {output_path}")
    return overlap_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Overlap the RNA m6A SNP list with one GWAS and apply BH FDR."
    )
    parser.add_argument("--snp-list", required=True, help="RNA_m6A_SNPlist.txt")
    parser.add_argument("--gwas", required=True, help="Raw GWAS summary statistics file")
    parser.add_argument("--trait", required=True, help="Trait identifier")
    parser.add_argument("--alpha", type=float, default=0.05, help="BH FDR level (default: 0.05)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV path")
    args = parser.parse_args()
    run_overlap_rna_m6a(args.snp_list, args.gwas, args.trait, args.output, args.alpha)


if __name__ == "__main__":
    main()
