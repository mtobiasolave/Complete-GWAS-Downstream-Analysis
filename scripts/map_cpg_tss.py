#!/usr/bin/env python3
"""
Step 5b: Map Overlapping DNA ASM CpGs to hg19 Coordinates and GENCODE v19 Gene TSS (< 5 kb).

1. Merges `cpg_coordinates_hg19.tsv` with `dna_asm_gwas_overlap.tsv` on `CPG_ID`.
2. Standardizes chromosome notation across tables (preserving `chr` prefix consistently).
3. Queries `gencode_v19_gene_tss.tsv` on the same chromosome, computing:
       DISTANCE_TO_TSS = |CPG_POS - TSS|
   and filtering for DISTANCE_TO_TSS <= tss_distance_cutoff (default 5,000 bp).
4. Retains multiple genes within 5 kb in tidy long format, and retains CpGs without any
   gene within 5 kb using GENE_SYMBOL = "NONE_WITHIN_5KB" and DISTANCE_TO_TSS = NA.

Output Columns:
    CPG_ID, CHR, CPG_POS, TAG_SNP, QUERY_SNP, RELATIONSHIP,
    GENE_SYMBOL, STRAND, TSS, DISTANCE_TO_TSS, BIOTYPE, <GWAS_P>, <GWAS_BETA>
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def standardize_chr(series: pd.Series) -> pd.Series:
    """Standardize chromosome labels to consistently include the 'chr' prefix (e.g., 'chr16')."""
    s = series.astype(str).str.strip()
    is_na = series.isna() | s.isin(["", "nan", "None", "NA", "<NA>"])
    normalized = s.apply(
        lambda val: val if val.lower().startswith("chr") else f"chr{val}"
    )
    # Ensure canonical lowercase 'chr' prefix
    normalized = normalized.str.replace(r"^CHR", "chr", regex=True)
    return normalized.where(~is_na, other=pd.NA)


def map_cpg_to_tss(
    overlap_path: str,
    cpg_coords_path: str,
    gene_tss_path: str,
    output_path: str,
    cutoff: int = 5000,
) -> pd.DataFrame:
    overlap_df = pd.read_csv(overlap_path, sep="\t", low_memory=False)
    cpg_df = pd.read_csv(cpg_coords_path, sep="\t", low_memory=False)
    tss_df = pd.read_csv(gene_tss_path, sep="\t", low_memory=False)

    base_overlap_cols = ["CPG_ID", "TAG_SNP", "QUERY_SNP", "RELATIONSHIP"]
    for col in base_overlap_cols:
        if col not in overlap_df.columns:
            raise ValueError(f"Overlap file {overlap_path} missing required column: {col}")
        overlap_df[col] = overlap_df[col].astype(str).str.strip()

    gwas_cols = [c for c in overlap_df.columns if c not in base_overlap_cols]

    # Clean and standardize CpG coordinate table
    cpg_df["CPG_ID"] = cpg_df["CPG_ID"].astype(str).str.strip()
    cpg_df["CHR"] = standardize_chr(cpg_df["CHR"])
    cpg_df["CPG_POS"] = pd.to_numeric(cpg_df["CPG_POS"], errors="coerce").astype("Int64")
    cpg_df = cpg_df[["CPG_ID", "CHR", "CPG_POS"]].drop_duplicates(subset=["CPG_ID"])

    # Clean and standardize GENCODE v19 TSS table
    tss_df["CHR"] = standardize_chr(tss_df["CHR"])
    tss_df["GENE_SYMBOL"] = tss_df["GENE_SYMBOL"].astype(str).str.strip()
    tss_df["STRAND"] = tss_df["STRAND"].astype(str).str.strip()
    tss_df["TSS"] = pd.to_numeric(tss_df["TSS"], errors="coerce").astype("Int64")
    tss_df["BIOTYPE"] = tss_df["BIOTYPE"].astype(str).str.strip()
    tss_df = tss_df.dropna(subset=["CHR", "TSS"]).drop_duplicates(
        subset=["CHR", "GENE_SYMBOL", "STRAND", "TSS", "BIOTYPE"]
    )

    # Extract unique CpGs from overlap table and merge coordinates
    unique_cpgs = pd.DataFrame({"CPG_ID": overlap_df["CPG_ID"].dropna().unique()})
    cpg_mapped = pd.merge(unique_cpgs, cpg_df, on="CPG_ID", how="left")

    # Vectorized per-chromosome TSS distance calculation on unique CpGs
    matched_chunks = []
    valid_cpgs = cpg_mapped.dropna(subset=["CHR", "CPG_POS"])

    for chrom, cpg_chr_group in valid_cpgs.groupby("CHR", sort=False):
        tss_chr_group = tss_df[tss_df["CHR"] == chrom]
        if tss_chr_group.empty:
            continue

        merged_chr = pd.merge(cpg_chr_group, tss_chr_group, on="CHR", how="inner")
        dist = (merged_chr["CPG_POS"] - merged_chr["TSS"]).abs()
        within_cutoff = merged_chr[dist <= cutoff].copy()
        if not within_cutoff.empty:
            within_cutoff["DISTANCE_TO_TSS"] = dist[dist <= cutoff].astype("Int64")
            matched_chunks.append(within_cutoff)

    cpg_tss_cols = [
        "CPG_ID",
        "CHR",
        "CPG_POS",
        "GENE_SYMBOL",
        "STRAND",
        "TSS",
        "DISTANCE_TO_TSS",
        "BIOTYPE",
    ]

    if matched_chunks:
        matched_cpg_tss = pd.concat(matched_chunks, ignore_index=True)[cpg_tss_cols]
        matched_cpg_tss = matched_cpg_tss.sort_values(
            by=["CPG_ID", "DISTANCE_TO_TSS", "GENE_SYMBOL"]
        )
    else:
        matched_cpg_tss = pd.DataFrame(columns=cpg_tss_cols)

    # Retain CpGs with no gene within cutoff as GENE_SYMBOL = "NONE_WITHIN_5KB" and DISTANCE_TO_TSS = NA
    matched_cpg_ids = set(matched_cpg_tss["CPG_ID"].unique())
    unmatched_cpgs = cpg_mapped[~cpg_mapped["CPG_ID"].isin(matched_cpg_ids)].copy()

    if not unmatched_cpgs.empty:
        unmatched_cpgs["GENE_SYMBOL"] = "NONE_WITHIN_5KB"
        unmatched_cpgs["STRAND"] = pd.NA
        unmatched_cpgs["TSS"] = pd.Series(pd.NA, index=unmatched_cpgs.index, dtype="Int64")
        unmatched_cpgs["DISTANCE_TO_TSS"] = pd.Series(
            pd.NA, index=unmatched_cpgs.index, dtype="Int64"
        )
        unmatched_cpgs["BIOTYPE"] = pd.NA
        unmatched_cpgs = unmatched_cpgs[cpg_tss_cols]
        all_cpg_tss = pd.concat([matched_cpg_tss, unmatched_cpgs], ignore_index=True)
    else:
        all_cpg_tss = matched_cpg_tss

    # Ensure nullable integer types are preserved after concat
    all_cpg_tss["CPG_POS"] = pd.to_numeric(all_cpg_tss["CPG_POS"], errors="coerce").astype("Int64")
    all_cpg_tss["TSS"] = pd.to_numeric(all_cpg_tss["TSS"], errors="coerce").astype("Int64")
    all_cpg_tss["DISTANCE_TO_TSS"] = pd.to_numeric(
        all_cpg_tss["DISTANCE_TO_TSS"], errors="coerce"
    ).astype("Int64")

    # Merge CpG-to-TSS annotations back with full SNP overlap table
    final_df = pd.merge(overlap_df, all_cpg_tss, on="CPG_ID", how="left")

    ordered_cols = [
        "CPG_ID",
        "CHR",
        "CPG_POS",
        "TAG_SNP",
        "QUERY_SNP",
        "RELATIONSHIP",
        "GENE_SYMBOL",
        "STRAND",
        "TSS",
        "DISTANCE_TO_TSS",
        "BIOTYPE",
    ] + gwas_cols

    final_df = final_df[ordered_cols]

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    final_df.to_csv(output_path, sep="\t", index=False, na_rep="NA")

    n_total_cpgs = len(unique_cpgs)
    n_cpgs_with_gene = len(matched_cpg_ids)
    n_unique_genes = (
        matched_cpg_tss["GENE_SYMBOL"].nunique() if not matched_cpg_tss.empty else 0
    )
    print(
        f"[map_cpg_tss] Mapped {n_total_cpgs:,} unique CpGs "
        f"({n_cpgs_with_gene:,} with >=1 gene TSS within {cutoff:,} bp; "
        f"{n_total_cpgs - n_cpgs_with_gene:,} NONE_WITHIN_5KB; "
        f"{n_unique_genes:,} unique proximal genes)."
    )
    print(f"[map_cpg_tss] Wrote {len(final_df):,} tidy rows to {output_path}")
    return final_df


def main() -> None:
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        overlap_path = smk.input.overlap
        cpg_coords_path = smk.input.cpg_coordinates
        gene_tss_path = smk.input.gene_tss
        output_path = smk.output[0]
        cutoff = int(getattr(smk.params, "tss_distance_cutoff", 5000))
        map_cpg_to_tss(overlap_path, cpg_coords_path, gene_tss_path, output_path, cutoff)
        return

    parser = argparse.ArgumentParser(
        description="Map overlapping DNA ASM CpG sites to hg19 coordinates and GENCODE v19 TSS (< 5 kb)."
    )
    parser.add_argument(
        "--overlap",
        required=True,
        help="Path to results/methylation/dna_asm_gwas_overlap.tsv",
    )
    parser.add_argument(
        "--cpg-coords",
        required=True,
        help="Path to data/annotations/cpg_coordinates_hg19.tsv",
    )
    parser.add_argument(
        "--gene-tss",
        required=True,
        help="Path to data/annotations/gencode_v19_gene_tss.tsv",
    )
    parser.add_argument(
        "--cutoff",
        type=int,
        default=5000,
        help="Maximum CpG-to-TSS distance in bp (default: 5000)",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output TSV path (results/methylation/dna_asm_cpg_gene_tss_matches.tsv)",
    )
    args = parser.parse_args()
    map_cpg_to_tss(
        args.overlap, args.cpg_coords, args.gene_tss, args.output, args.cutoff
    )


if __name__ == "__main__":
    main()

