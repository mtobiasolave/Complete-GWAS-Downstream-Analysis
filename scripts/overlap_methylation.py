#!/usr/bin/env python3
"""
Step 5 Helper: Overlap Methylation SNP Lists (DNA_ASM, RNA_m6A) with Raw GWAS Summary Statistics.
Cross-platform (Windows & Linux Slurm HPC) replacement for shell awk/sort pipeline.

Usage:
    python scripts/overlap_methylation.py <snp_list> <gwas_file> <out_overlap_full> <out_overlap_ids>
"""

import os
import sys


def find_snp_column(header_line: str) -> int:
    """
    Identify 0-based SNP column index by priority:
    RSID preferred over generic SNP / MARKERNAME / ID.
    Returns -1 if no matching column is found.
    """
    cols = header_line.strip().split()
    rsid_col = -1
    fallback_col = -1
    for idx, col_name in enumerate(cols):
        h = col_name.strip().upper()
        if h == "RSID":
            rsid_col = idx
            break
        if h in ("SNP", "MARKERNAME", "ID") and fallback_col == -1:
            fallback_col = idx
    return rsid_col if rsid_col >= 0 else fallback_col


def main() -> None:
    if len(sys.argv) != 5:
        sys.stderr.write(
            "Usage: python scripts/overlap_methylation.py "
            "<snp_list> <gwas_file> <out_overlap_full> <out_overlap_ids>\n"
        )
        sys.exit(1)

    snp_list_path = sys.argv[1]
    gwas_path = sys.argv[2]
    overlap_full_path = sys.argv[3]
    overlap_ids_path = sys.argv[4]

    os.makedirs(os.path.dirname(overlap_full_path) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(overlap_ids_path) or ".", exist_ok=True)

    # 1. Load candidate methylation rsIDs into a set
    cand_snps = set()
    with open(snp_list_path, "r", encoding="utf-8", errors="replace") as f_snp:
        for line in f_snp:
            parts = line.strip().split()
            if parts:
                cand_snps.add(parts[0])

    # 2. Stream GWAS file, detect SNP column, and write matching rows
    matched_snps = set()
    with open(gwas_path, "r", encoding="utf-8", errors="replace") as f_gwas, open(
        overlap_full_path, "w", encoding="utf-8", newline="\n"
    ) as f_out:
        header_line = f_gwas.readline()
        if not header_line:
            sys.stderr.write(f"ERROR: Empty GWAS file: {gwas_path}\n")
            sys.exit(1)

        snp_col = find_snp_column(header_line)
        if snp_col < 0:
            sys.stderr.write(
                f"ERROR: Could not find valid SNP/RSID column in {gwas_path}\n"
            )
            sys.exit(1)

        f_out.write(header_line.rstrip("\r\n") + "\n")
        max_split = snp_col + 1

        if snp_col == 0:
            for line in f_gwas:
                parts = line.split(None, 1)
                if parts and parts[0] in cand_snps:
                    matched_snps.add(parts[0])
                    f_out.write(line.rstrip("\r\n") + "\n")
        else:
            for line in f_gwas:
                parts = line.split(None, max_split)
                if len(parts) > snp_col and parts[snp_col] in cand_snps:
                    snp_id = parts[snp_col]
                    matched_snps.add(snp_id)
                    f_out.write(line.rstrip("\r\n") + "\n")

    # 3. Write unique sorted rsIDs (single-column, no header)
    with open(overlap_ids_path, "w", encoding="utf-8", newline="\n") as f_ids:
        for snp_id in sorted(matched_snps):
            f_ids.write(f"{snp_id}\n")

    print(
        f"[overlap_methylation] {os.path.basename(snp_list_path)} x {os.path.basename(gwas_path)}: "
        f"{len(matched_snps):,} unique overlapping SNPs written."
    )


if __name__ == "__main__":
    main()

