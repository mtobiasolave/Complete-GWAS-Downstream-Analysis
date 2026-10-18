#!/usr/bin/env python3
"""
Step 5 Helper: Overlap Methylation SNP Lists (DNA_ASM, RNA_m6A) with Raw GWAS Summary Statistics.
Cross-platform (Windows & Linux Slurm HPC) replacement for shell awk/sort pipeline.

Usage:
    python scripts/overlap_methylation.py <snp_list> <gwas_file> <out_overlap_full> <out_overlap_ids>
"""

import os
import sys

from common_utils import extract_matching_gwas_rows


def main() -> None:
    if len(sys.argv) != 5:
        sys.stderr.write(
            "Usage: python scripts/overlap_methylation.py "
            "<snp_list> <gwas_file> <out_overlap_full> <out_overlap_ids>\n"
        )
        sys.exit(1)

    snp_list_path, gwas_path, overlap_full_path, overlap_ids_path = sys.argv[1:5]

    os.makedirs(os.path.dirname(overlap_full_path) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(overlap_ids_path) or ".", exist_ok=True)

    # 1. Load candidate methylation rsIDs into a set
    cand_snps = set()
    with open(snp_list_path, "r", encoding="utf-8", errors="replace") as f_snp:
        for line in f_snp:
            parts = line.strip().split()
            if parts:
                cand_snps.add(parts[0])

    # 2. Stream GWAS file (plain or .gz) and extract matching rows
    try:
        header_line, matched_lines, matched_snps = extract_matching_gwas_rows(gwas_path, cand_snps)
    except ValueError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        sys.exit(1)

    with open(overlap_full_path, "w", encoding="utf-8", newline="\n") as f_out:
        f_out.write(header_line + "\n")
        for line in matched_lines:
            f_out.write(line + "\n")

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
