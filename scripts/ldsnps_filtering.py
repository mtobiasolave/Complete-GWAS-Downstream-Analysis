"""Process DNA_ASM_Tag_LD_snps.xlsx into the long-format DNA_ASM_mapping_table.tsv reference.

Retained as a standalone reference preprocessing utility.
"""

import argparse
import pandas as pd


def build_dna_asm_mapping(excel_file: str, output_tsv: str, nrows: int = None) -> pd.DataFrame:
    """Build CpG-to-TAG_SNP/QUERY_SNP mapping table from the source Excel workbook."""
    df = pd.read_excel(
        excel_file,
        sheet_name="All",
        usecols=["IlmnID.1", "SNP.1", "SNPstagged", "tagSNP"],
        nrows=nrows,
    )

    records = []
    for row in df.itertuples(index=False):
        cpg_raw, snp_y_raw, tagged_raw, tag_snp_raw = row
        cpg = str(cpg_raw).strip() if pd.notna(cpg_raw) else None
        snp_y = str(snp_y_raw).strip() if pd.notna(snp_y_raw) else None
        tagged = str(tagged_raw).strip() if pd.notna(tagged_raw) else None
        tag_snp = str(tag_snp_raw).strip() if pd.notna(tag_snp_raw) else None

        if tagged is None and tag_snp is None:
            continue

        if tagged and tagged.lower() == "únic snp":
            if snp_y:
                records.append({
                    "CPG_ID": cpg,
                    "TAG_SNP": snp_y,
                    "QUERY_SNP": snp_y,
                    "RELATIONSHIP": "TAG",
                })
            continue

        if tagged:
            assigned_tag = tag_snp if tag_snp else snp_y
            rel = "TAG" if tagged == assigned_tag else "LD_PROXY"
            records.append({
                "CPG_ID": cpg,
                "TAG_SNP": assigned_tag,
                "QUERY_SNP": tagged,
                "RELATIONSHIP": rel,
            })

    result_df = pd.DataFrame(records)
    if result_df.empty:
        return result_df

    # Guarantee every (CPG_ID, TAG_SNP) pair includes a self-referencing TAG row
    has_self_tag = (result_df["QUERY_SNP"] == result_df["TAG_SNP"]) & (
        result_df["RELATIONSHIP"] == "TAG"
    )
    pairs_with_self = set(
        zip(result_df.loc[has_self_tag, "CPG_ID"], result_df.loc[has_self_tag, "TAG_SNP"])
    )
    all_pairs = (
        result_df[["CPG_ID", "TAG_SNP"]].drop_duplicates().itertuples(index=False, name=None)
    )
    self_tags = [
        {"CPG_ID": cpg, "TAG_SNP": tag, "QUERY_SNP": tag, "RELATIONSHIP": "TAG"}
        for cpg, tag in all_pairs
        if (cpg, tag) not in pairs_with_self
    ]

    if self_tags:
        result_df = pd.concat([result_df, pd.DataFrame(self_tags)], ignore_index=True)

    result_df = result_df.drop_duplicates().sort_values(
        by=["CPG_ID", "TAG_SNP", "RELATIONSHIP"], ascending=[True, True, False]
    )
    result_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Generated {len(result_df)} unique mapping rows saved to {output_tsv}")
    return result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build DNA_ASM_mapping_table.tsv from DNA_ASM_Tag_LD_snps.xlsx."
    )
    parser.add_argument(
        "--input", default="DNA_ASM_Tag_LD_snps.xlsx", help="Path to source Excel file."
    )
    parser.add_argument(
        "--output",
        default="data/methylation/DNA_ASM_mapping_table.tsv",
        help="Output TSV path.",
    )
    parser.add_argument("--nrows", type=int, default=None, help="Optional row limit for testing.")
    args = parser.parse_args()
    build_dna_asm_mapping(args.input, args.output, nrows=args.nrows)