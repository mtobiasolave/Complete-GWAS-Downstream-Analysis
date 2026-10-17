# This script processes the DNA_ASM_Tag_LD_snps.xlsx file to create a mapping of CpG IDs to their associated tag SNPs and query LD SNPs, along with their relationships (TAG or LD_PROXY). The output is saved as a TSV file for further analysis.
# The output is DNA_ASM_mapping_table.tsv.

import pandas as pd
import numpy as np

def build_dna_asm_mapping(excel_file, output_tsv, nrows=None):
    # Load columns: X (IlmnID.1), Y (SNP.1), AC (SNPstagged), AD (tagSNP)
    # usecols by name matching the Excel sheet 'All'
    df = pd.read_excel(
        excel_file,
        sheet_name='All',
        usecols=['IlmnID.1', 'SNP.1', 'SNPstagged', 'tagSNP'],
        nrows=nrows
    )

    records = []

    for _, row in df.iterrows():
        cpg = str(row['IlmnID.1']).strip() if pd.notna(row['IlmnID.1']) else None
        snp_y = str(row['SNP.1']).strip() if pd.notna(row['SNP.1']) else None
        tagged = str(row['SNPstagged']).strip() if pd.notna(row['SNPstagged']) else None
        tag_snp = str(row['tagSNP']).strip() if pd.notna(row['tagSNP']) else None

        # Rule 2: Skip row if both tagSNP and SNPstagged are missing
        if tagged is None and tag_snp is None:
            continue

        # Case A: 'únic SNP'
        if tagged and tagged.lower() == 'únic snp':
            if snp_y:
                records.append({
                    'CPG_ID': cpg,
                    'TAG_SNP': snp_y,
                    'QUERY_SNP': snp_y,
                    'RELATIONSHIP': 'TAG'
                })
            continue

        # Case B: Standard tagged SNP entries
        if tagged:
            assigned_tag = tag_snp if tag_snp else snp_y
            rel = 'TAG' if tagged == assigned_tag else 'LD_PROXY'
            records.append({
                'CPG_ID': cpg,
                'TAG_SNP': assigned_tag,
                'QUERY_SNP': tagged,
                'RELATIONSHIP': rel
            })

    result_df = pd.DataFrame(records)

    # Rule 3 guarantee: Ensure every (CPG_ID, TAG_SNP) contains a self-referencing TAG row
    self_tags = []
    for (cpg, tag), group in result_df.groupby(['CPG_ID', 'TAG_SNP']):
        if not ((group['QUERY_SNP'] == tag) & (group['RELATIONSHIP'] == 'TAG')).any():
            self_tags.append({
                'CPG_ID': cpg,
                'TAG_SNP': tag,
                'QUERY_SNP': tag,
                'RELATIONSHIP': 'TAG'
            })
    
    if self_tags:
        result_df = pd.concat([result_df, pd.DataFrame(self_tags)], ignore_index=True)

    # Deduplicate exact identical rows (e.g. repeated across tissues)
    result_df = result_df.drop_duplicates()

    # Sort for clarity
    result_df = result_df.sort_values(by=['CPG_ID', 'TAG_SNP', 'RELATIONSHIP'], ascending=[True, True, False])

    # Save to TSV
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Generated {len(result_df)} unique mapping rows saved to {output_tsv}")
    return result_df

# For test run on first 100 rows:
# build_dna_asm_mapping('DNA_ASM_Tag_LD_snps.xlsx', 'DNA_ASM_mapping_sample_100.tsv', nrows=100)

# For the entire file:
build_dna_asm_mapping('DNA_ASM_Tag_LD_snps.xlsx', './DNA_ASM_mapping_table.tsv')