#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')  # Headless cluster rendering
import matplotlib.pyplot as plt
import seaborn as sns

COMPARISONS_DIR = "results/methylation/comparisons"
FIGURES_DIR = "results/figures"

def load_annotated_data(input_file):
    """Loads the annotated TSV using strict tab separation, skipping comments."""
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    return pd.read_csv(input_file, sep='\t', comment='#')

def extract_traits(df):
    """Extracts trait names based on FDR columns."""
    fdr_cols = [c for c in df.columns if c.startswith("FDR_")]
    return [c[len("FDR_"):] for c in fdr_cols]

def find_trait_columns(df, trait):
    """Dynamically locates P, FDR, and Effect columns for a trait."""
    p_col = f"P_{trait}" if f"P_{trait}" in df.columns else None
    fdr_col = f"FDR_{trait}" if f"FDR_{trait}" in df.columns else None

    eff_col = None
    eff_type = "BETA"
    candidates = [
        f"Effect_{trait} (BETA)",
        f"Effect_{trait} (OR)",
        f"Effect_{trait}"
    ]
    for c in candidates:
        if c in df.columns:
            eff_col = c
            if "(OR)" in c.upper():
                eff_type = "OR"
            break
            
    if not eff_col:
        for c in df.columns:
            if c.startswith(f"Effect_{trait}"):
                eff_col = c
                if "(OR)" in c.upper():
                    eff_type = "OR"
                break

    return p_col, fdr_col, eff_col, eff_type

def clean_trait_name(trait):
    """Strips redundant pipeline prefixes for cleaner Y-axis labels."""
    for prefix in ["m6A_", "DNA_ASM_", "RNA_m6A_", "RNAm6A_"]:
        if trait.startswith(prefix):
            return trait[len(prefix):]
    return trait

def create_heatmap(df, output_img, metric="effect", cap_percentile=90):
    traits = extract_traits(df)
    if not traits:
        print("Error: No trait columns found in table.")
        return

    # 1. Format clean composite X-axis labels (rsID + wrapped gene names)
    x_labels = []
    for _, row in df.iterrows():
        rsid = row["SNP"]
        genes = set()
        ens = str(row.get("Ensembl_Genes", ""))
        gtex = str(row.get("GTEx_Brain_eQTL_Genes", ""))
        
        if ens and ens not in ["None", "nan"]:
            genes.update([g.strip() for g in ens.split(",") if g.strip()])
        if gtex and gtex not in ["None", "nan"]:
            genes.update([g.strip() for g in gtex.split(",") if g.strip()])
            
        gene_list = sorted(genes)[:3]
        if gene_list:
            gene_str = "\n".join(gene_list)
            x_labels.append(f"{rsid}\n({gene_str})")
        else:
            x_labels.append(f"{rsid}\n(No gene)")

    # Clean Y-axis labels
    y_labels = [clean_trait_name(t) for t in traits]

    # 2. Extract values and significance markers
    matrix_values = np.zeros((len(traits), len(df)))
    annot_matrix = np.empty((len(traits), len(df)), dtype=object)

    for col_idx, (_, row) in enumerate(df.iterrows()):
        for row_idx, trait in enumerate(traits):
            p_col, fdr_col, eff_col, eff_type = find_trait_columns(df, trait)

            # Extract Effect Size
            eff_val = np.nan
            if eff_col and pd.notna(row[eff_col]):
                try:
                    raw_val = float(row[eff_col])
                    if eff_type == "OR":
                        eff_val = np.log(raw_val) if raw_val > 0 else np.nan
                    else:
                        eff_val = raw_val
                except (ValueError, TypeError):
                    eff_val = np.nan

            # Extract P-value and FDR
            try:
                fdr_val = float(row[fdr_col]) if fdr_col and pd.notna(row[fdr_col]) else 1.0
            except (ValueError, TypeError):
                fdr_val = 1.0

            try:
                p_val = float(row[p_col]) if p_col and pd.notna(row[p_col]) else 1.0
            except (ValueError, TypeError):
                p_val = 1.0

            # Determine plotted matrix value
            if metric == "zscore":
                if pd.notna(p_val) and p_val > 0 and pd.notna(eff_val):
                    safe_p = max(p_val, 1e-300)
                    z = abs(norm.ppf(safe_p / 2.0))
                    matrix_values[row_idx, col_idx] = z if eff_val >= 0 else -z
                else:
                    matrix_values[row_idx, col_idx] = 0.0
            else:
                matrix_values[row_idx, col_idx] = eff_val

            # Assign significance symbols
            if fdr_val < 0.05:
                annot_matrix[row_idx, col_idx] = "★"
            elif p_val < 0.05:
                annot_matrix[row_idx, col_idx] = "•"
            else:
                annot_matrix[row_idx, col_idx] = ""

    # Diagnostic output
    valid_vals = matrix_values[~np.isnan(matrix_values)]
    print("\n--- Diagnostic Value Check ---")
    print(f"Plotting Metric: {metric.upper()}")
    print(f"Value range: {valid_vals.min():.4f} to {valid_vals.max():.4f}")

    # 3. Dynamic color scaling with robust percentile clipping
    if metric == "zscore":
        vlimit = max(4.0, float(np.nanpercentile(np.abs(valid_vals), 95)))
        cbar_title = "Signed Z-Score: sign(β) × |Z|"
    else:
        vlimit = float(np.nanpercentile(np.abs(valid_vals), cap_percentile))
        vlimit = max(vlimit, 0.08)
        cbar_title = f"Harmonized Effect Size: Beta / ln(OR)\n(Capped at ±{vlimit:.2f})"

    print(f"Colorbar dynamic range set to: [-{vlimit:.3f}, +{vlimit:.3f}]\n")

    # 4. Dimensioning and Plotting
    fig_width = max(9, len(df) * 2.2 + 3)
    fig_height = max(7, len(traits) * 0.45 + 2.5)
    
    plt.figure(figsize=(fig_width, fig_height))

    ax = sns.heatmap(
        matrix_values,
        xticklabels=x_labels,
        yticklabels=y_labels,
        annot=annot_matrix,
        fmt="",
        cmap="vlag",
        center=0,
        vmin=-vlimit,
        vmax=vlimit,
        linewidths=0.75,
        linecolor="#eaeaea",
        cbar_kws={
            "label": cbar_title,
            "shrink": 0.75
        },
        annot_kws={"fontsize": 13, "color": "black", "weight": "bold"}
    )

    ax.set_title("Cross-Trait Pleiotropy & Brain Gene Architecture", fontsize=15, pad=18, weight="bold")
    ax.set_xlabel("Shared Methylation Variants & Implicated Genes (Ensembl/GTEx)", fontsize=11, labelpad=12, weight="bold")
    ax.set_ylabel("GWAS Traits", fontsize=11, labelpad=12, weight="bold")

    plt.xticks(rotation=0, ha="center", fontsize=9)
    plt.yticks(fontsize=9.5)

    plt.figtext(
        0.5, -0.05,
        "Significance:  ★ FDR q < 0.05  |  • Nominal P < 0.05  |  Empty = Non-significant\nEffect Direction: Red (+) = Risk-increasing / Higher trait score  |  Blue (-) = Protective / Lower trait score",
        wrap=True, horizontalalignment='center', fontsize=9.5, style='italic', color="#333333"
    )

    plt.tight_layout()
    
    # Save to destination directory (creates directory automatically if missing)
    os.makedirs(os.path.dirname(output_img) if os.path.dirname(output_img) else ".", exist_ok=True)
    plt.savefig(output_img, dpi=300, bbox_inches="tight")
    pdf_out = os.path.splitext(output_img)[0] + ".pdf"
    plt.savefig(pdf_out, bbox_inches="tight")
    plt.close()

    print(f"Heatmap successfully generated:")
    print(f"  • PNG: {output_img}")
    print(f"  • PDF: {pdf_out}\n")
    return output_img, pdf_out

def interactive_mode():
    files = sorted(glob.glob(os.path.join(COMPARISONS_DIR, "*_annotated.tsv")))
    if not files:
        print(f"No '*_annotated.tsv' files found in '{COMPARISONS_DIR}/'.")
        sys.exit(1)

    print("==================================================")
    print("      Cross-Trait Pleiotropy Heatmap Generator    ")
    print("==================================================")
    print("\nAvailable annotated comparison files:")
    for idx, f in enumerate(files, 1):
        print(f"  [{idx:2d}] {os.path.basename(f)}")

    choice = input(f"\nSelect file to plot (1-{len(files)}): ").strip()
    try:
        selected_file = files[int(choice) - 1]
    except (ValueError, IndexError):
        print("Invalid selection.")
        sys.exit(1)

    print("\nChoose Color Scale Metric:")
    print("  [1] Harmonized Effect Size with Robust Capping (Beta / ln(OR)) - Recommended")
    print("  [2] Signed Z-Scores (Standardized Statistical Evidence)")
    metric_choice = input("Select metric [1 or 2, default: 1]: ").strip()
    metric = "zscore" if metric_choice == "2" else "effect"

    # Route default output to results/figures/
    base_filename = os.path.basename(selected_file).replace(".tsv", f"_{metric}_heatmap.png")
    default_out = os.path.join(FIGURES_DIR, base_filename)
    
    custom_out = input(f"Output image path [default: {default_out}]: ").strip()
    out_img = custom_out if custom_out else default_out

    df = load_annotated_data(selected_file)
    create_heatmap(df, out_img, metric=metric)

def main():
    parser = argparse.ArgumentParser(description="Generate a cross-trait pleiotropy heatmap from annotated TSV.")
    parser.add_argument("--input", "-i", help="Path to input annotated TSV")
    parser.add_argument("--output", "-o", help="Path for output plot image (.png)")
    parser.add_argument("--metric", choices=["effect", "zscore"], default="effect", help="Plotting metric")
    parser.add_argument("--cap-percentile", type=int, default=90, help="Percentile to cap effect colorbar (default: 90)")
    
    args = parser.parse_args()

    if args.input:
        base_filename = os.path.basename(args.input).replace(".tsv", f"_{args.metric}_heatmap.png")
        default_out = os.path.join(FIGURES_DIR, base_filename)
        out_img = args.output or default_out
        
        df = load_annotated_data(args.input)
        create_heatmap(df, out_img, metric=args.metric, cap_percentile=args.cap_percentile)
    else:
        interactive_mode()

if __name__ == "__main__":
    main()