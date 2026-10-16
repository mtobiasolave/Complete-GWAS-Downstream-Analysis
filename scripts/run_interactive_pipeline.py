#!/usr/bin/env python3
"""
End-to-End Interactive Pipeline Runner for Epigenetic Hub Variant Analysis.
Chains:
  1. Multi-trait overlap comparison & candidate-space FDR correction (Storey q / BH)
  2. Consequence & brain eQTL gene annotation (Ensembl VEP + GTEx v10) with persistent caching
  3. Publication-grade cross-trait pleiotropy heatmaps (PDF & PNG)
"""

import os
import sys
import argparse

# Ensure scripts directory is in sys.path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from compare_overlap_snps import (
    discover_available_files,
    compare_and_filter,
    COMPARISONS_DIR,
    OVERLAP_DIR
)
from annotate_snps_genes import (
    annotate_file,
    DEFAULT_CACHE_FILE
)
from plot_pleiotropy_heatmap import (
    load_annotated_data,
    create_heatmap,
    FIGURES_DIR
)

def run_pipeline_interactive():
    records = discover_available_files()
    if not records:
        print(f"\n[!] Error: No '*_overlap_fullheaders.txt' files found in '{OVERLAP_DIR}/'.")
        print("    Please run the Tier 1 methylation overlap Snakemake step first:")
        print("    snakemake --profile profiles/slurm methylation\n")
        sys.exit(1)

    print("\n" + "=" * 60)
    print("   Epigenetic Pleiotropy Pipeline: Guided End-to-End Runner")
    print("=" * 60)

    # -------------------------------------------------------------
    # Step 1: Selection & Overlap Filtering
    # -------------------------------------------------------------
    print("\n>>> STEP 1/3: Candidate Epigenetic Overlap & FDR Filtering")

    available_meth = sorted(list({r["meth"] for r in records}))
    print("\nAvailable Methylation Types:")
    for idx, m in enumerate(available_meth, 1):
        print(f"  [{idx}] {m}")
    print(f"  [{len(available_meth) + 1}] All Types combined (Recommended)")

    meth_choice = input(f"\nSelect Methylation Type [1-{len(available_meth) + 1}, default: {len(available_meth) + 1}]: ").strip()
    try:
        meth_idx = int(meth_choice) - 1
        chosen_meth = [available_meth[meth_idx]] if meth_idx < len(available_meth) else available_meth
    except (ValueError, IndexError):
        chosen_meth = available_meth

    filtered_records = [r for r in records if r["meth"] in chosen_meth]
    available_traits = sorted(list({r["trait"] for r in filtered_records}))

    print(f"\nAvailable GWAS Traits ({len(available_traits)} found):")
    for idx, t in enumerate(available_traits, 1):
        print(f"  [{idx:2d}] {t}")

    trait_choice = input("\nEnter trait numbers (e.g. '1, 2, 4') or 'all' [default: all]: ").strip().lower()
    if trait_choice in ["all", ""]:
        selected_traits = available_traits
    else:
        try:
            indices = [int(x.strip()) - 1 for x in trait_choice.split(",") if x.strip()]
            selected_traits = [available_traits[i] for i in indices if 0 <= i < len(available_traits)]
        except ValueError:
            print("Invalid trait selection. Aborting.")
            sys.exit(1)

    if len(selected_traits) < 2:
        print("Error: You must select at least 2 traits to compare pleiotropy.")
        sys.exit(1)

    selected_records = [r for r in filtered_records if r["trait"] in selected_traits]

    print("\nMultiple-Testing FDR Correction Method:")
    print("  [1] Storey's q-value (Spline-estimated pi0, matches R qvalue) [Recommended]")
    print("  [2] Benjamini-Hochberg (Standard BH FDR)")
    fdr_sel = input("Select FDR Method [1 or 2, default: 1]: ").strip()
    chosen_method = "bh" if fdr_sel == "2" else "storey"

    total_traits = len(selected_traits)
    print(f"\nSignificance Recurrence Threshold (Comparing {total_traits} traits):")
    default_k = min(3, total_traits)
    k_in = input(f"Enter minimum significant traits k (1 to {total_traits}) [default: {default_k}]: ").strip()
    try:
        min_k = int(k_in)
    except ValueError:
        min_k = default_k

    default_tsv = os.path.join(COMPARISONS_DIR, f"common_{total_traits}_k{min_k}_{chosen_method}_sig_snps.tsv")
    print(f"\nTarget Output TSV: {default_tsv}")

    # Execute Step 1
    cmp_output, final_df = compare_and_filter(selected_records, default_tsv, min_k=min_k, method=chosen_method)
    
    if final_df.empty or len(final_df) == 0:
        print("[!] No SNPs met the significance threshold. Try lowering k or adjusting traits.")
        sys.exit(0)

    # -------------------------------------------------------------
    # Step 2: Functional Annotation (Ensembl VEP + GTEx Brain eQTLs)
    # -------------------------------------------------------------
    print("\n" + "-" * 60)
    print(">>> STEP 2/3: Brain Functional & Positional Annotation")
    print("-" * 60)
    print(f"Annotating {len(final_df)} SNPs via Ensembl VEP & GTEx Portal v10 (13 Brain Regions)...")
    print(f"Persistent Cache: {DEFAULT_CACHE_FILE}")

    annotated_tsv = cmp_output.replace(".tsv", "_annotated.tsv")
    annotate_file(
        input_path=cmp_output,
        output_path=annotated_tsv,
        create_summary=True,
        max_workers=5,
        cache_path=DEFAULT_CACHE_FILE,
        use_cache=True
    )

    # -------------------------------------------------------------
    # Step 3: Publication-Grade Pleiotropy Heatmap
    # -------------------------------------------------------------
    print("\n" + "-" * 60)
    print(">>> STEP 3/3: Publication Pleiotropy Heatmap Generation")
    print("-" * 60)
    print("Select Heatmap Color Metric:")
    print("  [1] Harmonized Effect Size: Beta / ln(OR) with robust capping [Recommended]")
    print("  [2] Signed Standardized Z-Score: sign(Beta) * |Z|")
    metric_choice = input("Select metric [1 or 2, default: 1]: ").strip()
    chosen_metric = "zscore" if metric_choice == "2" else "effect"

    base_plot_name = os.path.basename(annotated_tsv).replace(".tsv", f"_{chosen_metric}_heatmap.png")
    output_png = os.path.join(FIGURES_DIR, base_plot_name)

    annotated_df = load_annotated_data(annotated_tsv)
    out_png, out_pdf = create_heatmap(annotated_df, output_png, metric=chosen_metric)

    print("\n" + "=" * 60)
    print("             ALL 3 PIPELINE STEPS COMPLETED!            ")
    print("=" * 60)
    print(f"1. Comparison Candidate Table : {cmp_output}")
    print(f"2. Annotated Data Table       : {annotated_tsv}")
    print(f"3. Summary Effector Table     : {annotated_tsv.replace('_annotated.tsv', '_summary.tsv')}")
    print(f"4. Publication Vector PDF     : {out_pdf}")
    print(f"5. Publication 300 DPI PNG    : {out_png}")
    print("\nYou can now open and explore these results interactively:")
    print("   streamlit run dashboard/app.py\n")

def run_pipeline_headless(args):
    records = discover_available_files()
    if not records:
        print(f"Error: No '*_overlap_fullheaders.txt' files found in '{OVERLAP_DIR}/'.")
        sys.exit(1)

    meth_filter = ["DNA_ASM", "RNA_m6A", "RNAm6A"] if args.meth == "both" else ([args.meth] if args.meth != "RNAm6A" else ["RNA_m6A", "RNAm6A"])
    filtered = [r for r in records if r["meth"] in meth_filter]

    if args.traits == ["all"] or not args.traits:
        selected_records = filtered
    else:
        selected_records = [r for r in filtered if r["trait"] in args.traits]

    total_traits = len({r["trait"] for r in selected_records})
    default_tsv = os.path.join(COMPARISONS_DIR, f"common_{total_traits}_k{args.min_k}_{args.fdr_method}_sig_snps.tsv")
    cmp_out = args.output or default_tsv

    print(f"--- Running Step 1: Comparison ({total_traits} traits, k>={args.min_k}, {args.fdr_method}) ---")
    cmp_output, final_df = compare_and_filter(selected_records, cmp_out, min_k=args.min_k, method=args.fdr_method)

    if final_df.empty:
        print("No significant variants found.")
        return

    print("\n--- Running Step 2: Annotation (Ensembl + GTEx Brain) ---")
    annotated_tsv = cmp_output.replace(".tsv", "_annotated.tsv")
    annotate_file(
        input_path=cmp_output,
        output_path=annotated_tsv,
        create_summary=True,
        max_workers=args.threads,
        cache_path=args.cache_file,
        use_cache=not args.no_cache
    )

    print(f"\n--- Running Step 3: Heatmap Plotting ({args.metric}) ---")
    base_plot_name = os.path.basename(annotated_tsv).replace(".tsv", f"_{args.metric}_heatmap.png")
    output_png = os.path.join(FIGURES_DIR, base_plot_name)
    annotated_df = load_annotated_data(annotated_tsv)
    create_heatmap(annotated_df, output_png, metric=args.metric)
    print("\nPipeline execution complete.")

def main():
    parser = argparse.ArgumentParser(description="End-to-End Epigenetic Pleiotropy Pipeline Runner.")
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive guided runner")
    parser.add_argument("--traits", nargs="+", help="Traits to compare, or 'all'")
    parser.add_argument("--meth", choices=["DNA_ASM", "RNA_m6A", "RNAm6A", "both"], default="both", help="Methylation panel")
    parser.add_argument("--fdr-method", choices=["storey", "bh"], default="storey", help="FDR adjustment method")
    parser.add_argument("--min-k", type=int, default=3, help="Minimum significant traits recurrence threshold")
    parser.add_argument("--metric", choices=["effect", "zscore"], default="effect", help="Heatmap color metric")
    parser.add_argument("--threads", type=int, default=5, help="Number of concurrent API threads")
    parser.add_argument("--cache-file", default=DEFAULT_CACHE_FILE, help="Path to API JSON cache")
    parser.add_argument("--no-cache", action="store_true", help="Bypass API cache")
    parser.add_argument("--output", "-o", help="Output TSV path for Step 1")

    args = parser.parse_args()

    if args.interactive or len(sys.argv) == 1:
        run_pipeline_interactive()
    else:
        run_pipeline_headless(args)

if __name__ == "__main__":
    main()

