#!/usr/bin/env python3
"""
End-to-End Interactive & CLI Pipeline Runner for Epigenetic Hub Variant Analysis.
Supports selecting RNA m6A, DNA ASM, or Both modalities across traits:
  - RNA m6A: Multi-trait overlap comparison + FDR correction (Storey q / BH) ->
             Ensembl VEP & GTEx v10 Brain eQTL annotation ->
             Publication-grade cross-trait pleiotropy heatmaps (PDF & PNG).
  - DNA ASM: Per-trait CpG TSS (< 5 kb) & GTEx v10 Brain eQTL comparison ->
             Deduplicated per-SNP FDR correction (Storey q / BH) & cross-trait recurrence summary.
"""

import argparse
import os
import sys

# Ensure scripts directory is in sys.path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from annotate_snps_genes import DEFAULT_CACHE_FILE, annotate_file
from compare_overlap_snps import (
    COMPARISONS_DIR,
    OVERLAP_DIR,
    compare_and_filter,
    discover_available_files,
)
from interactive_cross_trait_asm import (
    DEFAULT_OUTPUT_FILE as DEFAULT_ASM_OUTPUT,
    discover_trait_comparison_files,
    run_cross_trait_asm_analysis,
)
from plot_pleiotropy_heatmap import FIGURES_DIR, create_heatmap, load_annotated_data


def execute_rna_m6a_branch(
    selected_records: list,
    min_k: int,
    chosen_method: str,
    chosen_metric: str = "effect",
    output_tsv: str = None,
    threads: int = 8,
    cache_file: str = DEFAULT_CACHE_FILE,
    use_cache: bool = True,
):
    total_traits = len({r["trait"] for r in selected_records})
    default_tsv = os.path.join(
        COMPARISONS_DIR,
        f"rna_m6a_{total_traits}_k{min_k}_{chosen_method}_sig_snps.tsv",
    )
    target_tsv = output_tsv or default_tsv

    print("\n" + "-" * 60)
    print(f">>> [RNA m6A] Step 1/3: Candidate Overlap & {chosen_method.upper()} FDR Filtering")
    print("-" * 60)
    cmp_output, final_df = compare_and_filter(
        selected_records, target_tsv, min_k=min_k, method=chosen_method
    )

    if final_df.empty or len(final_df) == 0:
        print("[!] No RNA m6A SNPs met the significance threshold. Try lowering k.")
        return None, None, None, None

    print("\n" + "-" * 60)
    print(">>> [RNA m6A] Step 2/3: Brain Functional & Positional Annotation (Ensembl VEP + GTEx v10)")
    print("-" * 60)
    annotated_tsv = cmp_output.replace(".tsv", "_annotated.tsv")
    annotate_file(
        input_path=cmp_output,
        output_path=annotated_tsv,
        create_summary=True,
        max_workers=threads,
        cache_path=cache_file,
        use_cache=use_cache,
    )

    print("\n" + "-" * 60)
    print(f">>> [RNA m6A] Step 3/3: Publication Pleiotropy Heatmap ({chosen_metric})")
    print("-" * 60)
    base_plot_name = os.path.basename(annotated_tsv).replace(
        ".tsv", f"_{chosen_metric}_heatmap.png"
    )
    output_png = os.path.join(FIGURES_DIR, base_plot_name)
    annotated_df = load_annotated_data(annotated_tsv)
    out_png, out_pdf = create_heatmap(annotated_df, output_png, metric=chosen_metric)
    return cmp_output, annotated_tsv, out_pdf, out_png


def execute_dna_asm_branch(
    selected_traits: list,
    min_k: int,
    chosen_method: str,
    threshold: float = 0.05,
    output_tsv: str = DEFAULT_ASM_OUTPUT,
):
    print("\n" + "-" * 60)
    print(
        f">>> [DNA ASM] Cross-Trait {chosen_method.upper()} FDR & CpG TSS (< 5 kb) == GTEx Brain Concordance"
    )
    print("-" * 60)
    wide_df, _, _, _ = run_cross_trait_asm_analysis(
        input_dir=OVERLAP_DIR,
        method=chosen_method,
        threshold=threshold,
        min_traits=min_k,
        output_path=output_tsv,
        selected_traits=selected_traits,
        save_output=True,
    )
    return output_tsv, wide_df


def run_pipeline_interactive():
    records = discover_available_files()
    asm_files = discover_trait_comparison_files(OVERLAP_DIR)

    if not records and not asm_files:
        print(f"\n[!] Error: No methylation overlap files found in '{OVERLAP_DIR}/'.")
        print("    Please run Snakemake first:")
        print("    snakemake --profile profiles/slurm all\n")
        sys.exit(1)

    print("\n" + "=" * 65)
    print("   Epigenetic Pleiotropy Pipeline: Unified Interactive Runner")
    print("=" * 65)

    print("\nSelect Methylation Modality to Analyze:")
    print("  [1] RNA_m6A (Ensembl VEP + GTEx v10 Brain eQTLs + Pleiotropy Heatmap)")
    print("  [2] DNA_ASM (CpG TSS < 5 kb + GTEx v10 Brain eQTL Concordance)")
    print("  [3] Both Modalities (Run both RNA_m6A and DNA_ASM across traits) [Recommended]")

    mod_choice = input("\nSelect Modality [1, 2, or 3, default: 3]: ").strip()
    if mod_choice == "1":
        run_rna, run_dna = True, False
    elif mod_choice == "2":
        run_rna, run_dna = False, True
    else:
        run_rna, run_dna = True, True

    rna_records = [r for r in records if r["meth"] in ("RNA_m6A", "RNAm6A")]
    available_traits = sorted(
        list({r["trait"] for r in rna_records} | set(asm_files.keys()))
    )

    print(f"\nAvailable GWAS Traits ({len(available_traits)} found):")
    for idx, t in enumerate(available_traits, 1):
        print(f"  [{idx:2d}] {t}")

    trait_choice = (
        input("\nEnter trait numbers (e.g. '1, 2, 3') or 'all' [default: all]: ")
        .strip()
        .lower()
    )
    if trait_choice in ["all", ""]:
        selected_traits = available_traits
    else:
        try:
            indices = [int(x.strip()) - 1 for x in trait_choice.split(",") if x.strip()]
            selected_traits = [
                available_traits[i] for i in indices if 0 <= i < len(available_traits)
            ]
        except ValueError:
            print("Invalid trait selection. Aborting.")
            sys.exit(1)

    if len(selected_traits) < 1:
        print("Error: You must select at least 1 trait.")
        sys.exit(1)

    print("\nMultiple-Testing FDR Correction Method:")
    print("  [1] Storey's q-value (Spline-estimated pi0, matches R qvalue) [Recommended]")
    print("  [2] Benjamini-Hochberg (Standard BH FDR)")
    fdr_sel = input("Select FDR Method [1 or 2, default: 1]: ").strip()
    chosen_method = "bh" if fdr_sel == "2" else "storey"

    total_traits = len(selected_traits)
    print(f"\nSignificance Recurrence Threshold (Comparing {total_traits} traits):")
    default_k = min(2, total_traits)
    k_in = input(
        f"Enter minimum significant traits k (1 to {total_traits}, or 'all') [default: {default_k}]: "
    ).strip()
    if k_in.lower() == "all":
        min_k = total_traits
    else:
        try:
            min_k = int(k_in)
        except ValueError:
            min_k = default_k

    chosen_metric = "effect"
    if run_rna:
        print("\nSelect RNA m6A Heatmap Color Metric:")
        print("  [1] Harmonized Effect Size: Beta / ln(OR) with robust capping [Recommended]")
        print("  [2] Signed Standardized Z-Score: sign(Beta) * |Z|")
        metric_choice = input("Select metric [1 or 2, default: 1]: ").strip()
        chosen_metric = "zscore" if metric_choice == "2" else "effect"

    if run_rna:
        selected_rna_records = [
            r for r in rna_records if r["trait"] in selected_traits
        ]
        if selected_rna_records:
            execute_rna_m6a_branch(
                selected_rna_records,
                min_k=min_k,
                chosen_method=chosen_method,
                chosen_metric=chosen_metric,
            )
        else:
            print("[!] No RNA_m6A overlap files found for the selected traits.")

    if run_dna:
        if asm_files:
            execute_dna_asm_branch(
                selected_traits=selected_traits,
                min_k=min_k,
                chosen_method=chosen_method,
                output_tsv=DEFAULT_ASM_OUTPUT,
            )
        else:
            print("[!] No DNA ASM comparison files found in results/methylation/.")

    print("\n" + "=" * 65)
    print("             INTERACTIVE PIPELINE COMPLETED!                 ")
    print("=" * 65)
    print("Explore both RNA m6A and DNA ASM results interactively in the dashboard:")
    print("   streamlit run dashboard/app.py --server.port 8502\n")


def run_pipeline_headless(args):
    records = discover_available_files()
    asm_files = discover_trait_comparison_files(OVERLAP_DIR)

    run_rna = args.meth in ("RNA_m6A", "RNAm6A", "both")
    run_dna = args.meth in ("DNA_ASM", "both")

    rna_records = [r for r in records if r["meth"] in ("RNA_m6A", "RNAm6A")]
    all_traits = sorted(list({r["trait"] for r in rna_records} | set(asm_files.keys())))

    if not args.traits or args.traits == ["all"]:
        selected_traits = all_traits
    else:
        selected_traits = [t for t in args.traits if t in all_traits]

    if run_rna:
        selected_rna_records = [
            r for r in rna_records if r["trait"] in selected_traits
        ]
        if selected_rna_records:
            execute_rna_m6a_branch(
                selected_rna_records,
                min_k=args.min_k,
                chosen_method=args.fdr_method,
                chosen_metric=args.metric,
                output_tsv=args.output if not run_dna else None,
                threads=args.threads,
                cache_file=args.cache_file,
                use_cache=not args.no_cache,
            )

    if run_dna:
        asm_out = args.output if (args.output and not run_rna) else DEFAULT_ASM_OUTPUT
        execute_dna_asm_branch(
            selected_traits=selected_traits,
            min_k=args.min_k,
            chosen_method=args.fdr_method,
            output_tsv=asm_out,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Unified Epigenetic Pleiotropy Pipeline Runner (RNA m6A & DNA ASM)."
    )
    parser.add_argument(
        "--interactive", "-i", action="store_true", help="Launch interactive guided runner"
    )
    parser.add_argument("--traits", nargs="+", help="Traits to compare, or 'all'")
    parser.add_argument(
        "--meth",
        choices=["DNA_ASM", "RNA_m6A", "RNAm6A", "both"],
        default="both",
        help="Methylation modality: 'RNA_m6A', 'DNA_ASM', or 'both' (default: both)",
    )
    parser.add_argument(
        "--fdr-method",
        choices=["storey", "bh"],
        default="storey",
        help="FDR adjustment method: 'storey' or 'bh'",
    )
    parser.add_argument(
        "--min-k",
        type=int,
        default=2,
        help="Minimum significant traits recurrence threshold (default: 2)",
    )
    parser.add_argument(
        "--metric",
        choices=["effect", "zscore"],
        default="effect",
        help="RNA m6A heatmap color metric",
    )
    parser.add_argument(
        "--threads", type=int, default=8, help="Number of concurrent API threads"
    )
    parser.add_argument(
        "--cache-file", default=DEFAULT_CACHE_FILE, help="Path to API JSON cache"
    )
    parser.add_argument("--no-cache", action="store_true", help="Bypass API cache")
    parser.add_argument("--output", "-o", help="Output TSV path")

    args = parser.parse_args()

    if args.interactive or len(sys.argv) == 1:
        run_pipeline_interactive()
    else:
        run_pipeline_headless(args)


if __name__ == "__main__":
    main()
