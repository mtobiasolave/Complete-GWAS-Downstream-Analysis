#!/usr/bin/env python3
import os
import sys
import glob
import argparse
from functools import reduce
import pandas as pd
import numpy as np
try:
    from statsmodels.stats.multitest import multipletests
except ImportError:
    def multipletests(pvals, alpha=0.05, method='fdr_bh'):
        """Lightweight vectorized Benjamini-Hochberg fallback if statsmodels is not installed."""
        p = np.asarray(pvals, dtype=float)
        n = len(p)
        if n == 0:
            return np.array([]), np.array([]), None, None
        order = np.argsort(p)
        ranked_p = p[order]
        q_vals = n * ranked_p / np.arange(1, n + 1)
        q_vals = np.minimum.accumulate(q_vals[::-1])[::-1]
        q_vals = np.clip(q_vals, 0.0, 1.0)
        out_q = np.empty(n)
        out_q[order] = q_vals
        reject = out_q <= alpha
        return reject, out_q, None, None

from scipy.interpolate import UnivariateSpline

OVERLAP_DIR = "results/methylation"
COMPARISONS_DIR = os.path.join(OVERLAP_DIR, "comparisons")

def storey_qvalues(p_values, pi0=None):
    """
    Storey q-value calculation replicating the Bioconductor 'qvalue' package in R.
    Estimates pi0 across a grid of lambdas and fits a cubic smoothing spline at lambda=1.
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([]), 1.0

    if pi0 is None:
        # Tuning grid lambda in [0.05, 0.90] matching R qvalue default
        lambdas = np.arange(0.05, 0.95, 0.05)
        pi0_lambdas = np.array([np.sum(p >= l) / (n * (1.0 - l)) for l in lambdas])
        
        # Fit natural cubic spline to evaluate pi0 at lambda = 1.0
        try:
            spline = UnivariateSpline(lambdas, pi0_lambdas, k=3, s=None)
            pi0 = float(spline(1.0))
        except Exception:
            # Fallback if spline fitting cannot converge on small distributions
            pi0 = float(np.mean(pi0_lambdas))
            
        # Bound pi0 in (0, 1]
        pi0 = min(1.0, max(0.01, pi0))

    # Rank p-values and calculate empirical FDR
    order = np.argsort(p)
    ranked_p = p[order]
    
    q_vals = pi0 * n * ranked_p / np.arange(1, n + 1)
    
    # Enforce monotonicity: cumulative minimum from right to left
    q_vals = np.minimum.accumulate(q_vals[::-1])[::-1]
    q_vals = np.clip(q_vals, 0.0, 1.0)
    
    # Restore original index order
    out_q = np.empty(n)
    out_q[order] = q_vals
    return out_q, pi0

def discover_available_files():
    """Scans results/methylation/ for *_overlap_fullheaders.txt files."""
    pattern = os.path.join(OVERLAP_DIR, "*_overlap_fullheaders.txt")
    files = glob.glob(pattern)
    
    records = []
    for f in files:
        basename = os.path.basename(f).replace("_overlap_fullheaders.txt", "")
        if basename.startswith("DNA_ASM_"):
            meth = "DNA_ASM"
            trait = basename[len("DNA_ASM_"):]
        elif basename.startswith("RNA_m6A_"):
            meth = "RNA_m6A"
            trait = basename[len("RNA_m6A_"):]
        elif basename.startswith("RNAm6A_"):
            meth = "RNA_m6A"
            trait = basename[len("RNAm6A_"):]
        else:
            parts = basename.split("_", 1)
            meth = parts[0]
            trait = parts[1] if len(parts) > 1 else ""
            
        records.append({
            "filepath": f,
            "meth": meth,
            "trait": trait
        })
    return records

def identify_columns(df, filepath):
    """Dynamically detects SNP, P-value, and Effect (Beta/OR) columns."""
    cols = {c.upper(): c for c in df.columns}
    
    snp_candidates = ["SNP", "RSID", "MARKERNAME", "ID", "SNPID"]
    snp_col = next((cols[c] for c in snp_candidates if c in cols), None)
    if not snp_col:
        raise ValueError(f"Could not find SNP column in {filepath}")

    p_candidates = ["P", "PVAL", "P_VALUE", "PVALUE", "P-VALUE", "P.VALUE"]
    p_col = next((cols[c] for c in p_candidates if c in cols), None)
    if not p_col:
        raise ValueError(f"Could not find P-value column in {filepath}")

    beta_candidates = ["BETA", "B", "EFFECT", "EFF", "ESTIMATE"]
    or_candidates = ["OR", "ODDS_RATIO", "ODDS-RATIO", "ODDSRATIO"]
    
    beta_col = next((cols[c] for c in beta_candidates if c in cols), None)
    or_col = next((cols[c] for c in or_candidates if c in cols), None)

    if beta_col:
        eff_col = beta_col
        eff_type = "BETA"
    elif or_col:
        eff_col = or_col
        eff_type = "OR"
    else:
        eff_col = None
        eff_type = "UNKNOWN"

    return snp_col, p_col, eff_col, eff_type

def process_trait_file(filepath, trait_name, method="storey"):
    """Loads a fullheader overlap file, runs FDR correction, and standardizes effect direction."""
    df = pd.read_csv(filepath, sep=r'\s+', engine='python')
    snp_col, p_col, eff_col, eff_type = identify_columns(df, filepath)

    df[p_col] = pd.to_numeric(df[p_col], errors='coerce')
    df = df.dropna(subset=[snp_col, p_col]).copy()

    # FDR correction across candidate methylation SNPs
    if method == "storey":
        qvals, est_pi0 = storey_qvalues(df[p_col])
        print(f"    [{trait_name}] Storey pi0 estimated at: {est_pi0:.4f}")
    else:
        _, qvals, _, _ = multipletests(df[p_col], alpha=0.05, method='fdr_bh')

    df['FDR_q'] = qvals

    if eff_col:
        df[eff_col] = pd.to_numeric(df[eff_col], errors='coerce')
        if eff_type == "OR":
            df['Direction'] = df[eff_col].apply(lambda x: "+" if x > 1.0 else ("-" if x < 1.0 else "0"))
        else:
            df['Direction'] = df[eff_col].apply(lambda x: "+" if x > 0.0 else ("-" if x < 0.0 else "0"))
        effect_summary = df[eff_col]
    else:
        df['Direction'] = "NA"
        effect_summary = np.nan

    res = pd.DataFrame({
        "SNP": df[snp_col],
        f"P_{trait_name}": df[p_col],
        f"FDR_{trait_name}": df['FDR_q'],
        f"Effect_{trait_name} ({eff_type})": effect_summary,
        f"Dir_{trait_name}": df['Direction']
    })
    
    return res.drop_duplicates(subset=["SNP"])

def generate_header(selected_records, min_k, method="storey"):
    """Generates descriptive header indicating the recurrence threshold."""
    meth_types = sorted(list({r["meth"] for r in selected_records}))
    traits = sorted(list({r["trait"] for r in selected_records}))

    meth_str = " & ".join(["m6A RNA" if m in ["RNAm6A", "RNA_m6A"] else ("DNA ASM" if m == "DNA_ASM" else m) for m in meth_types])
    total = len(traits)
    method_str = "Storey q" if method == "storey" else "BH FDR"

    if min_k == total:
        criteria = f"{method_str} <= 0.05 in ALL {total} traits"
    elif min_k == 1:
        criteria = f"{method_str} <= 0.05 in >= 1 trait"
    else:
        criteria = f"{method_str} <= 0.05 in >= {min_k} of {total} traits"

    return f"# Common {meth_str} methylation SNPs for {total} GWAS ({criteria})"

def compare_and_filter(selected_records, output_file, min_k=2, method="storey"):
    """
    Groups candidate records per trait, computes FDR correction, merges traits (outer join),
    filters by recurrence threshold (>= min_k), and assesses effect concordance.
    """
    # 1. Group records by unique trait
    trait_groups = {}
    for r in selected_records:
        trait_groups.setdefault(r["trait"], []).append(r)
        
    unique_traits = sorted(list(trait_groups.keys()))
    dfs = []
    print(f"\n--- Processing and Applying FDR Correction ({method.upper()}) ---")

    for trait_name in unique_traits:
        t_records = trait_groups[trait_name]
        raw_trait_dfs = []
        for r in t_records:
            df_part = pd.read_csv(r["filepath"], sep=r'\s+', engine='python')
            raw_trait_dfs.append(df_part)

        combined_trait_df = pd.concat(raw_trait_dfs, ignore_index=True)
        snp_col, p_col, eff_col, eff_type = identify_columns(combined_trait_df, t_records[0]["filepath"])

        combined_trait_df[p_col] = pd.to_numeric(combined_trait_df[p_col], errors='coerce')
        combined_trait_df = combined_trait_df.dropna(subset=[snp_col, p_col]).drop_duplicates(subset=[snp_col]).copy()

        # FDR correction across candidate methylation SNPs for this trait
        if method == "storey":
            qvals, est_pi0 = storey_qvalues(combined_trait_df[p_col])
            sig_snps = (qvals <= 0.05).sum()
            print(f"  • {trait_name}: {len(combined_trait_df):,} candidate SNPs | {sig_snps:,} sig at Storey q <= 0.05 (pi0: {est_pi0:.4f})")
        else:
            _, qvals, _, _ = multipletests(combined_trait_df[p_col], alpha=0.05, method='fdr_bh')
            sig_snps = (qvals <= 0.05).sum()
            print(f"  • {trait_name}: {len(combined_trait_df):,} candidate SNPs | {sig_snps:,} sig at BH FDR <= 0.05")

        combined_trait_df['FDR_q'] = qvals

        if eff_col:
            combined_trait_df[eff_col] = pd.to_numeric(combined_trait_df[eff_col], errors='coerce')
            if eff_type == "OR":
                direction = combined_trait_df[eff_col].apply(lambda x: "+" if x > 1.0 else ("-" if x < 1.0 else "0"))
            else:
                direction = combined_trait_df[eff_col].apply(lambda x: "+" if x > 0.0 else ("-" if x < 0.0 else "0"))
            effect_summary = combined_trait_df[eff_col]
        else:
            direction = "NA"
            effect_summary = np.nan

        res = pd.DataFrame({
            "SNP": combined_trait_df[snp_col],
            f"P_{trait_name}": combined_trait_df[p_col],
            f"FDR_{trait_name}": combined_trait_df['FDR_q'],
            f"Effect_{trait_name} ({eff_type})": effect_summary,
            f"Dir_{trait_name}": direction
        })
        dfs.append(res.drop_duplicates(subset=["SNP"]))

    if not dfs:
        print("\n[!] No valid datasets to compare.")
        return output_file, pd.DataFrame()

    # Outer join to retain candidate variants across phenotypes
    merged = reduce(lambda left, right: pd.merge(left, right, on="SNP", how="outer"), dfs)
    print(f"\nTotal unique candidate SNPs evaluated across phenotypes: {len(merged):,}")

    fdr_cols = [f"FDR_{t}" for t in unique_traits if f"FDR_{t}" in merged.columns]
    dir_cols = [f"Dir_{t}" for t in unique_traits if f"Dir_{t}" in merged.columns]

    if merged.empty or not fdr_cols:
        print(f"\n[!] WARNING: 0 candidate SNPs overlap.")
        final_df = pd.DataFrame(columns=["SNP", "Num_Sig_Traits", "Significant_In", "Concordance"])
    else:
        sig_matrix = (merged[fdr_cols] <= 0.05)
        merged["Num_Sig_Traits"] = sig_matrix.sum(axis=1)

        final_df = merged[merged["Num_Sig_Traits"] >= min_k].copy()

        if final_df.empty:
            print(f"\n[!] Notice: 0 SNPs reached FDR <= 0.05 in >= {min_k} phenotypes.")
            final_df["Significant_In"] = pd.Series(dtype=object)
            final_df["Concordance"] = pd.Series(dtype=object)
        else:
            def label_significance(row):
                sig_list = [t for t in unique_traits if pd.notna(row.get(f"FDR_{t}")) and row.get(f"FDR_{t}") <= 0.05]
                if len(sig_list) == len(unique_traits):
                    return "All traits"
                return ", ".join(sig_list)

            final_df["Significant_In"] = final_df.apply(label_significance, axis=1)

            def check_concordance(row):
                dirs = [str(row[d]).strip() for d in dir_cols if pd.notna(row.get(d)) and str(row[d]).strip() in ["+", "-"]]
                if not dirs:
                    return "NA"
                if len(set(dirs)) == 1:
                    return f"Concordant ({dirs[0]})"
                return "Discordant"

            final_df["Concordance"] = final_df.apply(check_concordance, axis=1)

            final_df["min_FDR"] = final_df[fdr_cols].min(axis=1)
            final_df = final_df.sort_values(
                by=["Num_Sig_Traits", "min_FDR"],
                ascending=[False, True]
            ).drop(columns=["min_FDR"])

    lead_cols = ["SNP", "Num_Sig_Traits", "Significant_In", "Concordance"]
    other_cols = [c for c in final_df.columns if c not in lead_cols]
    final_df = final_df[lead_cols + other_cols]

    header_line = generate_header(selected_records, min_k, method=method)
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else ".", exist_ok=True)
    
    with open(output_file, "w") as f:
        f.write(f"{header_line}\n")
    final_df.to_csv(output_file, sep="\t", index=False, mode="a")

    print("\n--- Summary ---")
    print(f"Header: {header_line}")
    print(f"SNPs saved (significant in >= {min_k} traits): {len(final_df):,}")
    print(f"Saved to: {output_file}\n")
    return output_file, final_df

def interactive_mode(records):
    """Interactive CLI menu."""
    if not records:
        print(f"No '*_overlap_fullheaders.txt' files found in '{OVERLAP_DIR}/'.")
        sys.exit(1)

    print("==================================================")
    print("      Methylation Overlap FDR & Pleiotropy Tool   ")
    print("==================================================")
    
    # 1. Methylation Type
    available_meth = sorted(list({r["meth"] for r in records}))
    print("\nAvailable Methylation Types:")
    for idx, m in enumerate(available_meth, 1):
        print(f"  [{idx}] {m}")
    print(f"  [{len(available_meth)+1}] All Types combined")

    meth_choice = input(f"\nSelect Methylation Type (1-{len(available_meth)+1}): ").strip()
    try:
        meth_idx = int(meth_choice) - 1
        chosen_meth = [available_meth[meth_idx]] if meth_idx < len(available_meth) else available_meth
    except (ValueError, IndexError):
        chosen_meth = available_meth

    filtered_records = [r for r in records if r["meth"] in chosen_meth]
    available_traits = sorted(list({r["trait"] for r in filtered_records}))

    # 2. Traits
    print("\nAvailable GWAS Traits:")
    for idx, t in enumerate(available_traits, 1):
        print(f"  [{idx:2d}] {t}")

    trait_choice = input("\nEnter trait numbers separated by commas (e.g., '1, 3') or 'all': ").strip().lower()
    if trait_choice == "all" or not trait_choice:
        selected_traits = available_traits
    else:
        try:
            indices = [int(x.strip()) - 1 for x in trait_choice.split(",") if x.strip()]
            selected_traits = [available_traits[i] for i in indices if 0 <= i < len(available_traits)]
        except ValueError:
            print("Invalid input. Aborting.")
            sys.exit(1)

    if len(selected_traits) < 2:
        print("Error: You must select at least 2 traits to compare.")
        sys.exit(1)

    selected_records = [r for r in filtered_records if r["trait"] in selected_traits]

    # 3. Significance Method Selection
    print("\nFDR Method:")
    print("  [1] Storey q-value (Matches R Bioconductor 'qvalue') [Default]")
    print("  [2] Benjamini-Hochberg (Standard fdr_bh)")
    fdr_sel = input("Select method [1 or 2; default: 1]: ").strip()
    chosen_method = "bh" if fdr_sel == "2" else "storey"

    # 4. Filter Criterion Prompt
    total_traits = len(selected_traits)
    print(f"\nSignificance Filtering Strategy (Comparing {total_traits} traits):")
    print("  [1] Inclusive (Significant in >= 1 trait)")
    print("  [2] Recurrent Hubs (e.g. Significant in >= 3 traits) - Recommended")
    print(f"  [3] Strict (Significant in ALL {total_traits} traits)")
    
    choice = input("Select strategy [1, 2, or 3; default: 2]: ").strip()
    if choice == "1":
        min_k = 1
    elif choice == "3":
        min_k = total_traits
    else:
        default_k = min(3, total_traits)
        k_in = input(f"Enter minimum number of significant traits k (2 to {total_traits}) [default: {default_k}]: ").strip()
        try:
            min_k = int(k_in)
        except ValueError:
            min_k = default_k

    # 5. Output path
    default_out = os.path.join(COMPARISONS_DIR, f"common_{total_traits}_k{min_k}_{chosen_method}_sig_snps.tsv")
    custom_out = input(f"\nOutput file path [default: {default_out}]: ").strip()
    out_file = custom_out if custom_out else default_out

    compare_and_filter(selected_records, out_file, min_k=min_k, method=chosen_method)

def main():
    parser = argparse.ArgumentParser(description="Find shared significant methylation SNPs across traits.")
    parser.add_argument("--traits", nargs="+", help="Specific trait names to compare, or 'all'")
    parser.add_argument("--meth", choices=["DNA_ASM", "RNAm6A", "RNA_m6A", "both"], default="both", help="Methylation type")
    parser.add_argument("--fdr-method", choices=["storey", "bh"], default="storey",
                        help="FDR calculation method: 'storey' (Storey q-value, default) or 'bh' (Benjamini-Hochberg)")
    parser.add_argument("--min-k", type=int, default=2, help="Minimum number of traits where FDR <= 0.05")
    parser.add_argument("--output", "-o", help="Path for output file (.tsv)")
    parser.add_argument("-i", "--interactive", action="store_true", help="Launch interactive menu")

    args = parser.parse_args()
    records = discover_available_files()

    if args.interactive or len(sys.argv) == 1:
        interactive_mode(records)
    else:
        if args.meth == "both":
            meth_filter = ["DNA_ASM", "RNA_m6A", "RNAm6A"]
        elif args.meth in ["RNAm6A", "RNA_m6A"]:
            meth_filter = ["RNA_m6A", "RNAm6A"]
        else:
            meth_filter = [args.meth]
        filtered = [r for r in records if r["meth"] in meth_filter]

        if args.traits == ["all"]:
            selected_records = filtered
        else:
            selected_records = [r for r in filtered if r["trait"] in args.traits]

        num_unique_traits = len({r["trait"] for r in selected_records})
        default_out = os.path.join(COMPARISONS_DIR, f"common_{num_unique_traits}_k{args.min_k}_{args.fdr_method}_sig_snps.tsv")
        out_file = args.output or default_out
        compare_and_filter(selected_records, out_file, min_k=args.min_k, method=args.fdr_method)

if __name__ == "__main__":
    main()