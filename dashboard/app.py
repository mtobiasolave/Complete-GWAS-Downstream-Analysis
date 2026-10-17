import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
import re
import os
import sys
import io
from pathlib import Path

# -----------------------------------------------------------------------------
# Configuration & Page Setup
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Pleiotropy & Functional Epigenomics Explorer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Cross-Trait Epigenetic Hub Explorer")
st.markdown(
    "Visualizing pleiotropic candidate variants ($m^6\\text{A}$ & DNA ASM) "
    "across psychiatric and behavioral phenotypes."
)

# -----------------------------------------------------------------------------
# Path Detection & Dynamic Dataset Discovery
# -----------------------------------------------------------------------------
current_path = Path(__file__).resolve()
PROJECT_ROOT = current_path.parent if (current_path.parent / "results").exists() else current_path.parent.parent
OVERLAP_DIR = str(PROJECT_ROOT / "results/methylation")
COMPARISONS_DIR = PROJECT_ROOT / "results/methylation/comparisons"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

# Add scripts directory to path for importing comparison and annotation tools
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

try:
    from compare_overlap_snps import discover_available_files, compare_and_filter
    from annotate_snps_genes import annotate_file, DEFAULT_CACHE_FILE
    from interactive_cross_trait_asm import (
        discover_trait_comparison_files,
        run_cross_trait_asm_analysis,
    )
except ImportError:
    discover_available_files = None
    compare_and_filter = None
    annotate_file = None
    discover_trait_comparison_files = None
    run_cross_trait_asm_analysis = None
    DEFAULT_CACHE_FILE = str(PROJECT_ROOT / "results/methylation/.api_cache.json")


def discover_available_datasets(directory: Path) -> dict:
    """
    Scans the results/methylation/comparisons directory for RNA m6A summary files
    (*_summary.tsv) and pairs each with its corresponding *_annotated.tsv file.
    """
    if not directory.exists():
        return {}

    summary_files = sorted(directory.glob("*_summary.tsv"))
    dataset_map = {}

    for s_path in summary_files:
        name_lower = s_path.name.lower()
        if "dna_asm" in name_lower or "cross_trait_asm" in name_lower:
            continue

        # Match companion annotated file
        a_path = s_path.parent / s_path.name.replace("_summary.tsv", "_annotated.tsv")
        if a_path.exists():
            # Check header comment if present to exclude pure DNA ASM runs
            try:
                with open(s_path, "r", encoding="utf-8") as fh:
                    first_line = fh.readline()
                if first_line.startswith("#") and "DNA ASM" in first_line and "m6A" not in first_line:
                    continue
            except OSError:
                pass

            stem_name = s_path.name.replace("_sig_snps_summary.tsv", "").replace("_summary.tsv", "")

            # Format clean display label for RNA m6A runs
            if "_storey" in stem_name.lower():
                label = f"RNA m6A — Storey's q-value ({stem_name})"
            elif "_bh" in stem_name.lower():
                label = f"RNA m6A — Benjamini-Hochberg BH ({stem_name})"
            else:
                label = f"RNA m6A ({stem_name})"

            dataset_map[label] = (str(s_path), str(a_path))

    return dataset_map


# -----------------------------------------------------------------------------
# Data Harmonization Engine
# -----------------------------------------------------------------------------
def parse_numeric(series: pd.Series) -> pd.Series:
    """Robust numeric cleaner: handles comma decimals, scientific notation, and dirty strings."""
    if series is None or series.empty:
        return pd.Series(dtype=float)
    cleaned_str = series.astype(str).str.strip().str.replace(",", ".", regex=False)
    extracted = cleaned_str.str.extract(r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")[0]
    return pd.to_numeric(extracted, errors="coerce")


def reshape_wide_to_long(df_wide: pd.DataFrame):
    """
    Dynamically identifies trait-specific columns, strips parenthetical metric tags
    to cleanly match P_, FDR_, and Dir_ columns, and derives standardized Z-scores.
    """
    effect_cols = [c for c in df_wide.columns if re.search(r"^effect_", c, re.IGNORECASE)]
    long_dfs = []
    debug_matching = {}

    for col in effect_cols:
        is_or = bool(re.search(r"\(or\)", col, re.IGNORECASE))

        # 1. Strip 'Effect_' prefix
        stem = re.sub(r"^effect_", "", col, flags=re.IGNORECASE).strip()

        # 2. Extract unit type if present
        unit_match = re.search(r"\((.*?)\)", stem)
        unit_str = unit_match.group(1).strip() if unit_match else ("OR" if is_or else "BETA")

        # 3. Strip parentheses completely to get the base trait key
        raw_trait_key = re.sub(r"\s*\(.*?\).*", "", stem).strip()

        # 4. Clean display label for plots (strip methylation prefixes)
        clean_trait = re.sub(r"^(m6A_|DNA_ASM_|ASM_|RNA_m6A_|RNAm6A_)", "", raw_trait_key, flags=re.IGNORECASE)
        display_label = f"{clean_trait} ({unit_str})"

        raw_values = parse_numeric(df_wide[col])
        if is_or:
            beta_equiv = np.log(raw_values.where(raw_values > 0))
        else:
            beta_equiv = raw_values

        col_pool = list(df_wide.columns)

        # Exact and fuzzy matching for P, FDR, and Dir columns
        p_col = next((c for c in col_pool if c == f"P_{raw_trait_key}"), None)
        if not p_col:
            p_candidates = [c for c in col_pool if c.startswith("P_") and clean_trait.lower() in c.lower()]
            p_col = p_candidates[0] if p_candidates else None

        fdr_col = next((c for c in col_pool if c == f"FDR_{raw_trait_key}"), None)
        if not fdr_col:
            fdr_candidates = [c for c in col_pool if c.startswith("FDR_") and clean_trait.lower() in c.lower()]
            fdr_col = fdr_candidates[0] if fdr_candidates else None

        dir_col = next((c for c in col_pool if c == f"Dir_{raw_trait_key}"), None)
        if not dir_col:
            dir_candidates = [c for c in col_pool if c.startswith("Dir_") and clean_trait.lower() in c.lower()]
            dir_col = dir_candidates[0] if dir_candidates else None

        p_series = parse_numeric(df_wide[p_col]) if p_col else pd.Series(np.nan, index=df_wide.index)
        fdr_series = parse_numeric(df_wide[fdr_col]) if fdr_col else pd.Series(np.nan, index=df_wide.index)
        dir_series = df_wide[dir_col].astype(str).str.strip() if dir_col else pd.Series("+", index=df_wide.index)

        # Derive Signed Z-Score: |Z| = norm.isf(P/2)
        p_safe = p_series.clip(lower=1e-300, upper=0.999999)
        z_abs = stats.norm.isf(p_safe / 2.0)
        z_abs = np.nan_to_num(z_abs, nan=np.nan, posinf=37.0, neginf=0.0)

        # Directional sign assignment
        sign_mult = np.where(dir_series.str.startswith("-"), -1.0, 1.0)
        if dir_col is None:
            sign_mult = np.where(beta_equiv < 0, -1.0, 1.0)

        z_scores = sign_mult * z_abs
        z_scores = np.where(p_series.isna(), np.nan, z_scores)

        # Trait-Normalized Beta: standardized effect across variants per trait
        std_beta = np.nanstd(beta_equiv)
        z_normalized_beta = (beta_equiv / std_beta) if (std_beta > 0) else beta_equiv

        debug_matching[display_label] = {
            "Matched_Key": raw_trait_key,
            "Effect_Col": col,
            "P_Col": p_col,
            "FDR_Col": fdr_col,
            "Dir_Col": dir_col,
            "Valid_P_Count": int(p_series.notna().sum())
        }

        sub_df = pd.DataFrame({
            "SNP": df_wide["SNP"],
            "Trait": display_label,
            "Beta_equiv": beta_equiv,
            "Z_score": z_scores,
            "Z_norm_beta": z_normalized_beta,
            "Raw_Effect": raw_values,
            "Metric_Type": unit_str,
            "P_val": p_series,
            "FDR": fdr_series,
            "Direction": dir_series
        })
        long_dfs.append(sub_df)

    if not long_dfs:
        return pd.DataFrame(columns=["SNP", "Trait", "Beta_equiv", "Z_score", "Z_norm_beta", "FDR", "P_val"]), debug_matching

    return pd.concat(long_dfs, ignore_index=True), debug_matching


@st.cache_data
def load_data(summary_path: str, annotated_path: str):
    df_summary = pd.read_csv(summary_path, sep="\t", comment="#")
    df_wide = pd.read_csv(annotated_path, sep="\t", comment="#")

    # Clean whitespace in column headers
    df_summary.columns = [c.strip() for c in df_summary.columns]
    df_wide.columns = [c.strip() for c in df_wide.columns]

    if "Num_Sig_Traits" in df_summary.columns:
        df_summary["Num_Sig_Traits"] = pd.to_numeric(df_summary["Num_Sig_Traits"], errors="coerce").fillna(0).astype(int)

    df_long, debug_info = reshape_wide_to_long(df_wide)
    return df_summary, df_wide, df_long, debug_info


def generate_publication_figure_bytes(df, metric="effect", cap_percentile=90):
    """
    Renders a publication-grade heatmap matching scripts/plot_pleiotropy_heatmap.py
    and returns high-resolution vector PDF and 300 DPI PNG bytes in memory.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fdr_cols = [c for c in df.columns if c.startswith("FDR_")]
    traits = [c[len("FDR_"):] for c in fdr_cols]
    if not traits:
        return None, None

    def clean_t_name(t):
        for prefix in ["m6A_", "DNA_ASM_", "RNA_m6A_", "RNAm6A_"]:
            if t.startswith(prefix):
                return t[len(prefix):]
        return t

    def find_t_cols(d, t):
        p_col = f"P_{t}" if f"P_{t}" in d.columns else None
        fdr_col = f"FDR_{t}" if f"FDR_{t}" in d.columns else None
        eff_col = None
        eff_type = "BETA"
        candidates = [f"Effect_{t} (BETA)", f"Effect_{t} (OR)", f"Effect_{t}"]
        for c in candidates:
            if c in d.columns:
                eff_col = c
                if "(OR)" in c.upper():
                    eff_type = "OR"
                break
        if not eff_col:
            for c in d.columns:
                if c.startswith(f"Effect_{t}"):
                    eff_col = c
                    if "(OR)" in c.upper():
                        eff_type = "OR"
                    break
        return p_col, fdr_col, eff_col, eff_type

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
            x_labels.append(f"{rsid}\n({chr(10).join(gene_list)})")
        else:
            x_labels.append(f"{rsid}\n(No gene)")

    y_labels = [clean_t_name(t) for t in traits]
    matrix_values = np.zeros((len(traits), len(df)))
    annot_matrix = np.empty((len(traits), len(df)), dtype=object)

    for col_idx, (_, row) in enumerate(df.iterrows()):
        for row_idx, trait in enumerate(traits):
            p_col, fdr_col, eff_col, eff_type = find_t_cols(df, trait)
            eff_val = np.nan
            if eff_col and pd.notna(row[eff_col]):
                try:
                    raw_val = float(row[eff_col])
                    eff_val = np.log(raw_val) if (eff_type == "OR" and raw_val > 0) else raw_val
                except (ValueError, TypeError):
                    eff_val = np.nan

            try:
                fdr_val = float(row[fdr_col]) if fdr_col and pd.notna(row[fdr_col]) else 1.0
            except (ValueError, TypeError):
                fdr_val = 1.0

            try:
                p_val = float(row[p_col]) if p_col and pd.notna(row[p_col]) else 1.0
            except (ValueError, TypeError):
                p_val = 1.0

            if metric == "zscore":
                if pd.notna(p_val) and p_val > 0 and pd.notna(eff_val):
                    safe_p = max(p_val, 1e-300)
                    z = abs(stats.norm.ppf(safe_p / 2.0))
                    matrix_values[row_idx, col_idx] = z if eff_val >= 0 else -z
                else:
                    matrix_values[row_idx, col_idx] = 0.0
            else:
                matrix_values[row_idx, col_idx] = eff_val

            if fdr_val < 0.05:
                annot_matrix[row_idx, col_idx] = "★"
            elif p_val < 0.05:
                annot_matrix[row_idx, col_idx] = "•"
            else:
                annot_matrix[row_idx, col_idx] = ""

    valid_vals = matrix_values[~np.isnan(matrix_values)]
    if len(valid_vals) == 0:
        valid_vals = np.array([0.0])

    if metric == "zscore":
        vlimit = max(4.0, float(np.nanpercentile(np.abs(valid_vals), 95))) if np.any(np.abs(valid_vals) > 0) else 4.0
        cbar_title = "Signed Z-Score: sign(β) × |Z|"
    else:
        vlimit = float(np.nanpercentile(np.abs(valid_vals), cap_percentile)) if np.any(np.abs(valid_vals) > 0) else 0.5
        vlimit = max(vlimit, 0.08)
        cbar_title = f"Harmonized Effect Size: Beta / ln(OR)\n(Capped at ±{vlimit:.2f})"

    fig_width = max(9, len(df) * 2.2 + 3)
    fig_height = max(7, len(traits) * 0.45 + 2.5)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    try:
        import seaborn as sns
        sns.heatmap(
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
            cbar_kws={"label": cbar_title, "shrink": 0.75},
            annot_kws={"fontsize": 13, "color": "black", "weight": "bold"},
            ax=ax
        )
    except ImportError:
        c = ax.imshow(matrix_values, cmap="coolwarm", vmin=-vlimit, vmax=vlimit, aspect="auto")
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels)
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels)
        plt.colorbar(c, ax=ax, label=cbar_title, shrink=0.75)
        for r_i in range(len(traits)):
            for c_i in range(len(df)):
                sym = annot_matrix[r_i, c_i]
                if sym:
                    ax.text(c_i, r_i, sym, ha="center", va="center", color="black", fontsize=13, weight="bold")

    ax.set_title("RNA m6A Cross-Trait Pleiotropy & Brain Gene Architecture", fontsize=15, pad=18, weight="bold")
    ax.set_xlabel("Shared RNA m6A Variants & Implicated Genes (Ensembl VEP / GTEx Brain)", fontsize=11, labelpad=12, weight="bold")
    ax.set_ylabel("GWAS Traits", fontsize=11, labelpad=12, weight="bold")

    plt.xticks(rotation=0, ha="center", fontsize=9)
    plt.yticks(fontsize=9.5)

    plt.figtext(
        0.5, -0.05,
        "Significance:  ★ FDR q < 0.05  |  • Nominal P < 0.05  |  Empty = Non-significant\nEffect Direction: Red (+) = Risk-increasing / Higher trait score  |  Blue (-) = Protective / Lower trait score",
        wrap=True, horizontalalignment='center', fontsize=9.5, style='italic', color="#333333"
    )
    plt.tight_layout()

    pdf_buf = io.BytesIO()
    fig.savefig(pdf_buf, format="pdf", bbox_inches="tight")
    pdf_bytes = pdf_buf.getvalue()
    pdf_buf.close()

    png_buf = io.BytesIO()
    fig.savefig(png_buf, format="png", dpi=300, bbox_inches="tight")
    png_bytes = png_buf.getvalue()
    png_buf.close()

    plt.close(fig)
    return pdf_bytes, png_bytes


# -----------------------------------------------------------------------------
# Main Tabs: RNA m6A Explorer, DNA ASM Functional Convergence, & RNA m6A Generator
# -----------------------------------------------------------------------------
tab_explore, tab_asm_convergence, tab_generate = st.tabs([
    "🧪 RNA m6A Methylation Explorer",
    "🧬 DNA ASM Functional Convergence",
    "⚙️ Run New RNA m6A Comparison & Annotation",
])

available_datasets = discover_available_datasets(COMPARISONS_DIR)

# -----------------------------------------------------------------------------
# TAB 1: RNA m6A METHYLATION EXPLORER (Ensembl VEP ∩ GTEx Brain eQTLs)
# -----------------------------------------------------------------------------
with tab_explore:
    st.subheader("🧪 RNA m6A Cross-Trait Pleiotropy & Functional Convergence")
    st.markdown(
        "Integrates **RNA $N^6$-methyladenosine ($m^6\\text{A}$)** candidate variants across GWAS phenotypes, "
        "evaluating functional convergence between **Ensembl VEP transcript consequences** "
        "and **GTEx v10 human brain cis-eQTL target genes (13 regions)**."
    )

    if not available_datasets:
        st.info(
            "👋 **Welcome!** No pre-computed RNA $m^6\\text{A}$ comparison runs were found in `results/methylation/comparisons/`.\n\n"
            "Please switch to the **⚙️ Run New RNA m6A Comparison & Annotation** tab above to generate your first RNA $m^6\\text{A}$ multi-trait analysis dataset!"
        )
    else:
        # Sidebar: Dataset Selection
        st.sidebar.header("🧪 RNA m6A Dataset & FDR")

        selected_dataset_label = st.sidebar.selectbox(
            "Select RNA m6A Analysis Run:",
            options=list(available_datasets.keys()),
            index=0,
            help="Toggle between Storey's q-value and Benjamini-Hochberg (BH) RNA m6A runs."
        )

        selected_summary_path, selected_annotated_path = available_datasets[selected_dataset_label]
        df_summary, df_wide, df_long, debug_info = load_data(selected_summary_path, selected_annotated_path)

        # Sidebar: Filter & Thresholds
        st.sidebar.markdown("---")
        st.sidebar.header("RNA m6A Filter & Thresholds")

        max_observed_k = int(df_summary["Num_Sig_Traits"].max()) if not df_summary.empty else 3
        default_k = min(2, max(1, max_observed_k))

        min_k = st.sidebar.slider(
            "Minimum Trait Recurrence (k)",
            min_value=1,
            max_value=max(2, max_observed_k),
            value=default_k,
            step=1
        )

        require_gene = st.sidebar.checkbox(
            "Filter to Convergent Genes Only (Ensembl == GTEx Brain)",
            value=False
        )

        st.sidebar.markdown("---")
        st.sidebar.header("RNA m6A Visualization Settings")

        heatmap_metric = st.sidebar.radio(
            "Heatmap Metric",
            options=[
                "Standardized Z-Score (from P-value)",
                "Trait-Normalized Effect (Z-scaled Beta)",
                "Raw Effect [Beta / ln(OR)]"
            ],
            index=0
        )

        # Apply Variant Filters
        mask = df_summary["Num_Sig_Traits"] >= min_k
        if require_gene:
            mask &= (~df_summary["Common_Genes"].isin(["None", "No match", np.nan]))

        filtered_summary = df_summary[mask].copy()
        valid_snps = filtered_summary["SNP"].tolist()

        filtered_long = df_long[df_long["SNP"].isin(valid_snps)].copy()
        filtered_wide = df_wide[df_wide["SNP"].isin(valid_snps)].copy()

        # KPI Overview Cards
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Selected Hub Variants", len(filtered_summary))
        col2.metric("Min Recurrence (k)", f"≥ {min_k}")
        unique_genes = filtered_summary[~filtered_summary["Common_Genes"].isin(["None", "No match", np.nan])]["Common_Genes"].nunique() if not filtered_summary.empty else 0
        col3.metric("Prioritized Effector Genes", unique_genes)
        col4.metric("Analyzed GWAS Phenotypes", df_long["Trait"].nunique())

        st.markdown("---")

        # Main Visualization: Interactive Pleiotropy Heatmap
        st.subheader("Cross-Phenotype Effect Architecture")

        if "from P-value" in heatmap_metric:
            active_val_col = "Z_score"
            colorbar_title = "Signed Z-Score"
            export_metric = "zscore"
        elif "Trait-Normalized" in heatmap_metric:
            active_val_col = "Z_norm_beta"
            colorbar_title = "Normalized Effect (Z)"
            export_metric = "effect"
        else:
            active_val_col = "Beta_equiv"
            colorbar_title = "Effect (Beta / ln(OR))"
            export_metric = "effect"

        if not filtered_long.empty and len(valid_snps) > 0:
            heatmap_df = filtered_long.pivot_table(
                index="Trait",
                columns="SNP",
                values=active_val_col,
                aggfunc="first"
            ).fillna(0)

            heatmap_df = heatmap_df.replace([np.inf, -np.inf], 0)

            gene_map = filtered_summary.set_index("SNP")["Common_Genes"].to_dict()
            col_labels = [f"<b>{snp}</b><br>({gene_map.get(snp, 'None')})" for snp in heatmap_df.columns]

            if active_val_col in ["Z_score", "Z_norm_beta"]:
                abs_vals = np.abs(heatmap_df.values)
                p95 = float(np.percentile(abs_vals[abs_vals > 0], 95)) if (abs_vals > 0).any() else 4.0
                vmax = max(min(p95, 10.0), 2.5)
            else:
                abs_vals = np.abs(heatmap_df.values)
                vmax = float(np.percentile(abs_vals, 90)) if abs_vals.size > 0 else 0.5
                vmax = max(vmax, 0.05)

            fig = px.imshow(
                heatmap_df,
                x=col_labels,
                y=heatmap_df.index,
                color_continuous_scale="RdBu_r",
                range_color=[-vmax, vmax],
                aspect="auto",
                labels=dict(x="Variant (Target Gene)", y="Trait", color=colorbar_title)
            )

            fig.update_layout(
                height=520,
                margin=dict(l=40, r=40, t=20, b=90),
                xaxis_tickangle=-45,
                coloraxis_colorbar=dict(title=colorbar_title)
            )
            st.plotly_chart(fig, width="stretch")

            # -------------------------------------------------------------
            # Feature: Direct Publication Figure Export
            # -------------------------------------------------------------
            with st.expander("🎨 Export Publication-Grade Heatmap (Vector PDF / 300 DPI PNG)", expanded=False):
                st.markdown(
                    "Generate vector and high-resolution publication figures matching "
                    "`plot_pleiotropy_heatmap.py`, with significance asterisks (★ FDR $q < 0.05$, • Nominal $P < 0.05$) "
                    "and composite variant annotations."
                )
                pub_col1, pub_col2 = st.columns([2, 2])
                with pub_col1:
                    pub_metric_sel = st.selectbox(
                        "Colorbar Scale Metric for Figure:",
                        ["Harmonized Effect Size (Beta / ln(OR))", "Signed Standardized Z-Score"],
                        index=0 if export_metric == "effect" else 1
                    )
                    chosen_pub_metric = "zscore" if "Z-Score" in pub_metric_sel else "effect"

                with pub_col2:
                    if st.button("🖼️ Render High-Res Figures", type="secondary"):
                        with st.spinner("Rendering publication-quality figure in memory..."):
                            pdf_bytes, png_bytes = generate_publication_figure_bytes(filtered_wide, metric=chosen_pub_metric)
                            st.session_state["pub_pdf"] = pdf_bytes
                            st.session_state["pub_png"] = png_bytes
                            st.session_state["pub_name"] = Path(selected_annotated_path).stem

                if "pub_pdf" in st.session_state and st.session_state["pub_pdf"] is not None:
                    d_col1, d_col2 = st.columns(2)
                    d_col1.download_button(
                        label="📥 Download Vector PDF",
                        data=st.session_state["pub_pdf"],
                        file_name=f"{st.session_state['pub_name']}_{chosen_pub_metric}.pdf",
                        mime="application/pdf"
                    )
                    d_col2.download_button(
                        label="📥 Download 300 DPI PNG",
                        data=st.session_state["pub_png"],
                        file_name=f"{st.session_state['pub_name']}_{chosen_pub_metric}.png",
                        mime="image/png"
                    )

            if active_val_col == "Z_score":
                st.caption(
                    "ℹ️ **Signed Z-score view active:** Standardized statistical effect sizes derived from association P-values "
                    "($\vert{}Z\vert{} = \\Phi^{-1}(1 - P/2)$). Positive (Red) = Risk; Negative (Blue) = Protective."
                )
            elif active_val_col == "Z_norm_beta":
                st.caption(
                    "ℹ️ **Trait-Normalized view active:** Standardizes effect sizes across variants per phenotype ($\\beta / \\sigma_{\\text{trait}}$), "
                    "resolving cohort sample size ($N$) disparities while maintaining directional signs."
                )
        else:
            st.info("No variants match the selected filter criteria. Try lowering the recurrence threshold (k).")

        # Data Table: Convergent Target Summary
        st.markdown("---")
        st.subheader("Prioritized Master Hubs & Effector Genes")

        col_tbl_left, col_tbl_right = st.columns([3, 1])

        with col_tbl_left:
            st.markdown("Convergent genes identified via **Ensembl VEP positional overlap $\\cap$ GTEx v10 Brain cis-eQTLs**:")

        with col_tbl_right:
            csv_data = filtered_summary.sort_values(by="Num_Sig_Traits", ascending=False).to_csv(sep="\t", index=False)
            st.download_button(
                label="📥 Download Table (TSV)",
                data=csv_data,
                file_name=f"pleiotropic_hubs_k{min_k}.tsv",
                mime="text/tab-separated-values"
            )

        st.dataframe(
            filtered_summary.sort_values(by="Num_Sig_Traits", ascending=False),
            width="stretch"
        )

        # Single-Variant Inspector
        if not filtered_summary.empty and len(valid_snps) > 0:
            st.markdown("---")
            st.subheader("🔍 Single-Variant Inspector")

            selected_snp = st.selectbox("Select Variant to Inspect:", valid_snps)
            snp_summary = filtered_summary[filtered_summary["SNP"] == selected_snp].iloc[0]
            snp_traits = filtered_long[filtered_long["SNP"] == selected_snp].copy()

            col_left, col_right = st.columns([1, 1.4])

            with col_left:
                st.markdown(f"### `{selected_snp}`")
                st.markdown(f"- **Convergent Effector Gene:** `{snp_summary.get('Common_Genes', 'None')}`")
                st.markdown(f"- **Significant Traits ($k$):** `{snp_summary.get('Num_Sig_Traits', 0)}`")
                if "Significant_In" in snp_summary:
                    st.markdown(f"- **Active In:** {snp_summary['Significant_In']}")
                if "Concordance" in snp_summary:
                    st.markdown(f"- **Directional Concordance:** `{snp_summary['Concordance']}`")

                st.markdown("**Per-Trait Metrics Summary:**")
                display_cols = ["Trait", "Metric_Type", "Raw_Effect", "Beta_equiv", "Z_score", "P_val", "FDR"]
                available_display_cols = [c for c in display_cols if c in snp_traits.columns]

                st.dataframe(
                    snp_traits[available_display_cols]
                    .sort_values(by="FDR", ascending=True)
                    .reset_index(drop=True),
                    height=260,
                    width="stretch"
                )

            with col_right:
                if not snp_traits.empty:
                    bar_metric = active_val_col
                    if snp_traits[bar_metric].isna().all():
                        bar_metric = "Beta_equiv"
                        bar_title = "Beta / ln(OR) (Fallback)"
                    else:
                        bar_title = colorbar_title

                    bar_fig = px.bar(
                        snp_traits,
                        x="Trait",
                        y=bar_metric,
                        color=bar_metric,
                        color_continuous_scale="RdBu_r",
                        title=f"Cross-Trait Profile for {selected_snp} ({bar_title})"
                    )
                    bar_fig.add_hline(y=0, line_dash="dash", line_color="gray")
                    bar_fig.update_layout(height=380, xaxis_tickangle=-45, showlegend=False)
                    st.plotly_chart(bar_fig, width="stretch")

        # Diagnostics Sidebar Expander
        with st.sidebar.expander("🛠️ Column Matcher Diagnostics", expanded=False):
            st.caption("Inspect auto-detected columns per trait:")
            st.json(debug_info)

# -----------------------------------------------------------------------------
# TAB 2: DNA ASM FUNCTIONAL CONVERGENCE (CpG TSS < 5 kb ∩ GTEx Brain eQTLs)
# -----------------------------------------------------------------------------
with tab_asm_convergence:
    st.subheader("🧬 DNA ASM Cross-Trait Recurrence & Functional Convergence")
    st.markdown(
        "Integrates **DNA Allele-Specific Methylation (ASM)** Tag & LD Proxy variants across GWAS traits, "
        "evaluating convergence between **proximal GENCODE v19 gene promoters ($< 5\\text{ kb}$ from CpG)** "
        "and **GTEx v10 human brain cis-eQTL target genes (13 regions)**."
    )

    if discover_trait_comparison_files is None or run_cross_trait_asm_analysis is None:
        st.error("Could not import `scripts/interactive_cross_trait_asm.py`.")
    else:
        asm_trait_files = discover_trait_comparison_files(OVERLAP_DIR)
        if not asm_trait_files:
            st.warning(
                "No `dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv` files found in `results/methylation/`.\n\n"
                "Run `snakemake --cores 2 methylation` (or `--profile profiles/slurm methylation`) first."
            )
        else:
            asm_traits = sorted(asm_trait_files.keys())
            n_asm_traits = len(asm_traits)

            ctrl_col1, ctrl_col2, ctrl_col3, ctrl_col4 = st.columns([1.4, 1.1, 1.0, 1.2])
            with ctrl_col1:
                asm_fdr_label = st.radio(
                    "Multiple-Testing Method:",
                    options=["Storey's q-value (pi0-adjusted)", "Benjamini-Hochberg FDR"],
                    index=0,
                    key="asm_fdr_method",
                )
                asm_method = "storey" if "Storey" in asm_fdr_label else "bh"

            with ctrl_col2:
                asm_k = st.slider(
                    "Min Trait Recurrence (k):",
                    min_value=1,
                    max_value=max(1, n_asm_traits),
                    value=min(2, n_asm_traits),
                    step=1,
                    key="asm_min_k",
                )

            with ctrl_col3:
                asm_alpha = st.number_input(
                    "FDR / q Threshold:",
                    min_value=0.0001,
                    max_value=0.25,
                    value=0.05,
                    step=0.01,
                    format="%.3f",
                    key="asm_alpha",
                )

            with ctrl_col4:
                asm_only_concordant = st.checkbox(
                    "Only Concordant Effector Genes (TSS < 5 kb == GTEx Brain eQTL)",
                    value=False,
                    key="asm_only_conc",
                )

            @st.cache_data(show_spinner="Computing cross-trait DNA ASM FDR & functional convergence...")
            def cached_asm_analysis(input_dir, method, threshold, min_k):
                return run_cross_trait_asm_analysis(
                    input_dir=input_dir,
                    method=method,
                    threshold=threshold,
                    min_traits=min_k,
                    output_path=str(PROJECT_ROOT / "results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv"),
                    save_output=True,
                )

            asm_wide, asm_sig_any, asm_all_long, asm_trait_stats = cached_asm_analysis(
                OVERLAP_DIR, asm_method, float(asm_alpha), int(asm_k)
            )

            display_wide = asm_wide.copy()
            if asm_only_concordant and not display_wide.empty:
                display_wide = display_wide[display_wide["CONCORDANT_GENES"] != "None"].copy()

            # KPI Metrics Row
            kpi1, kpi2, kpi3, kpi4 = st.columns(4)
            kpi1.metric(f"Recurrent ASM SNPs (k ≥ {asm_k})", f"{len(display_wide):,}")
            n_conc_snps = (
                int((display_wide["CONCORDANT_GENES"] != "None").sum())
                if not display_wide.empty
                else 0
            )
            kpi2.metric("Concordant SNPs (TSS == GTEx)", f"{n_conc_snps:,}")
            unique_conc_genes = set()
            if not display_wide.empty:
                for g_str in display_wide["CONCORDANT_GENES"]:
                    if g_str and g_str != "None":
                        for g in str(g_str).split(","):
                            if g.strip():
                                unique_conc_genes.add(g.strip())
            kpi3.metric("Unique Convergent Brain Genes", f"{len(unique_conc_genes):,}")
            kpi4.metric("Evaluated GWAS Traits", f"{n_asm_traits}")

            st.markdown("---")

            # Row 1: Recurrence Barplot & Per-Trait Summary Table
            viz_col1, viz_col2 = st.columns([1.2, 1.2])
            with viz_col1:
                st.markdown("#### 📊 Variant Recurrence Distribution Across Traits")
                if not asm_sig_any.empty:
                    rec_counts = (
                        asm_sig_any.assign(
                            Concordance_Status=np.where(
                                asm_sig_any["CONCORDANT_GENES"] != "None",
                                "Concordant (TSS < 5 kb == GTEx)",
                                "Other ASM Hit",
                            )
                        )
                        .groupby(["N_TRAITS", "Concordance_Status"], as_index=False)
                        .size()
                        .rename(columns={"size": "SNP_Count"})
                    )
                    rec_fig = px.bar(
                        rec_counts,
                        x="N_TRAITS",
                        y="SNP_Count",
                        color="Concordance_Status",
                        barmode="stack",
                        text_auto=True,
                        color_discrete_map={
                            "Concordant (TSS < 5 kb == GTEx)": "#d95f02",
                            "Other ASM Hit": "#7570b3",
                        },
                        labels={
                            "N_TRAITS": "Number of Significant Traits (k)",
                            "SNP_Count": "Number of Unique DNA ASM SNPs",
                        },
                    )
                    rec_fig.update_layout(
                        height=340,
                        xaxis=dict(dtick=1),
                        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                        margin=dict(l=20, r=20, t=40, b=20),
                    )
                    st.plotly_chart(rec_fig, width="stretch")
                else:
                    st.info("No variants reached significance at the current threshold.")

            with viz_col2:
                st.markdown("#### 📋 Per-Trait Candidate FDR Summary")
                if asm_trait_stats:
                    st.dataframe(pd.DataFrame(asm_trait_stats), width="stretch", height=280)

            st.markdown("---")

            # Row 2: Cross-Trait Effect Heatmap
            st.markdown(f"#### 🔥 Cross-Trait Effect Size Heatmap (`BETA`, $k \\ge {asm_k}$)")
            if not display_wide.empty:
                max_plot_snps = st.slider(
                    "Max Top Recurrent Variants to Display in Heatmap:",
                    min_value=10,
                    max_value=min(150, max(10, len(display_wide))),
                    value=min(40, len(display_wide)),
                    step=5,
                    key="asm_max_heatmap_snps",
                )
                top_snps_df = display_wide.head(max_plot_snps).copy()
                top_snp_list = top_snps_df["SNP"].tolist()

                # Pull complete BETA and ADJ_P across all traits from asm_all_long so even sub-threshold directions are visible
                sub_long = asm_all_long[asm_all_long["QUERY_SNP"].isin(top_snp_list)].copy()
                beta_piv = (
                    sub_long.pivot_table(
                        index="TRAIT", columns="QUERY_SNP", values="GWAS_BETA", aggfunc="first"
                    )
                    .reindex(columns=top_snp_list)
                    .fillna(0.0)
                )

                col_labels_asm = []
                for _, r_row in top_snps_df.iterrows():
                    snp_id = r_row["SNP"]
                    c_gene = r_row["CONCORDANT_GENES"]
                    t_gene = r_row["TSS_GENES_5KB"]
                    gene_tag = f"★{c_gene}" if c_gene != "None" else (t_gene if t_gene != "None" else "Intergenic")
                    if len(gene_tag) > 22:
                        gene_tag = gene_tag[:19] + "..."
                    col_labels_asm.append(f"<b>{snp_id}</b><br>({gene_tag})")

                abs_b = np.abs(beta_piv.values)
                b_max = float(np.percentile(abs_b[abs_b > 0], 95)) if (abs_b > 0).any() else 0.05
                b_max = max(b_max, 0.02)

                asm_hm_fig = px.imshow(
                    beta_piv,
                    x=col_labels_asm,
                    y=beta_piv.index,
                    color_continuous_scale="RdBu_r",
                    range_color=[-b_max, b_max],
                    aspect="auto",
                    labels=dict(
                        x="DNA ASM Variant (★ = Concordant TSS & Brain eQTL Gene)",
                        y="GWAS Trait",
                        color="GWAS BETA",
                    ),
                )
                asm_hm_fig.update_layout(
                    height=440,
                    xaxis_tickangle=-45,
                    margin=dict(l=40, r=40, t=20, b=100),
                )
                st.plotly_chart(asm_hm_fig, width="stretch")
            else:
                st.info("No variants match the current recurrence and concordance filters.")

            st.markdown("---")

            # Row 3: Interactive Concordance & Recurrence Table + One-Click Download
            tbl_hdr_left, tbl_hdr_right = st.columns([3, 1])
            with tbl_hdr_left:
                st.markdown("#### 🧬 Recurrent DNA ASM Variants, CpG Sites, Proximal TSS (< 5 kb) & GTEx Brain eQTLs")
            with tbl_hdr_right:
                tsv_export_str = display_wide.to_csv(sep="\t", index=False, na_rep="NA")
                st.download_button(
                    label="📥 Download Filtered TSV",
                    data=tsv_export_str,
                    file_name=f"cross_trait_asm_k{asm_k}_{asm_method}.tsv",
                    mime="text/tab-separated-values",
                    key="asm_download_btn",
                )

            st.dataframe(display_wide, width="stretch", height=420)

# -----------------------------------------------------------------------------
# TAB 3: IN-DASHBOARD RNA m6A COMPARISON GENERATOR
# -----------------------------------------------------------------------------
with tab_generate:
    st.subheader("⚙️ Generate New RNA m6A Comparison & Effector Gene Annotation")
    st.markdown(
        "Execute **RNA $m^6\\text{A}$** candidate variant aggregation directly from the dashboard. "
        "Filter by Storey's $q$-value or Benjamini-Hochberg (BH) FDR, set recurrence thresholds ($k$), "
        "and automatically run **Ensembl VEP & GTEx v10 brain eQTL** annotation (no CpG or TSS lookup needed)."
    )

    if discover_available_files is None:
        st.error("Error: Could not import `scripts/compare_overlap_snps.py`. Ensure the `scripts/` directory is present.")
    else:
        all_records = discover_available_files()
        active_records = [r for r in all_records if r["meth"] in ("RNA_m6A", "RNAm6A")]

        if not active_records:
            st.warning(f"No RNA $m^6\\text{{A}}$ overlap tables (`RNA_m6A_*_overlap_fullheaders.txt`) found in `{OVERLAP_DIR}`.")
            st.info(
                "Run the Tier 1 methylation overlap Snakemake step first:\n\n"
                "```bash\nsnakemake --profile profiles/slurm methylation\n```"
            )
        else:
            st.success(f"Discovered **{len(active_records)}** RNA $m^6\\text{{A}}$ GWAS trait overlap files in `{OVERLAP_DIR}`.")

            available_traits = sorted(list({r["trait"] for r in active_records}))

            fdr_method_label = st.radio(
                "FDR Correction Procedure:",
                options=["Storey's q-value (Cubic smoothing spline)", "Benjamini-Hochberg (fdr_bh)"],
                index=0,
                help="Storey q-value models the true proportion of null variants (pi0) across candidate RNA m6A SNPs."
            )
            chosen_fdr_method = "storey" if "Storey" in fdr_method_label else "bh"

            selected_traits = st.multiselect(
                "Select GWAS Phenotypes to Compare (RNA m6A):",
                options=available_traits,
                default=available_traits,
                help="Choose at least 2 phenotypes to assess RNA m6A pleiotropic recurrence."
            )

            num_traits_selected = len(selected_traits)

            if num_traits_selected < 2:
                st.warning("Please select at least 2 phenotypes to run a cross-trait RNA m6A comparison.")
            else:
                col_k1, col_k2 = st.columns([1, 1])
                with col_k1:
                    min_k_input = st.slider(
                        "Minimum Significant Phenotypes (k):",
                        min_value=1,
                        max_value=num_traits_selected,
                        value=min(2, num_traits_selected),
                        step=1,
                        help="Only RNA m6A variants reaching FDR <= 0.05 in at least k phenotypes will be retained."
                    )

                with col_k2:
                    default_out_name = f"rna_m6a_{num_traits_selected}_k{min_k_input}_{chosen_fdr_method}_sig_snps.tsv"
                    custom_out_name = st.text_input(
                        "Output TSV Filename:",
                        value=default_out_name,
                        help="Filename for the RNA m6A recurrence table in results/methylation/comparisons/"
                    )

                col_opt1, col_opt2 = st.columns(2)
                with col_opt1:
                    do_annotation = st.checkbox(
                        "Run Brain Effector Gene Annotation (Ensembl VEP + GTEx Brain v10)",
                        value=True,
                        help="Queries Ensembl VEP for transcript consequences and GTEx v10 for cis-eQTLs in 13 brain regions."
                    )
                with col_opt2:
                    use_api_cache = st.checkbox(
                        "Use Persistent Local API Cache (.api_cache.json)",
                        value=True,
                        help="Speeds up execution and prevents redundant queries for previously annotated variants."
                    )

                if st.button("🚀 Run RNA m6A Comparison & Effector Annotation", type="primary"):
                    target_records = [r for r in active_records if r["trait"] in selected_traits]
                    output_file_path = os.path.join(str(COMPARISONS_DIR), custom_out_name)

                    with st.status("Executing RNA m6A Epigenetic Comparison Pipeline...", expanded=True) as status_box:
                        st.write(f"1️⃣ Merging {len(target_records)} RNA m6A phenotype datasets & applying {chosen_fdr_method.upper()} FDR...")
                        cmp_out, final_df = compare_and_filter(
                            target_records,
                            output_file_path,
                            min_k=min_k_input,
                            method=chosen_fdr_method
                        )

                        if final_df.empty:
                            status_box.update(label="No RNA m6A variants met the recurrence criteria.", state="error")
                            st.error(f"0 variants reached FDR <= 0.05 in >= {min_k_input} phenotypes. Try lowering k.")
                        else:
                            st.write(f"✓ Found {len(final_df):,} pleiotropic RNA m6A variants meeting $k \\ge {min_k_input}$.")

                            if do_annotation and annotate_file is not None:
                                st.write("2️⃣ Querying Ensembl VEP & GTEx Portal v10 (13 Brain Tissues)...")
                                annotated_tsv = cmp_out.replace(".tsv", "_annotated.tsv")
                                annotate_file(
                                    input_path=cmp_out,
                                    output_path=annotated_tsv,
                                    create_summary=True,
                                    max_workers=8,
                                    cache_path=str(DEFAULT_CACHE_FILE),
                                    use_cache=use_api_cache
                                )
                                st.write("✓ Ensembl VEP & GTEx Brain annotation complete!")

                            # Clear Streamlit cache to load the new dataset
                            st.cache_data.clear()
                            status_box.update(label="RNA m6A Pipeline Completed Successfully!", state="complete", expanded=False)
                            st.success(
                                f"🎉 Run **{custom_out_name}** generated successfully! "
                                "Switch to the **🧪 RNA m6A Methylation Explorer** tab to inspect the interactive heatmap and effector targets."
                            )