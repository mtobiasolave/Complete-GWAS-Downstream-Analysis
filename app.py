"""
Methylation x GWAS Explorer: Streamlit Front-End for the DNA & RNA Methylation Pipeline.

Configures and launches Snakemake runs (local background process or Slurm controller job),
explores the cross-trait pleiotropy tables written to results/dna/ and results/rna/, and
runs and visualizes LDSC genetic correlations (results/ldsc_rg/).

Run from the repository root with snakemake9_env activated:
    streamlit run app.py
"""

import gzip
import json
import os
import re
import shutil
import signal
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import yaml

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from common_utils import (  # noqa: E402
    detect_n_column,
    discover_gwas,
    munge_column_flags,
    read_gwas_header,
)
from ldsc_rg_summary import (  # noqa: E402
    PAIR_SEP,
    RG_RESULT_COLUMNS,
    RG_TABLE_COLUMNS,
    add_fdr,
    pair_stem,
    parse_rg_log,
    rg_pairs,
)

RUNS_DIR = ROOT / "logs" / "runs"
DNA_TABLE = ROOT / "results" / "dna" / "pleiotropy_dna.tsv"
RNA_TABLE = ROOT / "results" / "rna" / "pleiotropy_rna.tsv"
MUNGED_DIR = ROOT / "results" / "munged"
MUNGE_LOG_DIR = ROOT / "logs" / "munge"
LDSC_WRITE_RE = re.compile(r"Writing summary statistics for [\d,]+ SNPs \(([\d,]+) with nonmissing")
LDSC_METRIC_RE = {
    "Mean chi²": re.compile(r"Mean chi\^2 = ([\d.]+)"),
    "λGC": re.compile(r"Lambda GC = ([\d.]+)"),
}
RG_DIR = ROOT / "results" / "ldsc_rg"
RG_STDOUT_DIR = RG_DIR / "stdout"
RG_MODES = ["All-vs-All Matrix", "Target vs. Panel (Hub & Spokes)"]
# Diverging blue <-> red with a neutral gray midpoint at rg = 0
RG_COLORSCALE = [[0.0, "#184f95"], [0.25, "#6da7ec"], [0.5, "#f0efec"], [0.75, "#ec8f8e"], [1.0, "#a32c2b"]]
RG_SIG_COLOR = "#2a78d6"
RG_NS_COLOR = "#898781"

TARGET_MODES = {
    "Both (DNA + RNA)": "all",
    "DNA Methylation Only": "dna",
    "RNA Methylation Only": "rna",
}
ENGINES = {
    "Local (Direct Snakemake)": "local",
    "HPC (Slurm / sbatch)": "slurm",
}
PROGRESS_RE = re.compile(r"(\d+) of (\d+) steps \((\d+(?:\.\d+)?)%\) done")
IS_WINDOWS = os.name == "nt"

COLUMN_CONFIG = {
    "SNP": st.column_config.TextColumn("rsID"),
    "N_TRAITS": st.column_config.NumberColumn("N traits (all)"),
    "N_SELECTED": st.column_config.NumberColumn("N traits (selected)"),
    "TRAITS_PRESENT": st.column_config.TextColumn("Traits present"),
    "CONCORDANT": st.column_config.CheckboxColumn("Concordant"),
    "CONCORDANT_GENES": st.column_config.TextColumn("Concordant gene(s)"),
    "GENES": st.column_config.TextColumn(
        "Gene(s) (Ensembl, overlapping)", help="Ensembl GRCh37 genes whose transcripts the SNP falls inside"
    ),
    "CPG_GENES_5KB": st.column_config.TextColumn("CpG gene(s) (TSS ≤ 5 kb)"),
    "GTEX_BRAIN_GENES": st.column_config.TextColumn("GTEx brain eQTL gene(s)"),
    "CPG_IDS": st.column_config.TextColumn("CpG ID(s)"),
    "GTEX_TISSUES": st.column_config.TextColumn("GTEx brain tissues"),
    "MIN_FDR_BH": st.column_config.NumberColumn("Min BH FDR", format="%.2e"),
    "MIN_FDR_BH_DNA": st.column_config.NumberColumn("Min BH FDR (DNA)", format="%.2e"),
    "MIN_FDR_BH_RNA": st.column_config.NumberColumn("Min BH FDR (RNA)", format="%.2e"),
}

st.set_page_config(page_title="Methylation x GWAS Explorer", page_icon="🧬", layout="wide")


# =====================================================================
# 1. Configuration & Trait Discovery
# =====================================================================

def load_config() -> dict:
    """Merge config.yaml with config_local.yaml (same precedence as the Snakefile)."""
    config = {}
    for name in ("config.yaml", "config_local.yaml"):
        path = ROOT / name
        if path.exists():
            config.update(yaml.safe_load(path.read_text(encoding="utf-8")) or {})
    return config


def available_traits(config: dict) -> dict:
    """Map trait name -> GWAS file for every summary statistics file the pipeline can use."""
    traits = discover_gwas(ROOT / config.get("raw_data_dir", "data/raw_data"))
    traits.update(config.get("gwas_sumstats") or {})
    return traits


# =====================================================================
# 2. Asynchronous Run Management (Local Process or Slurm Controller Job)
# =====================================================================

@st.cache_resource
def local_processes() -> dict:
    """Popen handles of local runs started by this server (run id -> Popen)."""
    return {}


def latest_run():
    runs = sorted(RUNS_DIR.glob("run_*.json"))
    return json.loads(runs[-1].read_text(encoding="utf-8")) if runs else None


def start_run(
    target: str, traits: list, engine: str, cores: int, extra_config: list = None, forcerun: str = None
) -> dict:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    log_rel = f"logs/runs/{run_id}.log"
    status_rel = f"logs/runs/{run_id}.exit"
    traits_arg = ",".join(traits)
    run = {
        "id": run_id,
        "target": target,
        "traits": traits,
        "engine": engine,
        "log": log_rel,
        "status_file": status_rel,
        "started": datetime.now().isoformat(timespec="seconds"),
        "forcerun": forcerun,
        "extra_config": extra_config or [],
    }
    pipeline_opts = [opt for item in extra_config or [] for opt in ("--extra-config", item)]
    if forcerun:
        pipeline_opts += ["--forcerun", forcerun]

    if engine == "local":
        cmd = [
            sys.executable, "scripts/run_pipeline.py",
            "--target", target, "--traits", traits_arg, "--engine", "local",
            "--cores", str(cores), "--status-file", status_rel, *pipeline_opts,
        ]
        # Detach from the Streamlit server so the run survives page reloads
        detach = (
            {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
            if IS_WINDOWS
            else {"start_new_session": True}
        )
        with open(ROOT / log_rel, "w", encoding="utf-8") as log_file:
            proc = subprocess.Popen(
                cmd, cwd=ROOT, stdout=log_file, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, **detach,
            )
        run["pid"] = proc.pid
        local_processes()[run_id] = proc
    else:
        cmd = [
            "sbatch", "--parsable", f"--output={log_rel}",
            "scripts/run_pipeline.sbatch", target, traits_arg, status_rel, *pipeline_opts,
        ]
        result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True)
        run["job_id"] = result.stdout.strip().split(";")[0]

    (RUNS_DIR / f"{run_id}.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    return run


def local_process_alive(run: dict) -> bool:
    proc = local_processes().get(run["id"])
    if proc is not None:
        return proc.poll() is None
    # Run started by an earlier server session: fall back to a PID check
    if IS_WINDOWS:
        out = subprocess.run(
            ["tasklist", "/FI", f"PID eq {run['pid']}", "/NH"], capture_output=True, text=True
        )
        return str(run["pid"]) in out.stdout
    try:
        os.kill(run["pid"], 0)
        return True
    except OSError:
        return False


def run_status(run: dict) -> str:
    """Return 'running', 'pending', 'succeeded', 'failed', 'cancelled' or 'unknown'."""
    status_path = ROOT / run["status_file"]
    if status_path.exists():
        code = status_path.read_text(encoding="utf-8").strip()
        return {"0": "succeeded", "cancelled": "cancelled"}.get(code, "failed")

    if run["engine"] == "slurm":
        try:
            out = subprocess.run(
                ["squeue", "-h", "-j", run["job_id"], "-o", "%T"], capture_output=True, text=True
            )
        except OSError:
            return "unknown"
        state = out.stdout.strip().upper()
        if state == "PENDING":
            return "pending"
        # A controller that left the queue without an exit code was killed (time limit, node failure)
        return "running" if state else "failed"

    return "running" if local_process_alive(run) else "failed"


def cancel_run(run: dict) -> None:
    if run["engine"] == "slurm":
        subprocess.run(["scancel", run["job_id"]], capture_output=True)
    elif IS_WINDOWS:
        subprocess.run(["taskkill", "/PID", str(run["pid"]), "/T", "/F"], capture_output=True)
        # A forced kill on Windows leaves the Snakemake lock behind
        subprocess.run([sys.executable, "-m", "snakemake", "--unlock"], cwd=ROOT, capture_output=True)
    else:
        try:
            os.killpg(run["pid"], signal.SIGTERM)
        except OSError:
            pass
    (ROOT / run["status_file"]).write_text("cancelled", encoding="utf-8")


def read_log_tail(path: Path, max_bytes: int = 64_000) -> str:
    if not path.exists():
        return ""
    with open(path, "rb") as f:
        f.seek(0, os.SEEK_END)
        f.seek(max(0, f.tell() - max_bytes))
        return f.read().decode("utf-8", errors="replace")


def run_monitor() -> None:
    """Live status, progress and log tail of the most recent run."""
    run = latest_run()
    if run is None:
        st.info("No pipeline runs yet. Configure a run in the sidebar.")
        return

    status = run_status(run)
    active = status in ("running", "pending")
    icon = {"running": "⏳", "pending": "🕒", "succeeded": "✅", "failed": "❌", "cancelled": "⛔"}.get(status, "❔")
    engine_label = "Slurm job " + run["job_id"] if run["engine"] == "slurm" else "local process"

    log_text = read_log_tail(ROOT / run["log"])
    matches = PROGRESS_RE.findall(log_text)
    done, total, pct = (int(matches[-1][0]), int(matches[-1][1]), float(matches[-1][2])) if matches else (0, 0, 0.0)
    if status == "succeeded":
        pct = 100.0

    with st.expander(
        f"{icon} Last run `{run['id']}` — {status} · target `{run['target']}` · {engine_label}",
        expanded=active or status == "failed",
    ):
        st.progress(pct / 100, text=f"{done} of {total} steps ({pct:.0f}%)" if total else "Waiting for Snakemake to start…")
        st.caption(f"Traits: {', '.join(run['traits'])} · Started {run['started']} · Log: `{run['log']}`")
        st.code("\n".join(log_text.splitlines()[-40:]) or "(log is empty)", language="log", height=300)
        if active and st.button("⛔ Cancel run", key=f"cancel_{run['id']}"):
            cancel_run(run)
            st.rerun(scope="app")

    # Reload the whole page once a run this session watched has finished, so result tables refresh
    if active:
        st.session_state["watching_run"] = run["id"]
    elif st.session_state.pop("watching_run", None) == run["id"]:
        st.rerun(scope="app")


# =====================================================================
# 3. Cached Table Loaders & In-Memory Filters
# =====================================================================

@st.cache_data(show_spinner=False)
def load_pleiotropy(path: str, mtime: float) -> pd.DataFrame:
    """Load a pleiotropy TSV; `mtime` keys the cache so a finished run invalidates it."""
    df = pd.read_csv(path, sep="\t")
    for col in df.columns:
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].fillna("").astype(str)
    if "CONCORDANT" in df.columns:
        df["CONCORDANT"] = df["CONCORDANT"].astype(str).str.lower().eq("true")
    df["TRAIT_LIST"] = df["TRAITS_PRESENT"].str.split(", ")
    return df


def table_or_none(path: Path):
    return load_pleiotropy(str(path), path.stat().st_mtime) if path.exists() else None


@st.cache_data(show_spinner=False)
def load_multiomics(dna_mtime: float, rna_mtime: float) -> pd.DataFrame:
    """SNPs BH-significant in both the DNA and RNA methylation tracks."""
    dna = load_pleiotropy(str(DNA_TABLE), dna_mtime)
    rna = load_pleiotropy(str(RNA_TABLE), rna_mtime)
    merged = pd.merge(
        dna[["SNP", "TRAIT_LIST", "CONCORDANT", "CONCORDANT_GENES", "CPG_GENES_5KB", "GTEX_BRAIN_GENES", "MIN_FDR_BH"]],
        rna[["SNP", "TRAIT_LIST", "MIN_FDR_BH"]],
        on="SNP",
        suffixes=("_DNA", "_RNA"),
    )
    merged["TRAITS_DNA"] = merged["TRAIT_LIST_DNA"].str.join(", ")
    merged["TRAITS_RNA"] = merged["TRAIT_LIST_RNA"].str.join(", ")
    merged["TRAIT_LIST"] = [
        sorted(set(d) | set(r)) for d, r in zip(merged["TRAIT_LIST_DNA"], merged["TRAIT_LIST_RNA"])
    ]
    merged["TRAITS_PRESENT"] = merged["TRAIT_LIST"].str.join(", ")
    return merged.drop(columns=["TRAIT_LIST_DNA", "TRAIT_LIST_RNA"])


def trait_sharing_filter(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """Count presence among user-selected traits and keep SNPs present in >= N of them."""
    all_traits = sorted({t for traits in df["TRAIT_LIST"] for t in traits})
    selected = set(st.multiselect(
        "Traits to count",
        all_traits,
        default=all_traits,
        key=f"{key}_traits",
        help="Traits with at least one BH-significant SNP in the latest run.",
    ))
    max_n = max(len(selected), 1)
    min_n = (
        st.slider("Present in ≥ N selected traits", 1, max_n, 1, key=f"{key}_min_n")
        if max_n > 1
        else 1
    )
    n_selected = df["TRAIT_LIST"].map(lambda traits: len(selected.intersection(traits)))
    return df.assign(N_SELECTED=n_selected)[n_selected >= min_n]


def search_filter(df: pd.DataFrame, query: str, columns: list) -> pd.DataFrame:
    if not query:
        return df
    mask = pd.Series(False, index=df.index)
    for col in columns:
        mask |= df[col].astype(str).str.contains(query.strip(), case=False, regex=False)
    return df[mask]


def show_table(df: pd.DataFrame, columns: list, file_name: str, key: str) -> None:
    sorted_df = df.sort_values(["N_SELECTED", columns[-1]], ascending=[False, True])
    st.dataframe(
        sorted_df[columns], hide_index=True, width="stretch", height=520, column_config=COLUMN_CONFIG
    )
    st.download_button(
        "⬇️ Export filtered table (TSV)",
        sorted_df[columns].to_csv(sep="\t", index=False).encode("utf-8"),
        file_name=file_name,
        mime="text/tab-separated-values",
        key=f"{key}_download",
    )


def missing_results(track: str, target: str) -> None:
    st.info(
        f"No {track} methylation results yet. Run the pipeline with "
        f"**Both (DNA + RNA)** or **{target}** from the sidebar."
    )


# =====================================================================
# 4. GWAS Sumstats Munging (LDSC munge_sumstats.py)
# =====================================================================

@st.cache_data(show_spinner=False)
def gwas_header(path: str, mtime: float) -> list:
    return read_gwas_header(path)


def munged_signature() -> tuple:
    """Names, sizes and mtimes of results/munged/ files; keys the inventory cache."""
    if not MUNGED_DIR.is_dir():
        return ()
    return tuple(
        (p.name, p.stat().st_size, p.stat().st_mtime) for p in sorted(MUNGED_DIR.iterdir()) if p.is_file()
    )


@st.cache_data(show_spinner=False)
def munged_inventory(signature: tuple) -> pd.DataFrame:
    """One row per results/munged/*.sumstats.gz with summary metrics parsed from its LDSC log."""
    rows = []
    for name, size, mtime in signature:
        if not name.endswith(".sumstats.gz"):
            continue
        trait = name[: -len(".sumstats.gz")]
        log_path = MUNGED_DIR / f"{trait}.log"
        log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
        written = LDSC_WRITE_RE.search(log_text)
        row = {
            "Trait": trait,
            "Size (MB)": round(size / 1e6, 2),
            "Last modified": datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M"),
            "SNPs with Z": int(written.group(1).replace(",", "")) if written else None,
        }
        for label, pattern in LDSC_METRIC_RE.items():
            match = pattern.search(log_text)
            row[label] = float(match.group(1)) if match else None
        rows.append(row)
    return pd.DataFrame(rows, columns=["Trait", "Size (MB)", "Last modified", "SNPs with Z", *LDSC_METRIC_RE])


def file_head(path: Path, n_lines: int = 5) -> str:
    open_fn = gzip.open if path.suffix == ".gz" else open
    with open_fn(path, "rt", encoding="utf-8", errors="replace") as f:
        return "".join(line for _, line in zip(range(n_lines), f))


def munge_trait_status(trait: str, run: dict, status: str) -> str:
    """Per-trait state within a munge run, from output/log mtimes relative to the run start."""
    started = datetime.fromisoformat(run["started"]).timestamp()
    output = MUNGED_DIR / f"{trait}.sumstats.gz"
    live_log = MUNGE_LOG_DIR / f"{trait}.log"
    if output.exists() and output.stat().st_mtime >= started:
        return "✅ Done"
    if live_log.exists() and live_log.stat().st_mtime >= started:
        if status in ("running", "pending"):
            return "⏳ Running"
        return "⛔ Cancelled" if status == "cancelled" else "❌ Failed"
    if output.exists() and not run.get("forcerun"):
        return "✅ Done (already up to date)"
    if status in ("running", "pending"):
        return "🕒 Queued"
    return "⛔ Cancelled" if status == "cancelled" else "❌ Failed"


def munge_monitor() -> None:
    """Per-trait status and live log viewer for the latest munging run."""
    run = latest_run()
    if run is None or run["target"] != "munge":
        st.caption("No munging run in progress or recorded as the latest run.")
        return

    status = run_status(run)
    st.markdown(f"**Latest munging run** `{run['id']}` — {status}")
    st.dataframe(
        pd.DataFrame([{"Trait": t, "Status": munge_trait_status(t, run, status)} for t in run["traits"]]),
        hide_index=True,
    )
    trait = st.selectbox("Log viewer", run["traits"], key="munge_log_trait")
    live_log = MUNGE_LOG_DIR / f"{trait}.log"
    log_path = live_log if live_log.exists() else MUNGED_DIR / f"{trait}.log"
    log_text = read_log_tail(log_path)
    st.caption(f"`{log_path.relative_to(ROOT).as_posix()}`")
    st.code(
        "\n".join(log_text.splitlines()[-60:]) or "(no log yet: trait is queued)",
        language="log",
        height=320,
    )


# =====================================================================
# 5. Genetic Correlation (LDSC ldsc.py --rg)
# =====================================================================

def rg_signature() -> tuple:
    """Names and mtimes of results/ldsc_rg/{trait1}__{trait2}.log files; keys the results cache."""
    if not RG_DIR.is_dir():
        return ()
    return tuple((p.name, p.stat().st_mtime) for p in sorted(RG_DIR.glob(f"*{PAIR_SEP}*.log")))


@st.cache_data(show_spinner=False)
def load_rg_results(signature: tuple) -> pd.DataFrame:
    """rg, SE, Z, P (plus h2 and intercepts) parsed from every computed pair's LDSC log."""
    rows = [parse_rg_log(RG_DIR / name) for name, _ in signature]
    return pd.DataFrame(rows, columns=RG_RESULT_COLUMNS)


def rg_view(results: pd.DataFrame, pairs: list) -> pd.DataFrame:
    """Computed pairs of the current selection, with BH FDR across them (the on-screen family)."""
    finished = results[~results["Note"].str.startswith(("incomplete", "unparseable"))]
    wanted = pd.DataFrame(pairs, columns=["Trait 1", "Trait 2"])
    return add_fdr(wanted.merge(finished.drop(columns="FDR"), on=["Trait 1", "Trait 2"]))


def rg_heatmap(view: pd.DataFrame, traits: list, alpha: float) -> go.Figure:
    """Symmetric rg matrix on a diverging -1..+1 scale; '*' marks BH FDR < alpha."""
    traits = sorted(traits)
    pos = {t: i for i, t in enumerate(traits)}
    n = len(traits)
    z = np.full((n, n), np.nan)
    text = np.full((n, n), "", dtype=object)
    hover = np.full((n, n), "Not computed yet", dtype=object)
    for i, trait in enumerate(traits):
        hover[i, i] = trait
    for row in view.to_dict("records"):
        t1, t2 = row["Trait 1"], row["Trait 2"]
        if pd.isna(row["rg"]):
            cell, tip = "NA", f"{t1} × {t2}<br>rg not estimated: {row['Note']}"
        else:
            cell = f"{row['rg']:.2f}" + ("*" if row["FDR"] < alpha else "")
            tip = (
                f"{t1} × {t2}<br>rg = {row['rg']:.3f} (SE {row['SE']:.3f})"
                f"<br>Z = {row['Z']:.2f} · P = {row['P-value']:.2e} · FDR = {row['FDR']:.2e}"
            )
        for a, b in ((pos[t1], pos[t2]), (pos[t2], pos[t1])):
            z[a, b], text[a, b], hover[a, b] = row["rg"], cell, tip

    fig = go.Figure(go.Heatmap(
        z=z, x=traits, y=traits, zmin=-1, zmax=1, colorscale=RG_COLORSCALE,
        text=text, texttemplate="%{text}", customdata=hover,
        hovertemplate="%{customdata}<extra></extra>", xgap=2, ygap=2,
        colorbar=dict(title="rg", tickvals=[-1, -0.5, 0, 0.5, 1]),
    ))
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(height=min(220 + 48 * n, 900), margin=dict(t=20, l=10, r=10, b=10))
    return fig


def rg_forest(view: pd.DataFrame, target: str, alpha: float) -> go.Figure:
    """rg of each panel trait with the target, 95% CI (rg ± 1.96 SE), sorted by rg."""
    df = view.dropna(subset=["rg"]).copy()
    df["Trait"] = np.where(df["Trait 1"] == target, df["Trait 2"], df["Trait 1"])
    df["CI low"], df["CI high"] = df["rg"] - 1.96 * df["SE"], df["rg"] + 1.96 * df["SE"]
    df = df.sort_values("rg")
    significant = df["FDR"] < alpha

    fig = go.Figure()
    for is_sig, name, color in ((False, "Not significant", RG_NS_COLOR), (True, f"BH FDR < {alpha:g}", RG_SIG_COLOR)):
        part = df[significant == is_sig]
        if part.empty:
            continue
        fig.add_trace(go.Scatter(
            x=part["rg"], y=part["Trait"], mode="markers", name=name,
            marker=dict(color=color, size=10),
            error_x=dict(type="data", array=1.96 * part["SE"], color=color, thickness=2, width=4),
            customdata=part[["SE", "CI low", "CI high", "Z", "P-value", "FDR"]].to_numpy(),
            hovertemplate=(
                "%{y}<br>rg = %{x:.3f} (SE %{customdata[0]:.3f})"
                "<br>95% CI %{customdata[1]:.3f} to %{customdata[2]:.3f}"
                "<br>Z = %{customdata[3]:.2f} · P = %{customdata[4]:.2e} · FDR = %{customdata[5]:.2e}<extra></extra>"
            ),
        ))
    fig.add_vline(x=0, line_dash="dash", line_width=1, line_color=RG_NS_COLOR)
    fig.update_yaxes(categoryorder="array", categoryarray=df["Trait"].tolist(), title=None)
    fig.update_xaxes(title=f"Genetic correlation with {target} (rg, 95% CI)")
    fig.update_layout(
        height=min(180 + 36 * len(df), 900),
        margin=dict(t=40, l=10, r=10, b=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
    )
    return fig


def rg_run_pairs(run: dict) -> list:
    """Trait pairs of an rg run, from the `rg=` Snakemake config it was launched with."""
    for item in run.get("extra_config", []):
        if item.startswith("rg="):
            rg_cfg = json.loads(item[len("rg="):])
            return rg_pairs(rg_cfg.get("traits") or [], rg_cfg.get("target") or None)
    return []


def rg_pair_status(stem: str, run: dict, status: str) -> str:
    """Per-pair state within an rg run, from output/log mtimes relative to the run start."""
    started = datetime.fromisoformat(run["started"]).timestamp()
    output = RG_DIR / f"{stem}.log"
    live_log = RG_STDOUT_DIR / f"{stem}.log"
    if output.exists() and output.stat().st_mtime >= started:
        return "✅ Done"
    if live_log.exists() and live_log.stat().st_mtime >= started:
        if status in ("running", "pending"):
            return "⏳ Running"
        return "⛔ Cancelled" if status == "cancelled" else "❌ Failed"
    if output.exists() and not run.get("forcerun"):
        return "✅ Done (already up to date)"
    if status in ("running", "pending"):
        return "🕒 Queued"
    return "⛔ Cancelled" if status == "cancelled" else "❌ Failed"


def rg_monitor() -> None:
    """Status spinner, per-pair progress and live LDSC log of the latest rg run."""
    run = latest_run()
    if run is None or run["target"] != "rg":
        st.caption("No genetic correlation run in progress or recorded as the latest run.")
        return

    status = run_status(run)
    stems = [pair_stem(*pair) for pair in rg_run_pairs(run)]
    states = [rg_pair_status(stem, run, status) for stem in stems]
    n_done = sum(state.startswith("✅") for state in states)
    box_state = "running" if status in ("running", "pending") else "complete" if status == "succeeded" else "error"
    with st.status(
        f"LDSC rg run `{run['id']}` — {status} · {n_done} of {len(stems)} pairs done",
        state=box_state,
        expanded=box_state != "complete",
    ):
        st.dataframe(
            pd.DataFrame({"Pair": [s.replace(PAIR_SEP, " × ") for s in stems], "Status": states}),
            hide_index=True,
        )
        if stems:
            stem = st.selectbox(
                "Log viewer", stems, format_func=lambda s: s.replace(PAIR_SEP, " × "), key="rg_log_pair"
            )
            live_log = RG_STDOUT_DIR / f"{stem}.log"
            log_path = live_log if live_log.exists() else RG_DIR / f"{stem}.log"
            st.caption(f"`{log_path.relative_to(ROOT).as_posix()}`")
            st.code(
                "\n".join(read_log_tail(log_path).splitlines()[-60:]) or "(no log yet: pair is queued)",
                language="log",
                height=320,
            )


# =====================================================================
# 6. Sidebar: Execution & Configuration Panel
# =====================================================================

config = load_config()
gwas_files = available_traits(config)
current_run = latest_run()
run_active = current_run is not None and run_status(current_run) in ("running", "pending")

with st.sidebar:
    st.header("⚙️ Run Pipeline")
    trait_names = sorted(gwas_files)
    configured = config.get("traits") or []
    default_traits = [t for t in configured if t in gwas_files] or trait_names
    selected_traits = st.multiselect(
        "GWAS traits",
        trait_names,
        default=default_traits,
        help=f"Summary statistics found in `{config.get('raw_data_dir', 'data/raw_data')}/`.",
    )
    mode_label = st.radio("Target mode", list(TARGET_MODES))

    has_sbatch = shutil.which("sbatch") is not None
    engine_label = st.radio(
        "Execution engine",
        list(ENGINES),
        index=1 if has_sbatch else 0,
        help="Auto-detected: HPC is preselected when `sbatch` is available on this machine.",
    )
    engine = ENGINES[engine_label]
    cores = 1
    if engine == "local":
        cpu_count = os.cpu_count() or 1
        cores = st.number_input("Local cores", 1, cpu_count, min(2, cpu_count))
    elif not has_sbatch:
        st.warning("`sbatch` was not found on this machine; HPC runs are unavailable here.")
    else:
        st.caption("Submits a controller job to `irbio01`, which submits each step via `profiles/slurm`.")

    launch_blocked = run_active or not selected_traits or (engine == "slurm" and not has_sbatch)
    if st.button("🚀 Run pipeline", type="primary", disabled=launch_blocked, width="stretch"):
        try:
            start_run(TARGET_MODES[mode_label], selected_traits, engine, int(cores))
            st.rerun()
        except (OSError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            st.error(f"Could not start the pipeline: {detail}")
    if run_active:
        st.caption("A run is in progress; results below show the last completed outputs.")

    st.divider()
    st.caption("Multiple testing: Benjamini-Hochberg FDR ≤ 5%, applied per trait across the reference SNPs found in its GWAS.")


# =====================================================================
# 7. Main Area: Run Monitor & Interactive Results
# =====================================================================

st.title("🧬 Methylation × GWAS Pleiotropy Explorer")
st.caption(
    "DNA ASM tag SNPs & LD proxies → CpG gene (TSS ≤ 5 kb) vs. GTEx v10 brain eQTL concordance · "
    "RNA m6A SNPs → overlapping Ensembl genes (GRCh37) · BH FDR 5% per trait"
)

st.fragment(run_monitor, run_every=3 if run_active else None)()

tab_dna, tab_rna, tab_multi, tab_munge, tab_rg = st.tabs(
    ["🧬 DNA Methylation", "🔬 RNA Methylation", "🔗 Multi-Omics Overlap", "🧹 Munge Sumstats", "📊 Genetic Correlation"]
)

with tab_dna:
    dna_df = table_or_none(DNA_TABLE)
    if dna_df is None:
        missing_results("DNA", "DNA Methylation Only")
    else:
        col_filters, col_search = st.columns([3, 2])
        with col_filters:
            dna_view = trait_sharing_filter(dna_df, "dna")
        with col_search:
            concordant_only = st.toggle("Only show concordant genes (CpG gene == GTEx brain eQTL)", key="dna_concordant")
            dna_query = st.text_input("Search rsID / gene / CpG", key="dna_search")
        if concordant_only:
            dna_view = dna_view[dna_view["CONCORDANT"]]
        dna_view = search_filter(
            dna_view, dna_query, ["SNP", "CONCORDANT_GENES", "CPG_GENES_5KB", "GTEX_BRAIN_GENES", "CPG_IDS"]
        )

        m1, m2, m3 = st.columns(3)
        m1.metric("SNPs shown", f"{len(dna_view):,}", f"of {len(dna_df):,}", delta_color="off")
        m2.metric("Concordant SNPs shown", f"{int(dna_view['CONCORDANT'].sum()):,}")
        m3.metric("Concordant genes shown", f"{dna_view['CONCORDANT_GENES'].str.split(', ').explode().replace('', pd.NA).nunique():,}")
        show_table(
            dna_view,
            ["SNP", "N_SELECTED", "TRAITS_PRESENT", "CONCORDANT", "CONCORDANT_GENES", "CPG_GENES_5KB",
             "GTEX_BRAIN_GENES", "CPG_IDS", "GTEX_TISSUES", "MIN_FDR_BH"],
            "pleiotropy_dna_filtered.tsv",
            "dna",
        )

with tab_rna:
    rna_df = table_or_none(RNA_TABLE)
    if rna_df is None:
        missing_results("RNA", "RNA Methylation Only")
    else:
        if "GENES" not in rna_df.columns:
            st.info("This RNA table predates the Ensembl gene step. Re-run **RNA Methylation Only** to add the gene column.")
            rna_df = rna_df.assign(GENES="")
        col_filters, col_search = st.columns([3, 2])
        with col_filters:
            rna_view = trait_sharing_filter(rna_df, "rna")
        with col_search:
            rna_query = st.text_input("Search rsID / gene", key="rna_search")
        rna_view = search_filter(rna_view, rna_query, ["SNP", "GENES"])

        st.metric("SNPs shown", f"{len(rna_view):,}", f"of {len(rna_df):,}", delta_color="off")
        show_table(
            rna_view,
            ["SNP", "GENES", "N_SELECTED", "TRAITS_PRESENT", "MIN_FDR_BH"],
            "pleiotropy_rna_filtered.tsv",
            "rna",
        )

with tab_multi:
    if not (DNA_TABLE.exists() and RNA_TABLE.exists()):
        st.info("Both DNA and RNA results are needed. Run the pipeline with **Both (DNA + RNA)**.")
    else:
        multi_df = load_multiomics(DNA_TABLE.stat().st_mtime, RNA_TABLE.stat().st_mtime)
        st.caption(
            "SNPs that are BH-significant in at least one trait in **both** tracks "
            "(DNA ASM tag/proxy SNP and RNA m6A SNP). Trait counts use the union of both tracks."
        )
        if multi_df.empty:
            st.warning("No SNP is significant in both the DNA and RNA methylation tracks.")
        else:
            col_filters, col_search = st.columns([3, 2])
            with col_filters:
                multi_view = trait_sharing_filter(multi_df, "multi")
            with col_search:
                multi_query = st.text_input("Search rsID / gene", key="multi_search")
            multi_view = search_filter(
                multi_view, multi_query, ["SNP", "CONCORDANT_GENES", "CPG_GENES_5KB", "GTEX_BRAIN_GENES"]
            )

            st.metric("Shared DNA ∩ RNA SNPs shown", f"{len(multi_view):,}", f"of {len(multi_df):,}", delta_color="off")
            show_table(
                multi_view,
                ["SNP", "N_SELECTED", "TRAITS_PRESENT", "TRAITS_DNA", "TRAITS_RNA", "CONCORDANT",
                 "CONCORDANT_GENES", "CPG_GENES_5KB", "GTEX_BRAIN_GENES", "MIN_FDR_BH_DNA", "MIN_FDR_BH_RNA"],
                "multiomics_overlap_filtered.tsv",
                "multi",
            )

with tab_munge:
    ldsc_cfg = config.get("ldsc") or {}
    munge_script = ROOT / ldsc_cfg.get("munge_script", "../ldsc/munge_sumstats.py")
    munged_traits = {t for t in trait_names if (MUNGED_DIR / f"{t}.sumstats.gz").exists()}
    st.caption(
        "Harmonize raw GWAS summary statistics with LDSC `munge_sumstats.py` into "
        "`results/munged/{trait}.sumstats.gz`. Runs with the execution engine selected in the sidebar."
    )
    if not munge_script.is_file():
        st.error(f"LDSC munge script not found: `{ldsc_cfg.get('munge_script')}` (config `ldsc.munge_script`).")
    elif "usecols=list(cname_translation.keys())" not in munge_script.read_text(encoding="utf-8", errors="replace"):
        st.warning(
            "This LDSC copy looks unpatched for Python 3. Apply the bundled fixes with "
            "`git -C ../ldsc apply \"$PWD/patches/ldsc_python3.patch\"` (run from the repository root)."
        )

    include_munged = st.checkbox("Include already munged traits (re-munge)", key="munge_include_done")
    munge_options = trait_names if include_munged else [t for t in trait_names if t not in munged_traits]
    munge_traits = st.multiselect(
        "Raw GWAS traits to munge",
        munge_options,
        key="munge_traits",
        placeholder="Choose traits" if munge_options else "All traits are already munged",
    )

    n_problems = {}
    n_override_out, n_col_override_out = {}, {}
    with st.expander("⚙️ Munging options", expanded=bool(munge_traits)):
        merge_alleles = st.text_input(
            "Reference alleles file (--merge-alleles)",
            value=ldsc_cfg.get("merge_alleles", "data/reference/w_hm3.snplist"),
            help="Path relative to the repository root.",
        )
        merge_ok = (ROOT / merge_alleles).is_file()
        if not merge_ok:
            st.error(f"Reference alleles file not found: `{merge_alleles}`")

        st.markdown(
            "**Sample size per trait.** Precedence: manual N > N column override > detected N column."
        )
        n_defaults = ldsc_cfg.get("n_override") or {}
        n_col_defaults = ldsc_cfg.get("n_col_override") or {}
        headers = {
            t: gwas_header(str(ROOT / gwas_files[t]), (ROOT / gwas_files[t]).stat().st_mtime)
            for t in munge_traits
        }
        n_table = pd.DataFrame(
            [
                {
                    "Trait": t,
                    "Detected N column": detect_n_column(headers[t]) or "— none —",
                    "N column override": str(n_col_defaults.get(t) or ""),
                    "Manual N": float(n_defaults[t]) if n_defaults.get(t) not in (None, "") else None,
                }
                for t in munge_traits
            ],
            columns=["Trait", "Detected N column", "N column override", "Manual N"],
        )
        edited = st.data_editor(
            n_table,
            hide_index=True,
            width="stretch",
            disabled=["Trait", "Detected N column"],
            column_config={
                "N column override": st.column_config.TextColumn(help="Header column holding the per-SNP sample size"),
                "Manual N": st.column_config.NumberColumn(min_value=1, step=1, format="%d", help="Constant N (--N)"),
            },
            key="munge_n_editor_" + "_".join(munge_traits),
        )

        resolved_flags = []
        for row in edited.to_dict("records"):
            trait = row["Trait"]
            manual_n = None if pd.isna(row["Manual N"]) else float(row["Manual N"])
            n_col = str(row["N column override"] or "").strip()
            try:
                flags = munge_column_flags(headers[trait], n_override=manual_n, n_col_override=n_col)
                resolved_flags.append(f"{trait}: {' '.join(flags)}")
            except ValueError as exc:
                n_problems[trait] = str(exc)
            n_override_out[trait] = manual_n if manual_n is not None else ""
            n_col_override_out[trait] = n_col
        if resolved_flags:
            st.code("\n".join(resolved_flags), language="text")
        for trait, message in n_problems.items():
            st.error(f"**{trait}**: {message}")

    munge_blocked = (
        run_active
        or not munge_traits
        or bool(n_problems)
        or not merge_ok
        or not munge_script.is_file()
        or (engine == "slurm" and not has_sbatch)
    )
    if st.button("🧹 Run Munging", type="primary", disabled=munge_blocked, key="munge_run"):
        overrides = {
            "merge_alleles": merge_alleles,
            "n_override": n_override_out,
            "n_col_override": n_col_override_out,
        }
        try:
            start_run(
                "munge",
                munge_traits,
                engine,
                int(cores),
                extra_config=[f"ldsc={json.dumps(overrides)}"],
                forcerun="munge_sumstats" if munged_traits.intersection(munge_traits) else None,
            )
            st.rerun()
        except (OSError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            st.error(f"Could not start munging: {detail}")
    if run_active:
        st.caption("Another run is in progress; munging can start once it finishes.")

    munge_active = run_active and current_run["target"] == "munge"
    st.fragment(munge_monitor, run_every=3 if munge_active else None)()

    st.divider()
    st.subheader("📦 Munged sumstats inventory")
    inventory = munged_inventory(munged_signature())
    if inventory.empty:
        st.info("No munged sumstats in `results/munged/` yet.")
    else:
        st.dataframe(
            inventory,
            hide_index=True,
            width="stretch",
            column_config={
                "SNPs with Z": st.column_config.NumberColumn(format="%d"),
                "Mean chi²": st.column_config.NumberColumn(format="%.3f"),
                "λGC": st.column_config.NumberColumn(format="%.3f"),
            },
        )
        col_trait, col_kind, col_toggle = st.columns([2, 3, 1])
        preview_trait = col_trait.selectbox("Preview trait", inventory["Trait"], key="munge_preview_trait")
        preview_kind = col_kind.radio(
            "File", [".sumstats.gz (first 5 lines)", ".log (full LDSC log)"], horizontal=True, key="munge_preview_kind"
        )
        if col_toggle.toggle("Show preview", key="munge_preview_on"):
            if preview_kind.startswith(".sumstats"):
                st.code(file_head(MUNGED_DIR / f"{preview_trait}.sumstats.gz"), language="text")
            else:
                log_path = MUNGED_DIR / f"{preview_trait}.log"
                if log_path.exists():
                    st.code(log_path.read_text(encoding="utf-8", errors="replace"), language="log", height=400)
                else:
                    st.warning(f"No LDSC log found for {preview_trait}.")

with tab_rg:
    ldsc_cfg = config.get("ldsc") or {}
    rg_script_rel = ldsc_cfg.get("ldsc_script", "../ldsc/ldsc.py")
    ld_ref_dirs = {ldsc_cfg.get("ld_ref_chr", "data/reference/LDscores/"), ldsc_cfg.get("w_ld_chr", "data/reference/LDscores/")}
    rg_script_ok = (ROOT / rg_script_rel).is_file()
    ld_ref_ok = all(any((ROOT / d).glob("*.l2.ldscore.gz")) for d in ld_ref_dirs)
    fdr_alpha = float(config.get("fdr_alpha", 0.05))
    st.caption(
        "Pairwise genetic correlation with LDSC `ldsc.py --rg` on munged summary statistics. Each pair runs "
        "once (traits in alphabetical order) as `results/ldsc_rg/{trait1}__{trait2}.log` and is reused across "
        "modes. Runs with the execution engine selected in the sidebar."
    )
    if not rg_script_ok:
        st.error(f"LDSC script not found: `{rg_script_rel}` (config `ldsc.ldsc_script`).")
    if not ld_ref_ok:
        st.error(
            f"No LD scores (`*.l2.ldscore.gz`) in {', '.join(f'`{d}`' for d in sorted(ld_ref_dirs))} "
            "(config `ldsc.ld_ref_chr` / `ldsc.w_ld_chr`)."
        )

    rg_inventory = munged_inventory(munged_signature())
    if rg_inventory.empty:
        st.info("No munged sumstats in `results/munged/` yet. Munge traits in the **🧹 Munge Sumstats** tab first.")
    else:
        rg_inventory["Raw GWAS"] = rg_inventory["Trait"].isin(gwas_files)
        rg_choices = rg_inventory.loc[rg_inventory["Raw GWAS"], "Trait"].tolist()
        mean_chi2 = dict(zip(rg_inventory["Trait"], rg_inventory["Mean chi²"]))

        def rg_label(trait: str) -> str:
            chi2 = mean_chi2.get(trait)
            return trait if chi2 is None or pd.isna(chi2) else f"{trait} · mean χ² {chi2:.3f}"

        def keep_valid(key: str, options: list) -> None:
            """Drop remembered selections that are no longer offered (e.g. a deleted munged file)."""
            if key in st.session_state:
                st.session_state[key] = [t for t in st.session_state[key] if t in options]

        with st.expander("🩺 Munged sumstats QC", expanded=True):
            st.dataframe(
                rg_inventory[["Trait", "SNPs with Z", "Mean chi²", "λGC", "Raw GWAS"]],
                hide_index=True,
                width="stretch",
                column_config={
                    "SNPs with Z": st.column_config.NumberColumn("SNPs retained", format="%d"),
                    "Mean chi²": st.column_config.NumberColumn("Mean χ²", format="%.3f"),
                    "λGC": st.column_config.NumberColumn(format="%.3f"),
                    "Raw GWAS": st.column_config.CheckboxColumn(
                        help="Traits whose raw GWAS file is no longer in the raw data folder cannot be selected."
                    ),
                },
            )

        rg_mode = st.radio("Analysis Mode", RG_MODES, horizontal=True, key="rg_mode")
        hub_mode = rg_mode == RG_MODES[1]
        if hub_mode:
            rg_target = st.selectbox("Primary / Target trait", rg_choices, key="rg_target", format_func=rg_label)
            panel_choices = [t for t in rg_choices if t != rg_target]
            keep_valid("rg_panel", panel_choices)
            rg_traits = st.multiselect(
                "Comparison panel", panel_choices, key="rg_panel", format_func=rg_label, placeholder="Choose traits"
            )
        else:
            rg_target = None
            keep_valid("rg_matrix_traits", rg_choices)
            col_all, col_clear, _ = st.columns([1, 1, 4])
            col_all.button(
                "Select all", key="rg_select_all", width="stretch",
                on_click=lambda: st.session_state.update(rg_matrix_traits=list(rg_choices)),
            )
            col_clear.button(
                "Clear all", key="rg_clear_all", width="stretch",
                on_click=lambda: st.session_state.update(rg_matrix_traits=[]),
            )
            rg_traits = st.multiselect(
                "Traits (N ≥ 2)", rg_choices, key="rg_matrix_traits", format_func=rg_label, placeholder="Choose traits"
            )

        selected_pairs = rg_pairs(rg_traits, rg_target)
        rg_results = rg_view(load_rg_results(rg_signature()), selected_pairs)
        n_computed = len(rg_results)
        m1, m2, m3 = st.columns(3)
        m1.metric("Regressions (K)" if hub_mode else "Pairwise combinations N(N−1)/2", len(selected_pairs))
        m2.metric("Already computed", n_computed)
        m3.metric("To run", len(selected_pairs) - n_computed)

        rerun_pairs = st.checkbox("Re-run already computed pairs", key="rg_rerun")
        rg_blocked = (
            run_active
            or not selected_pairs
            or not rg_script_ok
            or not ld_ref_ok
            or (engine == "slurm" and not has_sbatch)
            or (n_computed == len(selected_pairs) and not rerun_pairs)
        )
        if st.button("📊 Run Genetic Correlation", type="primary", disabled=rg_blocked, key="rg_run"):
            rg_run_cfg = {"traits": sorted(rg_traits), "target": rg_target or ""}
            try:
                start_run(
                    "rg",
                    sorted(set(rg_traits) | ({rg_target} if rg_target else set())),
                    engine,
                    int(cores),
                    extra_config=[f"rg={json.dumps(rg_run_cfg)}"],
                    forcerun="ldsc_rg" if rerun_pairs else None,
                )
                st.rerun()
            except (OSError, subprocess.CalledProcessError) as exc:
                detail = getattr(exc, "stderr", "") or str(exc)
                st.error(f"Could not start the genetic correlation run: {detail}")
        if run_active:
            st.caption("Another run is in progress; genetic correlation can start once it finishes.")
        elif selected_pairs and n_computed == len(selected_pairs) and not rerun_pairs:
            st.caption("Every pair of this selection is already computed; tick **Re-run** to recompute them.")

    rg_active = run_active and current_run["target"] == "rg"
    st.fragment(rg_monitor, run_every=3 if rg_active else None)()

    if not rg_inventory.empty:
        st.divider()
        st.subheader("📈 Results")
        if rg_results.empty:
            st.info("No computed pairs for this selection yet.")
        else:
            n_estimated = int(rg_results["rg"].notna().sum())
            st.caption(
                f"{n_computed} of {len(selected_pairs)} pairs computed · BH FDR across the {n_estimated} pairs "
                f"shown with an estimated rg · `*` = FDR < {fdr_alpha:g}. Pairs where LDSC could not estimate rg "
                "(e.g. h2 out of bounds) are listed with a note."
            )
            if n_estimated == 0:
                st.warning("No pair of this selection has an estimated rg; see the Note column.")
            elif hub_mode:
                st.plotly_chart(rg_forest(rg_results, rg_target, fdr_alpha), key="rg_forest")
            else:
                st.plotly_chart(rg_heatmap(rg_results, rg_traits, fdr_alpha), key="rg_heatmap")

            rg_table = rg_results.sort_values(["P-value", "Trait 1", "Trait 2"], na_position="last")
            st.dataframe(
                rg_table[RG_TABLE_COLUMNS + ["Note"]],
                hide_index=True,
                width="stretch",
                column_config={
                    "rg": st.column_config.NumberColumn(format="%.3f"),
                    "SE": st.column_config.NumberColumn(format="%.3f"),
                    "Z": st.column_config.NumberColumn(format="%.2f"),
                    "P-value": st.column_config.NumberColumn(format="%.2e"),
                    "FDR": st.column_config.NumberColumn("FDR (BH)", format="%.2e"),
                },
            )
            st.download_button(
                "⬇️ Export table (TSV, incl. h2 & intercepts)",
                rg_table[RG_RESULT_COLUMNS].to_csv(sep="\t", index=False, na_rep="NA").encode("utf-8"),
                file_name="genetic_correlation_hub.tsv" if hub_mode else "genetic_correlation_matrix.tsv",
                mime="text/tab-separated-values",
                key="rg_download",
            )
