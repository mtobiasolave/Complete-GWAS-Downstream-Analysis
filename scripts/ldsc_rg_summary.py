#!/usr/bin/env python3
"""
Genetic Correlation (rg): Parse LDSC --rg Logs and Apply BH FDR Across Trait Pairs.

Each pair is run once as results/ldsc_rg/{trait1}__{trait2}.log with the two traits in
alphabetical order (rg is symmetric), so a pair computed in one analysis mode is reused
by the other. Shared by the Snakefile (pair selection, rg_summary rule) and app.py.

Output columns:
    Trait 1, Trait 2, rg, SE, Z, P-value, FDR,
    h2 (Trait 1), h2 SE (Trait 1), h2 (Trait 2), h2 SE (Trait 2), gcov_int, gcov_int SE, Note

Usage:
    python scripts/ldsc_rg_summary.py --logs results/ldsc_rg/A__B.log results/ldsc_rg/A__C.log \
        --alpha 0.05 --output results/ldsc_rg/rg_summary.tsv
"""

import argparse
import itertools
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from common_utils import bh_fdr

PAIR_SEP = "__"
RG_TABLE_COLUMNS: List[str] = ["Trait 1", "Trait 2", "rg", "SE", "Z", "P-value", "FDR"]
RG_EXTRA_COLUMNS: List[str] = [
    "h2 (Trait 1)",
    "h2 SE (Trait 1)",
    "h2 (Trait 2)",
    "h2 SE (Trait 2)",
    "gcov_int",
    "gcov_int SE",
    "Note",
]
RG_RESULT_COLUMNS: List[str] = RG_TABLE_COLUMNS + RG_EXTRA_COLUMNS

SUMMARY_HEADER = "Summary of Genetic Correlation Results"
# Columns after p1/p2 in the LDSC summary row, in LDSC's order
SUMMARY_FIELDS = ["rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se"]
H2_RE = re.compile(r"Total Observed scale h2: (\S+) \((\S+)\)")
RG_NOTE_RE = re.compile(r"Genetic Correlation: nan \(nan\) \(([^)]*)\)")


def canonical_pair(trait_a: str, trait_b: str) -> Tuple[str, str]:
    """Alphabetical (trait1, trait2) order used for every rg output file."""
    return tuple(sorted((trait_a, trait_b)))


def pair_stem(trait_a: str, trait_b: str) -> str:
    """File stem of a pair's LDSC output, e.g. 'Anxiety__SUD2023'."""
    return PAIR_SEP.join(canonical_pair(trait_a, trait_b))


def rg_pairs(traits: Iterable[str], target: Optional[str] = None) -> List[Tuple[str, str]]:
    """All-vs-all pairs of `traits`, or `target` vs each of `traits` (hub & spokes)."""
    traits = sorted({t for t in traits if t})
    if target:
        return sorted({canonical_pair(target, t) for t in traits if t != target})
    return [canonical_pair(a, b) for a, b in itertools.combinations(traits, 2)]


def _to_float(token: str) -> float:
    try:
        return float(token)
    except (TypeError, ValueError):
        return np.nan


def parse_rg_log(path) -> Dict[str, object]:
    """Parse one LDSC --rg log for a single pair into a RG_RESULT_COLUMNS record.

    Trait names come from the file stem (trait1__trait2), not the sumstats paths LDSC prints.
    Values LDSC reports as NA (e.g. h2 out of bounds) are NaN with the reason in Note.
    """
    path = Path(path)
    trait1, _, trait2 = path.stem.partition(PAIR_SEP)
    record: Dict[str, object] = {col: np.nan for col in RG_RESULT_COLUMNS}
    record.update({"Trait 1": trait1, "Trait 2": trait2, "Note": ""})

    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == SUMMARY_HEADER)
        tokens = lines[start + 2].split()
    except (StopIteration, IndexError):
        record["Note"] = "incomplete log (no LDSC summary)"
        return record
    if len(tokens) < len(SUMMARY_FIELDS) + 2:
        record["Note"] = "unparseable LDSC summary row"
        return record

    values = dict(zip(SUMMARY_FIELDS, map(_to_float, tokens[-len(SUMMARY_FIELDS):])))
    record.update({
        "rg": values["rg"],
        "SE": values["se"],
        "Z": values["z"],
        "P-value": values["p"],
        "gcov_int": values["gcov_int"],
        "gcov_int SE": values["gcov_int_se"],
    })
    h2 = H2_RE.findall(text)
    if len(h2) >= 2:
        record["h2 (Trait 1)"], record["h2 SE (Trait 1)"] = map(_to_float, h2[0])
        record["h2 (Trait 2)"], record["h2 SE (Trait 2)"] = map(_to_float, h2[1])
    if np.isnan(values["rg"]):
        reason = RG_NOTE_RE.search(text)
        record["Note"] = " ".join(reason.group(1).split()) if reason else "rg not estimated"
    return record


def add_fdr(df: pd.DataFrame) -> pd.DataFrame:
    """BH FDR across the given pairs (the testing family); pairs without a valid P stay NaN."""
    df = df.copy()
    df["FDR"] = bh_fdr(df["P-value"].to_numpy(dtype=float)) if not df.empty else []
    return df


def summarize(log_paths: Sequence) -> pd.DataFrame:
    rows = [parse_rg_log(p) for p in log_paths]
    df = pd.DataFrame(rows, columns=RG_RESULT_COLUMNS)
    df = add_fdr(df).sort_values(["P-value", "Trait 1", "Trait 2"], na_position="last")
    return df[RG_RESULT_COLUMNS].reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize LDSC --rg logs with BH FDR across pairs.")
    parser.add_argument("--logs", nargs="*", default=[], help="results/ldsc_rg/{trait1}__{trait2}.log files")
    parser.add_argument("--alpha", type=float, default=0.05, help="BH FDR level (default: 0.05)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV (results/ldsc_rg/rg_summary.tsv)")
    args = parser.parse_args()

    summary = summarize(args.logs)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    n_valid = int(summary["rg"].notna().sum())
    n_sig = int((summary["FDR"] <= args.alpha).sum())
    print(
        f"[ldsc_rg_summary] {len(summary):,} pairs ({n_valid:,} with a valid rg) | "
        f"{n_sig:,} significant at BH FDR <= {args.alpha} -> {args.output}"
    )


if __name__ == "__main__":
    main()
