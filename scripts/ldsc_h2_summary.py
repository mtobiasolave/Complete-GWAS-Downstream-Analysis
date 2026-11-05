#!/usr/bin/env python3
"""
SNP Heritability (h2): Build LDSC --h2 Flags, Parse --h2 Logs, Summarize Across Traits.

Each trait runs once as results/ldsc_h2/{trait}.log. Shared by the Snakefile (flag building
and validation, h2_summary rule) and app.py (validation, settings comparison, log parsing).

Liability scale: LDSC prints either "Total Observed scale h2" or, when --samp-prev and
--pop-prev are both given, "Total Liability scale h2" (never both). The observed-scale value
is recovered exactly by dividing by LDSC's own conversion factor (h2_obs_to_liab with h2 = 1),
using the prevalences recorded in the log's Call block.

Z = h2 / SE (identical on both scales) and P-value = 1 - Phi(Z), one-sided for h2 > 0.
LDSC rounds h2 and SE to 4 decimals in its log, so Z and P inherit that rounding.

Output columns:
    Trait, Scale, h2 (observed), SE (observed), h2 (liability), SE (liability), Z, P-value,
    Pop prev (K), Samp prev (P), Intercept, Intercept SE, Lambda GC, Mean chi2, Ratio, Ratio SE,
    SNPs, Note

Usage:
    python scripts/ldsc_h2_summary.py --logs results/ldsc_h2/A.log results/ldsc_h2/B.log \
        --output results/ldsc_h2/h2_summary.tsv
"""

import argparse
import math
import os
import re
from pathlib import Path
from statistics import NormalDist
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

H2_TABLE_COLUMNS: List[str] = ["Trait", "Scale", "h2 (observed)", "SE (observed)", "h2 (liability)", "SE (liability)", "Z", "P-value"]
H2_EXTRA_COLUMNS: List[str] = [
    "Pop prev (K)",
    "Samp prev (P)",
    "Intercept",
    "Intercept SE",
    "Lambda GC",
    "Mean chi2",
    "Ratio",
    "Ratio SE",
    "SNPs",
    "Note",
]
H2_RESULT_COLUMNS: List[str] = H2_TABLE_COLUMNS + H2_EXTRA_COLUMNS

TOTAL_RE = re.compile(r"Total (Observed|Liability) scale h2: (\S+) \((\S+)\)")
LAMBDA_RE = re.compile(r"^Lambda GC: (\S+)", re.M)
MEAN_CHI2_RE = re.compile(r"^Mean Chi\^2: (\S+)", re.M)
INTERCEPT_RE = re.compile(r"^Intercept: (\S+) \((\S+)\)", re.M)
INTERCEPT_FIXED_RE = re.compile(r"^Intercept: constrained to ([\d.eE+-]+?)\.?\s*$", re.M)
RATIO_RE = re.compile(r"^Ratio: (\S+) \((\S+)\)", re.M)
RATIO_NOTE_RE = re.compile(r"^(Ratio < 0 \(usually indicates GC correction\)|Ratio: NA \(mean chi\^2 < 1\))", re.M)
SNPS_RE = re.compile(r"After merging with regression SNP LD, (\d+) SNPs remain")
CALL_FLAG_RE = re.compile(r"^--(samp-prev|pop-prev|intercept-h2) (\S+)", re.M)
MUNGE_N_RE = re.compile(r"^--(N-col|N) (\S+)", re.M)


def _to_float(token) -> float:
    try:
        return float(token)
    except (TypeError, ValueError):
        return np.nan


def _unset(value) -> bool:
    return value is None or (isinstance(value, str) and not value.strip()) or (
        isinstance(value, float) and np.isnan(value)
    )


# =====================================================================
# 1. Settings: Prevalences (Liability Scale) & Intercept Constraint
# =====================================================================

def validate_prevalence(pop_prev, samp_prev) -> Optional[Tuple[float, float]]:
    """Return (K, P) for the liability scale, None for the observed scale.

    Both prevalences or neither (co-dependent), each strictly between 0 and 1.
    """
    if _unset(pop_prev) and _unset(samp_prev):
        return None
    if _unset(pop_prev) or _unset(samp_prev):
        raise ValueError("Set both or neither of the population (--pop-prev) and sample (--samp-prev) prevalences.")
    pop, samp = float(pop_prev), float(samp_prev)
    for name, value in (("Population prevalence (--pop-prev)", pop), ("Sample prevalence (--samp-prev)", samp)):
        if not 0.0 < value < 1.0:
            raise ValueError(f"{name} must be strictly between 0 and 1, got {value:g}.")
    return pop, samp


def requested_h2_settings(trait: str, h2_cfg: dict) -> Dict[str, Optional[float]]:
    """Validated settings of one trait from config `h2` ({prevalence: {trait: {...}}, intercept_h2})."""
    prevalence = ((h2_cfg.get("prevalence") or {}).get(trait)) or {}
    liability = validate_prevalence(prevalence.get("pop_prev"), prevalence.get("samp_prev"))
    intercept = h2_cfg.get("intercept_h2")
    if not _unset(intercept):
        intercept = float(intercept)
        if not intercept > 0:
            raise ValueError(f"Fixed intercept (--intercept-h2) must be positive, got {intercept:g}.")
    else:
        intercept = None
    return {
        "pop_prev": liability[0] if liability else None,
        "samp_prev": liability[1] if liability else None,
        "intercept_h2": intercept,
    }


def h2_flags(trait: str, h2_cfg: dict) -> List[str]:
    """Extra ldsc.py --h2 flags: prevalences only for the liability scale, intercept if fixed."""
    settings = requested_h2_settings(trait, h2_cfg)
    flags: List[str] = []
    if settings["pop_prev"] is not None:
        flags += ["--samp-prev", repr(settings["samp_prev"]), "--pop-prev", repr(settings["pop_prev"])]
    if settings["intercept_h2"] is not None:
        flags += ["--intercept-h2", repr(settings["intercept_h2"])]
    return flags


def logged_h2_settings(log_text: str) -> Dict[str, Optional[float]]:
    """Settings an existing LDSC --h2 log was computed with, from its Call block."""
    flags = {name.replace("-", "_"): float(value) for name, value in CALL_FLAG_RE.findall(log_text)}
    return {key: flags.get(key) for key in ("pop_prev", "samp_prev", "intercept_h2")}


def liability_factor(samp_prev: float, pop_prev: float) -> float:
    """LDSC's observed-to-liability factor: K^2 (1-K)^2 / (P (1-P) phi(Phi^-1(1-K))^2)."""
    normal = NormalDist()
    threshold = normal.inv_cdf(1.0 - pop_prev)
    return pop_prev ** 2 * (1 - pop_prev) ** 2 / (samp_prev * (1 - samp_prev) * normal.pdf(threshold) ** 2)


def munge_n_source(munge_log: Path) -> Tuple[str, str]:
    """How a trait's munged N was set, from its munge_sumstats.py log Call.

    Returns (kind, detail) with kind in: effective, half_effective, total, manual, unknown.
    """
    text = munge_log.read_text(encoding="utf-8", errors="replace") if munge_log.exists() else ""
    match = MUNGE_N_RE.search(text)
    if not match:
        return "unknown", "no N flag found in the munge log"
    flag, value = match.groups()
    if flag == "N":
        return "manual", f"manual N = {_to_float(value):,.0f}"
    if "half" in value.lower():
        return "half_effective", f"column {value}"
    if "eff" in value.lower():
        return "effective", f"column {value}"
    return "total", f"column {value}"


# =====================================================================
# 2. LDSC --h2 Log Parsing
# =====================================================================

def parse_h2_log(path) -> Dict[str, object]:
    """Parse one LDSC --h2 log into a H2_RESULT_COLUMNS record (trait name from the file stem)."""
    path = Path(path)
    record: Dict[str, object] = {col: np.nan for col in H2_RESULT_COLUMNS}
    record.update({"Trait": path.stem, "Scale": "", "Note": ""})

    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    total = TOTAL_RE.search(text)
    if not total:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        errors = [line for line in lines if "Error" in line or "error" in line]
        record["Note"] = f"LDSC error: {errors[-1]}" if errors else "incomplete log (no LDSC summary)"
        return record

    scale, h2, se = total.group(1), _to_float(total.group(2)), _to_float(total.group(3))
    settings = logged_h2_settings(text)
    record["Scale"] = scale
    record["Pop prev (K)"] = settings["pop_prev"] if settings["pop_prev"] is not None else np.nan
    record["Samp prev (P)"] = settings["samp_prev"] if settings["samp_prev"] is not None else np.nan
    if scale == "Liability":
        factor = liability_factor(settings["samp_prev"], settings["pop_prev"])
        record.update({
            "h2 (liability)": h2,
            "SE (liability)": se,
            "h2 (observed)": h2 / factor,
            "SE (observed)": se / factor,
        })
    else:
        record.update({"h2 (observed)": h2, "SE (observed)": se})

    if se > 0:
        record["Z"] = h2 / se
        # One-sided upper tail via erfc: NormalDist().cdf(-Z) cancels to 0 for Z > ~8
        record["P-value"] = 0.5 * math.erfc(record["Z"] / math.sqrt(2.0))

    for column, pattern in (("Lambda GC", LAMBDA_RE), ("Mean chi2", MEAN_CHI2_RE)):
        match = pattern.search(text)
        record[column] = _to_float(match.group(1)) if match else np.nan
    free, fixed = INTERCEPT_RE.search(text), INTERCEPT_FIXED_RE.search(text)
    notes = []
    if free:
        record["Intercept"], record["Intercept SE"] = _to_float(free.group(1)), _to_float(free.group(2))
    elif fixed:
        record["Intercept"] = _to_float(fixed.group(1))
        notes.append(f"intercept constrained to {fixed.group(1)}")
    ratio, ratio_note = RATIO_RE.search(text), RATIO_NOTE_RE.search(text)
    if ratio:
        record["Ratio"], record["Ratio SE"] = _to_float(ratio.group(1)), _to_float(ratio.group(2))
    elif ratio_note:
        notes.append(ratio_note.group(1))
    snps = SNPS_RE.search(text)
    record["SNPs"] = int(snps.group(1)) if snps else np.nan
    record["Note"] = "; ".join(notes)
    return record


def summarize(log_paths: Sequence) -> pd.DataFrame:
    df = pd.DataFrame([parse_h2_log(p) for p in log_paths], columns=H2_RESULT_COLUMNS)
    return df.sort_values("Trait").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize LDSC --h2 logs (observed & liability scale).")
    parser.add_argument("--logs", nargs="*", default=[], help="results/ldsc_h2/{trait}.log files")
    parser.add_argument("--output", "-o", required=True, help="Output TSV (results/ldsc_h2/h2_summary.tsv)")
    args = parser.parse_args()

    summary = summarize(args.logs)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    summary.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    n_estimated = int(summary["h2 (observed)"].notna().sum())
    n_liability = int((summary["Scale"] == "Liability").sum())
    print(
        f"[ldsc_h2_summary] {len(summary):,} traits ({n_estimated:,} with an h2 estimate, "
        f"{n_liability:,} on the liability scale) -> {args.output}"
    )


if __name__ == "__main__":
    main()
