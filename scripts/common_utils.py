"""Shared utility functions and constants for the DNA & RNA Methylation x GWAS pipeline.

Consolidates repetitive operations across pipeline scripts:
- Project-root path resolution and GWAS summary statistics discovery
- Chromosome string normalization ('chr1' vs '1')
- GWAS header detection and streaming O(1)-memory SNP filtering
- Benjamini-Hochberg (BH) FDR correction of candidate-space GWAS overlaps
- Gene symbol set formatting helpers
- Resilient HTTP session creation and atomic JSON cache I/O
- LDSC munge_sumstats.py column-flag mapping (SNP, alleles, P, sample size)
"""

import gzip
import io
import json
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd

# Repository root (scripts/ lives one level below it); all default paths resolve from here
PROJECT_ROOT = Path(__file__).resolve().parents[1]


# =====================================================================
# 1. Canonical Column & Tissue Constants
# =====================================================================

HEADER_ALIASES: Dict[str, List[str]] = {
    "SNP": ["RSID", "SNP", "SNPID", "MARKERNAME", "ID", "RS_NUMBER", "VARIANT_ID"],
    "P": ["P", "PVAL", "P_VALUE", "PVALUE", "P.VALUE", "P-VALUE", "P-VAL", "P_BOLT_LMM", "P_WALD"],
    "BETA": ["BETA", "B", "EFFECT", "EFF", "ESTIMATE", "LOGOR", "LOG_ODDS"],
    "OR": ["OR", "ODDS_RATIO", "ODDS-RATIO", "ODDSRATIO"],
}

DEFAULT_BRAIN_TISSUES: List[str] = [
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
]

TISSUE_TYPO_FIXES: Dict[str, str] = {
    "Brain_Cerebellar_Hemispher": "Brain_Cerebellar_Hemisphere",
    "Brain_Amydgala": "Brain_Amygdala",
}

INVALID_TOKENS: Set[str] = {"", "NA", "NAN", "NONE", "NULL", ".", "<NA>", "NONE_WITHIN_5KB"}

# Recognized GWAS file suffixes, longest first so "_gwas.txt.gz" wins over ".txt.gz"
GWAS_SUFFIXES: Tuple[str, ...] = (
    "_gwas.txt.gz",
    "_gwas.tsv.gz",
    "_gwas.txt",
    "_gwas.tsv",
    ".txt.gz",
    ".tsv.gz",
    ".txt",
    ".tsv",
)


def normalize_tissues(tissues_arg: Optional[Sequence[str]] = None) -> List[str]:
    """Return a clean list of GTEx v10 brain tissue IDs, fixing common typos."""
    raw = list(tissues_arg) if tissues_arg else list(DEFAULT_BRAIN_TISSUES)
    cleaned: List[str] = []
    for t in raw:
        t_fixed = TISSUE_TYPO_FIXES.get(t.strip(), t.strip())
        if t_fixed and t_fixed not in cleaned:
            cleaned.append(t_fixed)
    return cleaned


# =====================================================================
# 2. GWAS Discovery, File & String Normalization Helpers
# =====================================================================

def discover_gwas(raw_dir: Union[str, Path] = "data/raw_data") -> Dict[str, str]:
    """Map each trait name to its GWAS file path (POSIX-style, as given) in `raw_dir`."""
    raw_dir = Path(raw_dir)
    discovered: Dict[str, str] = {}
    if not raw_dir.is_dir():
        return discovered
    for fpath in sorted(raw_dir.iterdir()):
        if fpath.name.startswith(".") or not fpath.is_file():
            continue
        for suf in GWAS_SUFFIXES:
            if fpath.name.endswith(suf):
                trait_name = fpath.name[: -len(suf)]
                if trait_name:
                    discovered[trait_name] = fpath.as_posix()
                break
    return discovered


def is_gzipped(filepath: str) -> bool:
    """Check if a file is gzip-compressed by inspecting its first two magic bytes."""
    with open(filepath, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def standardize_chr(chr_input: Union[pd.Series, str, int]) -> Union[pd.Series, str]:
    """Standardize chromosome identifiers to 'chr1'..'chr22', 'chrX', 'chrY', 'chrM'.

    Works on both a pandas Series (fully vectorized) and scalar strings/integers.
    """
    if isinstance(chr_input, pd.Series):
        s = chr_input.astype(str).str.strip()
        s = s.str.replace(r"\.0$", "", regex=True)
        lower = s.str.lower()
        stripped = np.where(lower.str.startswith("chr"), s.str.slice(3), s)
        stripped_s = pd.Series(stripped, index=chr_input.index, dtype="object")
        upper = stripped_s.str.upper()
        norm_body = np.where(
            upper.isin(["X", "Y", "M", "MT"]),
            np.where(upper == "MT", "M", upper),
            stripped_s,
        )
        res = "chr" + pd.Series(norm_body, index=chr_input.index, dtype="object")
        empty_mask = s.isin(["", "nan", "None", "NA", "<NA>"])
        if empty_mask.any():
            res = res.mask(empty_mask, "")
        return res

    val = str(chr_input).strip()
    if val.endswith(".0"):
        val = val[:-2]
    if not val or val.lower() in ("nan", "none", "na"):
        return ""
    if val.lower().startswith("chr"):
        val = val[3:]
    u = val.upper()
    if u == "MT":
        val = "M"
    elif u in ("X", "Y", "M"):
        val = u
    return f"chr{val}"


# =====================================================================
# 3. GWAS Header Detection & Streaming Extraction
# =====================================================================

def detect_column(columns: Iterable[str], col_type: str) -> Optional[str]:
    """Find the actual column name in a DataFrame matching known aliases for `col_type`."""
    col_map = {c.upper(): c for c in columns}
    for candidate in HEADER_ALIASES.get(col_type, []):
        if candidate.upper() in col_map:
            return col_map[candidate.upper()]
    return None


def find_snp_col_index(header_tokens: Sequence[str]) -> int:
    """Locate the 0-based index of the SNP identifier column in a GWAS header."""
    upper_map = {col.strip().upper(): idx for idx, col in enumerate(header_tokens)}
    for candidate in HEADER_ALIASES["SNP"]:
        if candidate in upper_map:
            return upper_map[candidate]
    raise ValueError(
        f"Could not find a valid SNP column in GWAS header: {list(header_tokens)}. "
        f"Expected one of {HEADER_ALIASES['SNP']}."
    )


def extract_matching_gwas_rows(
    gwas_path: str, target_snps: Set[str]
) -> Tuple[str, List[str], Set[str]]:
    """Stream a raw GWAS file (plain or .gz) and extract lines matching `target_snps`.

    Memory footprint is O(matched rows), avoiding loading multi-GB GWAS files into RAM.

    Returns:
        (header_line, matched_lines, matched_snps_set)
    """
    open_fn = gzip.open if is_gzipped(gwas_path) else open
    matched_lines: List[str] = []
    matched_snps: Set[str] = set()

    with open_fn(gwas_path, "rt", encoding="utf-8", errors="replace") as f:
        header_line = f.readline()
        if not header_line:
            raise ValueError(f"GWAS file is empty: {gwas_path}")

        header_cols = header_line.strip().split()
        snp_col_idx = find_snp_col_index(header_cols)
        max_split = snp_col_idx + 1

        if snp_col_idx == 0:
            for line in f:
                tokens = line.split(None, 1)
                if tokens and tokens[0] in target_snps:
                    matched_lines.append(line.rstrip("\r\n"))
                    matched_snps.add(tokens[0])
        else:
            for line in f:
                tokens = line.split(None, max_split)
                if len(tokens) > snp_col_idx and tokens[snp_col_idx] in target_snps:
                    snp_val = tokens[snp_col_idx]
                    matched_lines.append(line.rstrip("\r\n"))
                    matched_snps.add(snp_val)

    return header_line.rstrip("\r\n"), matched_lines, matched_snps


def extract_matching_gwas_df(gwas_path: str, target_snps: Set[str]) -> pd.DataFrame:
    """Stream a raw GWAS file and return a DataFrame of rows matching `target_snps`."""
    header_line, matched_lines, _ = extract_matching_gwas_rows(gwas_path, target_snps)
    if not matched_lines:
        return pd.DataFrame(columns=header_line.strip().split())
    buffer = io.StringIO(header_line + "\n" + "\n".join(matched_lines))
    return pd.read_csv(buffer, sep=r"\s+", engine="c", low_memory=False)


# =====================================================================
# 4. Benjamini-Hochberg FDR Correction of GWAS Overlaps
# =====================================================================

BH_COLUMNS: List[str] = ["SNP", "GWAS_P", "GWAS_BETA", "GWAS_FDR_BH", "BH_SIGNIFICANT"]


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg (1995) step-up adjusted p-values.

    Invalid entries (NaN or outside [0, 1]) are excluded from the testing family
    and returned as NaN in their original positions.
    """
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid_mask = ~np.isnan(p) & (p >= 0.0) & (p <= 1.0)
    p_valid = p[valid_mask]
    m = len(p_valid)
    if m == 0:
        return out

    order = np.argsort(p_valid, kind="mergesort")
    adj = p_valid[order] * m / np.arange(1, m + 1, dtype=float)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj_valid = np.empty(m, dtype=float)
    adj_valid[order] = np.clip(adj, 0.0, 1.0)
    out[valid_mask] = adj_valid
    return out


GWAS_LOOKUP_COLUMNS: List[str] = ["SNP", "GWAS_P", "GWAS_BETA"]


def gwas_lookup(gwas_path: str, candidate_snps: Set[str]) -> pd.DataFrame:
    """Extract a reference SNP set from one GWAS with harmonized P and effect size.

    Duplicate GWAS rows keep the first occurrence. Effect sizes are harmonized to
    BETA, converting OR to ln(OR) when no BETA column is present.

    Returns a DataFrame with GWAS_LOOKUP_COLUMNS, one row per overlapping SNP.
    """
    df_raw = extract_matching_gwas_df(gwas_path, candidate_snps)

    snp_col = detect_column(df_raw.columns, "SNP")
    p_col = detect_column(df_raw.columns, "P")
    beta_col = detect_column(df_raw.columns, "BETA")
    or_col = detect_column(df_raw.columns, "OR")
    if not snp_col or not p_col:
        raise ValueError(
            f"Could not detect SNP ({snp_col}) or P ({p_col}) column in {gwas_path}. "
            f"Columns found: {list(df_raw.columns)}"
        )

    if beta_col:
        beta_vals = pd.to_numeric(df_raw[beta_col], errors="coerce")
    elif or_col:
        or_vals = pd.to_numeric(df_raw[or_col], errors="coerce")
        beta_vals = np.log(or_vals.where(or_vals > 0))
    else:
        beta_vals = pd.Series(np.nan, index=df_raw.index)

    overlap = (
        pd.DataFrame({
            "SNP": df_raw[snp_col].astype(str).str.strip(),
            "GWAS_P": pd.to_numeric(df_raw[p_col], errors="coerce"),
            "GWAS_BETA": beta_vals,
        })
        .drop_duplicates(subset=["SNP"], keep="first")
        .reset_index(drop=True)
    )
    print(
        f"[gwas_lookup] {os.path.basename(gwas_path)}: {len(overlap):,} of {len(candidate_snps):,} "
        f"candidate SNPs found (P col='{p_col}', effect col='{beta_col or or_col}', "
        f"OR->ln(OR)={bool(not beta_col and or_col)})."
    )
    return overlap[GWAS_LOOKUP_COLUMNS]


def gwas_overlap_bh(gwas_path: str, candidate_snps: Set[str], alpha: float = 0.05) -> pd.DataFrame:
    """Overlap a reference SNP set with one GWAS and apply BH FDR within that trait.

    The testing family is every unique candidate SNP found in the GWAS with a valid
    p-value (see gwas_lookup for extraction and effect-size harmonization).

    Returns a DataFrame with BH_COLUMNS, one row per overlapping SNP.
    """
    overlap = gwas_lookup(gwas_path, candidate_snps)
    overlap["GWAS_FDR_BH"] = bh_fdr(overlap["GWAS_P"].to_numpy())
    overlap["BH_SIGNIFICANT"] = overlap["GWAS_FDR_BH"] <= alpha

    n_tested = int(overlap["GWAS_FDR_BH"].notna().sum())
    print(
        f"[gwas_overlap_bh] {os.path.basename(gwas_path)}: {n_tested:,} candidate SNPs tested | "
        f"{int(overlap['BH_SIGNIFICANT'].sum()):,} significant at BH FDR <= {alpha}."
    )
    return overlap[BH_COLUMNS]


# =====================================================================
# 5. Gene Set & Token Formatting Helpers
# =====================================================================

def join_unique_clean(values: Iterable[object]) -> str:
    """Split comma/semicolon-separated values, deduplicate, sort, and join with ', '."""
    tokens: Set[str] = set()
    for val in values:
        if pd.isna(val):
            continue
        for part in re.split(r"[,;]", str(val)):
            t = part.strip()
            if t and t.upper() not in INVALID_TOKENS:
                tokens.add(t)
    if not tokens:
        return "NA"
    return ", ".join(sorted(tokens))


# =====================================================================
# 6. HTTP Session & Atomic JSON Cache Helpers
# =====================================================================

def create_http_session(retries: int = 3, backoff_factor: float = 0.5, pool_size: int = 20):
    """Create a requests.Session configured with connection pooling and exponential backoff retries."""
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    session = requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "POST"]),
    )
    adapter = HTTPAdapter(
        pool_connections=pool_size, pool_maxsize=pool_size, max_retries=retry
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def load_json_cache(cache_path: Optional[str]) -> dict:
    """Load a JSON dictionary cache from disk if it exists, otherwise return {}."""
    if not cache_path or not os.path.exists(cache_path):
        return {}
    try:
        with open(cache_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            return data if isinstance(data, dict) else {}
    except Exception as exc:
        print(f"[WARNING] Could not read cache file {cache_path}: {exc}")
        return {}


def save_json_cache(cache_path: Optional[str], cache_data: dict) -> None:
    """Atomically merge and write a JSON dictionary cache to disk via a PID-tagged temp file."""
    if not cache_path or not cache_data:
        return
    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    tmp_path = f"{cache_path}.{os.getpid()}.tmp"
    try:
        existing = load_json_cache(cache_path)
        existing.update(cache_data)
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(existing, f, indent=2)
        os.replace(tmp_path, cache_path)
    except Exception as exc:
        print(f"[WARNING] Could not save cache file {cache_path}: {exc}")


# =====================================================================
# 7. LDSC munge_sumstats.py Column Mapping
# =====================================================================

MUNGE_ALIASES: Dict[str, List[str]] = {
    "--snp": ["SNP", "snp", "rsID", "rsid"],
    "--a1": ["A1", "a1", "Allele1", "effect_allele", "alt", "Aallele"],
    "--a2": ["A2", "a2", "Allele2", "other_allele", "ref", "Ballele"],
    "--p": ["P", "p", "pvalue", "p_value", "P_value", "P.value"],
}
MUNGE_N_ALIASES: List[str] = ["N", "n", "N_eff", "Neff", "N_EFF", "N_TOTAL", "Weight"]


def read_gwas_header(gwas_path: Union[str, Path]) -> List[str]:
    """Return the whitespace-split header of a plain or gzipped GWAS file."""
    open_fn = gzip.open if is_gzipped(str(gwas_path)) else open
    with open_fn(gwas_path, "rt", encoding="utf-8", errors="replace") as f:
        return f.readline().strip().split()


def detect_n_column(header: Sequence[str]) -> Optional[str]:
    """First sample-size column in `header` matching MUNGE_N_ALIASES, or None."""
    return next((c for c in MUNGE_N_ALIASES if c in header), None)


def munge_column_flags(
    header: Sequence[str], n_override: object = None, n_col_override: object = None
) -> List[str]:
    """Map a GWAS header to munge_sumstats.py column flags.

    Sample size precedence: manual N (--N) > N column override > auto-detected N column.
    Empty strings / None mean "not set". Raises ValueError with a readable message.
    """
    flags: List[str] = []
    for ldsc_arg, possible_names in MUNGE_ALIASES.items():
        match = next((c for c in possible_names if c in header), None)
        if match is None:
            raise ValueError(
                f"No column for '{ldsc_arg}' in header {list(header)}; expected one of {possible_names}."
            )
        flags += [ldsc_arg, match]

    if n_override not in (None, ""):
        n_value = float(n_override)
        if n_value <= 0:
            raise ValueError(f"Manual N must be positive, got {n_override}.")
        flags += ["--N", str(int(n_value)) if n_value.is_integer() else repr(n_value)]
    elif n_col_override not in (None, ""):
        if n_col_override not in header:
            raise ValueError(f"N column '{n_col_override}' is not in header {list(header)}.")
        flags += ["--N-col", str(n_col_override)]
    else:
        n_col = detect_n_column(header)
        if n_col is None:
            raise ValueError(
                f"No sample-size column in header {list(header)} (expected one of {MUNGE_N_ALIASES}). "
                f"Set a manual N (config ldsc.n_override or the app's Munge tab)."
            )
        flags += ["--N-col", n_col]
    return flags
