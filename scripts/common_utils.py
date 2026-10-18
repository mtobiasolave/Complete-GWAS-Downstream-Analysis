"""Shared utility functions and constants for the GWAS, LDSC & Methylation pipeline.

Consolidates repetitive operations across pipeline scripts:
- Chromosome string normalization ('chr1' vs '1')
- GWAS header detection and streaming O(1)-memory SNP filtering
- Multiple-testing corrections (Benjamini-Hochberg FDR and Storey's q-value)
- Gene symbol set parsing and concordance helpers
- Resilient HTTP session creation and atomic JSON cache I/O
- LDSC genetic correlation log parsing
"""

import gzip
import json
import os
import re
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple, Union

import numpy as np
import pandas as pd

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


# =====================================================================
# 1. Canonical Column & Tissue Constants
# =====================================================================

HEADER_ALIASES: Dict[str, List[str]] = {
    "SNP": ["RSID", "SNP", "SNPID", "MARKERNAME", "ID", "RS_NUMBER", "VARIANT_ID"],
    "P": ["P", "PVAL", "P_VALUE", "PVALUE", "P.VALUE", "P-VALUE", "P-VAL", "P_BOLT_LMM", "P_WALD"],
    "EFFECT": ["BETA", "OR", "B", "Z", "EFFECT", "EFF", "LOG_ODDS", "LOGOR", "ESTIMATE", "ZSCORE"],
    "BETA": ["BETA", "B", "EFFECT", "EFF", "ESTIMATE", "LOGOR", "LOG_ODDS"],
    "OR": ["OR", "ODDS_RATIO", "ODDS-RATIO", "ODDSRATIO"],
    "CHR": ["CHR", "CHROM", "#CHROM", "CHROMOSOME"],
    "BP": ["BP", "POS", "POSITION", "BASE_PAIR_LOCATION"],
    "A1": ["A1", "ALLELE1", "EFFECT_ALLELE", "ALT", "INC_ALLELE", "EA"],
    "A2": ["A2", "ALLELE2", "OTHER_ALLELE", "REF", "DEC_ALLELE", "NEA", "OA"],
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

INVALID_TOKENS: Set[str] = {"", "NA", "NAN", "NONE", "NULL", "."}


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
# 2. File & String Normalization Helpers
# =====================================================================

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


def clean_rsid(val: object) -> str:
    """Normalize an rsID string (strip alleles like rs123:A:G -> rs123, lowercase)."""
    if pd.isna(val):
        return ""
    s = str(val).strip()
    if not s or s.upper() in INVALID_TOKENS:
        return ""
    if ":" in s:
        s = s.split(":", 1)[0].strip()
    return s.lower()


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
    import io

    header_line, matched_lines, _ = extract_matching_gwas_rows(gwas_path, target_snps)
    if not matched_lines:
        return pd.DataFrame(columns=header_line.strip().split())
    buffer = io.StringIO(header_line + "\n" + "\n".join(matched_lines))
    return pd.read_csv(buffer, sep=r"\s+", engine="c", low_memory=False)



# =====================================================================
# 4. Multiple-Testing Correction (BH FDR & Storey's q-value)
# =====================================================================

def _bh_fallback(p: np.ndarray) -> np.ndarray:
    """Pure NumPy Benjamini-Hochberg step-up adjusted p-values."""
    m = len(p)
    if m == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    ranked_p = p[order]
    ranks = np.arange(1, m + 1, dtype=float)
    adj = ranked_p * m / ranks
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)
    out = np.empty(m, dtype=float)
    out[order] = adj
    return out


def storey_qvalues(
    pvals: Sequence[float], lambda_val: float = 0.5, return_pi0: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, float]]:
    """Compute Storey's (2002) q-values estimating pi0 at `lambda_val` (default 0.5).

    Preserves NaN positions and clips q-values to [0, 1].
    """
    p = np.asarray(pvals, dtype=float)
    q_out = np.full(p.shape, np.nan, dtype=float)
    valid_mask = ~np.isnan(p) & (p >= 0.0) & (p <= 1.0)
    p_valid = p[valid_mask]
    m = len(p_valid)
    if m == 0:
        return (q_out, 1.0) if return_pi0 else q_out

    pi0 = np.sum(p_valid > lambda_val) / (m * (1.0 - lambda_val))
    pi0 = float(np.clip(pi0, 0.0, 1.0))
    if pi0 <= 0.0:
        pi0 = 1.0

    sort_idx = np.argsort(p_valid)
    p_sorted = p_valid[sort_idx]
    ranks = np.arange(1, m + 1, dtype=float)

    q_sorted = (pi0 * m * p_sorted) / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    q_valid = np.empty(m, dtype=float)
    q_valid[sort_idx] = q_sorted
    q_out[valid_mask] = q_valid
    return (q_out, pi0) if return_pi0 else q_out


def compute_fdr(pvals: Sequence[float], method: str = "bh") -> np.ndarray:
    """Compute multiple-testing adjusted p-values ('bh'/'fdr_bh' or 'storey'/'qvalue')."""
    method_norm = method.strip().lower()
    if method_norm in ("storey", "qvalue"):
        return storey_qvalues(pvals)

    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid_mask = ~np.isnan(p) & (p >= 0.0) & (p <= 1.0)
    p_valid = p[valid_mask]
    if len(p_valid) == 0:
        return out

    if HAS_STATSMODELS:
        _, p_adj, _, _ = multipletests(p_valid, alpha=0.05, method="fdr_bh")
    else:
        p_adj = _bh_fallback(p_valid)

    out[valid_mask] = p_adj
    return out


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


def find_common_genes(genes_str_1: object, genes_str_2: object) -> str:
    """Return sorted comma-separated intersection of two gene lists (or 'None').

    Handles tissue suffixes such as 'MAPT(Brain_Cortex)' automatically.
    """
    if pd.isna(genes_str_1) or pd.isna(genes_str_2):
        return "None"
    s1 = str(genes_str_1).strip()
    s2 = str(genes_str_2).strip()
    if s1 in ("None", "Error", "", "NA") or s2 in ("None", "Error", "", "NA"):
        return "None"

    set1 = {
        g.split("(")[0].strip().upper(): g.split("(")[0].strip()
        for g in s1.split(",")
        if g.strip() and g.strip().upper() not in INVALID_TOKENS
    }
    set2 = {
        g.split("(")[0].strip().upper()
        for g in s2.split(",")
        if g.strip() and g.strip().upper() not in INVALID_TOKENS
    }
    common_keys = sorted(set(set1.keys()) & set2)
    return ", ".join(set1[k] for k in common_keys) if common_keys else "None"


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
# 7. LDSC Genetic Correlation Log Parser
# =====================================================================

def clean_ldsc_trait_name(filepath_str: str) -> str:
    """Extract clean trait name from LDSC munged filepath (e.g. 'results/munged/AD.sumstats.gz' -> 'AD')."""
    base = os.path.basename(str(filepath_str).strip())
    return base.replace(".sumstats.gz", "").replace(".sumstats", "")


def parse_ldsc_rg_logs(log_files: Iterable[str]) -> pd.DataFrame:
    """Parse 'Summary of Genetic Correlation Results' tables from LDSC .log files."""
    dfs: List[pd.DataFrame] = []

    for log_file in log_files:
        if not os.path.exists(log_file):
            continue
        with open(log_file, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()

        start_idx = None
        for idx, line in enumerate(lines):
            if "Summary of Genetic Correlation Results" in line:
                start_idx = idx + 1
                break

        if start_idx is None or start_idx >= len(lines):
            continue

        header = lines[start_idx].strip().split()
        data_rows: List[List[str]] = []
        for line in lines[start_idx + 1 :]:
            stripped = line.strip()
            if not stripped or stripped.startswith("Analysis finished") or stripped.startswith("Total time"):
                break
            cols = stripped.split()
            if len(cols) == len(header):
                data_rows.append(cols)

        if data_rows:
            dfs.append(pd.DataFrame(data_rows, columns=header))

    if not dfs:
        return pd.DataFrame()

    df = pd.concat(dfs, ignore_index=True)
    for col in ("rg", "se", "z", "p"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "p1" in df.columns and "p2" in df.columns:
        df["trait1"] = df["p1"].map(clean_ldsc_trait_name)
        df["trait2"] = df["p2"].map(clean_ldsc_trait_name)

    return df
