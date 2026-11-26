#!/usr/bin/env python3
"""
Shared pieces of Step B2 (gene annotation of BH-significant RNA m6A SNPs), used by both backends
(rna_genes_ucsc.py: UCSC hg19 / GENCODE V19; rna_genes_ensembl.py: Ensembl GRCh37 VEP).

1. Hardened HTTP client (ApiClient): a requests session that ignores proxy environment variables
   (trust_env = False: a Windows proxy broke TLS) and has no hidden urllib3 retries; explicit
   (connect, read) timeouts; optional request pacing; exponential backoff with full jitter that
   honors Retry-After on 429/500/502/503/504, connection errors and timeouts; a fail-fast server
   check that falls back to the environment proxy only when the direct connection fails (e.g. an
   HPC node that requires one) and logs it.
2. Output schema (identical for both backends) and the record -> row mapping.
3. One annotation cache per backend. Records carry their source in `release`; records of the
   other backend found in a cache file (e.g. the single cache used before both backends existed)
   are moved to that backend's cache in the same directory, so the two databases never mix or
   overwrite each other's annotations.
"""

import html
import json
import os
import random
import re
import threading
import time
from email.utils import parsedate_to_datetime

import requests
from requests.adapters import HTTPAdapter

from common_utils import load_json_cache

LOG = "[annotate_rna_genes]"

SOURCES = ("ucsc", "ensembl")
DEFAULT_SOURCE = "ucsc"
# Default cache file name of each backend (in results/.cache/)
CACHE_NAMES = {
    "ucsc": "ucsc_hg19_genes_cache.json",
    "ensembl": "ensembl_grch37_genes_cache.json",
}
UCSC_RELEASE_PREFIX = "UCSC"     # `release` of UCSC records: "UCSC hg19 GENCODE V19 + <dbSNP track>"

BACKOFF_BASE = 2.0
BACKOFF_CAP = 120.0
RETRYABLE_STATUS = {429, 500, 502, 503, 504}
CONNECTION_ERRORS = (
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
    requests.exceptions.ChunkedEncodingError,
)

PRIMARY_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
STATUS_LABELS = {
    "intergenic": "intergenic",
    "not_in_dbsnp": "not in dbSNP",          # UCSC: rsID in no hg19 dbSNP track
    "not_in_ensembl": "not in Ensembl",      # Ensembl: rsID unknown to VEP
    "no_primary_locus": "no primary-assembly locus",
}
OUTPUT_COLUMNS = [
    "SNP",
    "GENES",
    "ENSEMBL_GENE_IDS",
    "STATUS",
    "CHR",
    "POS",
    "CONSEQUENCES",
    "TRANSCRIPT_BIOTYPES",
    "ENSEMBL_RELEASE",
]


# =====================================================================
# 1. Network Layer
# =====================================================================

def ascii_safe(text) -> str:
    """Log-safe text: Windows consoles (cp1252) crash on characters such as U+2019 in error pages."""
    return str(text).encode("ascii", "replace").decode("ascii")


def describe_response(resp: requests.Response) -> str:
    """'HTTP <status>' plus a one-line message (JSON error, HTML page title, or text)."""
    # Error pages often omit the charset (requests then assumes ISO-8859-1): decode as UTF-8
    text = (resp.content or b"").decode("utf-8", "replace")
    message = ""
    if "json" in resp.headers.get("Content-Type", ""):
        try:
            body = resp.json()
            message = body.get("error", "") if isinstance(body, dict) else ""
        except ValueError:
            pass
    if not message:
        title = re.search(r"(?is)<title>(.*?)</title>", text)
        message = html.unescape(title.group(1)) if title else text
    message = ascii_safe(re.sub(r"\s+", " ", str(message)).strip()[:200])
    return f"HTTP {resp.status_code}" + (f": {message}" if message else "")


def retry_delay(attempt: int, resp: requests.Response = None) -> float:
    """Server Retry-After (seconds or HTTP date) if given, else capped exponential backoff with
    full jitter."""
    retry_after = resp.headers.get("Retry-After") if resp is not None else None
    if retry_after:
        try:
            return min(BACKOFF_CAP, max(0.0, float(retry_after)))
        except ValueError:
            try:
                return min(BACKOFF_CAP, max(0.0, parsedate_to_datetime(retry_after).timestamp() - time.time()))
            except (TypeError, ValueError):
                pass
    return random.uniform(0, min(BACKOFF_CAP, BACKOFF_BASE * 2 ** attempt))


def create_session(pool_size: int = 2) -> requests.Session:
    """Session that ignores proxy environment variables and has no hidden urllib3 retries
    (urllib3's own retries hide the HTTP status behind 'too many 500 error responses')."""
    session = requests.Session()
    session.trust_env = False
    session.proxies = {}
    adapter = HTTPAdapter(pool_connections=pool_size, pool_maxsize=pool_size, max_retries=0)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


class ApiClient:
    """Retrying, optionally paced REST client shared by the UCSC and Ensembl backends.
    Thread-safe for concurrent requests (pacing and counters are locked)."""

    name = "REST API"

    def __init__(self, server: str, timeout: tuple, attempts: int, min_interval: float = 0.0,
                 session: requests.Session = None, pool_size: int = 2):
        self.server = server.rstrip("/")
        self.timeout = timeout
        self.attempts = attempts
        self.min_interval = min_interval
        self.session = session or create_session(pool_size)
        self.requests_sent = 0
        self.last_failure_was_connection = False
        self._lock = threading.Lock()
        self._last_start = 0.0

    def _pace(self) -> None:
        with self._lock:
            wait = self._last_start + self.min_interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            self._last_start = time.monotonic()
            self.requests_sent += 1

    def request(self, method: str, url: str, label: str, attempts: int = None, **kwargs):
        """Send a request, retrying retryable statuses, connection errors and timeouts.
        Returns (response or None, last error description)."""
        attempts = attempts or self.attempts
        last_error = "no attempt made"
        for attempt in range(attempts):
            self._pace()
            started = time.time()
            resp = None
            try:
                resp = self.session.request(method, url, timeout=self.timeout, **kwargs)
                self.last_failure_was_connection = False
                if resp.status_code == 200:
                    return resp, ""
                last_error = describe_response(resp)
                retryable = resp.status_code in RETRYABLE_STATUS
            except CONNECTION_ERRORS as exc:
                last_error = ascii_safe(f"{type(exc).__name__}: {exc}")
                self.last_failure_was_connection = True
                retryable = True
            except requests.exceptions.RequestException as exc:
                last_error = ascii_safe(f"{type(exc).__name__}: {exc}")
                self.last_failure_was_connection = True
                retryable = False
            elapsed = time.time() - started
            if not retryable or attempt == attempts - 1:
                print(f"{LOG} {label} failed after {attempt + 1} attempt(s) ({elapsed:.1f}s): {last_error}", flush=True)
                break
            delay = retry_delay(attempt, resp)
            print(f"{LOG} {label} attempt {attempt + 1}/{attempts} ({elapsed:.1f}s): {last_error}; "
                  f"retrying in {delay:.1f}s", flush=True)
            time.sleep(delay)
        return None, last_error

    def get_json(self, path: str, params: dict, label: str, attempts: int = None, headers: dict = None):
        """JSON payload of a GET endpoint, or (None, error) on failure or an API 'error' field."""
        resp, error = self.request("GET", f"{self.server}{path}", label, attempts, params=params, headers=headers)
        if resp is None:
            return None, error
        try:
            payload = resp.json()
        except ValueError as exc:
            return None, f"invalid JSON ({exc})"
        if isinstance(payload, dict) and payload.get("error"):
            return None, ascii_safe(f"API error: {payload['error']}")
        return payload, ""

    def check(self, path: str, params: dict, is_ok, headers: dict = None,
              allow_proxy_fallback: bool = True) -> None:
        """Fail fast when the server does not answer `path` as expected. If the direct connection
        fails and a proxy is set in the environment, retry through that proxy (and keep it)."""

        def probe(label):
            payload, error = self.get_json(path, params, label, attempts=3, headers=headers)
            if payload is not None and not is_ok(payload):
                return False, ascii_safe(f"unexpected response: {str(payload)[:200]}")
            return payload is not None, error

        ok, error = probe("server check")
        if not ok and allow_proxy_fallback and self.last_failure_was_connection:
            env_proxies = requests.utils.get_environ_proxies(self.server)
            if env_proxies:
                print(f"{LOG} Direct connection to {self.server} failed ({error}); "
                      f"retrying through the proxy set in the environment.", flush=True)
                self.session.proxies = dict(env_proxies)
                ok, error = probe("server check (proxy)")
                if not ok:
                    self.session.proxies = {}
        if not ok:
            raise RuntimeError(
                f"{LOG} {self.name} {self.server} is unavailable ({error}). "
                "No SNPs were queried; re-run when the server is reachable."
            )


# =====================================================================
# 2. Output Rows
# =====================================================================

def annotation_row(snp: str, record: dict) -> dict:
    """Output row (OUTPUT_COLUMNS) of one cached annotation record, from either backend."""
    return {
        "SNP": snp,
        "GENES": ", ".join(record["genes"]) or STATUS_LABELS.get(record["status"], record["status"]),
        "ENSEMBL_GENE_IDS": ", ".join(record["gene_ids"]) or "NA",
        "STATUS": record["status"],
        "CHR": record["chr"] or "NA",
        "POS": record["pos"] or "NA",
        "CONSEQUENCES": ", ".join(record["consequences"]) or "NA",
        "TRANSCRIPT_BIOTYPES": ", ".join(record["biotypes"]) or "NA",
        "ENSEMBL_RELEASE": record["release"],
    }


# =====================================================================
# 3. Per-Backend Annotation Caches
# =====================================================================

def record_source(record: dict) -> str:
    """Backend that produced a cached record: 'ucsc' (release 'UCSC hg19 GENCODE V19 + ...') or
    'ensembl' (release = Ensembl release number(s), or 'NA')."""
    return "ucsc" if str(record.get("release", "")).startswith(UCSC_RELEASE_PREFIX) else "ensembl"


def default_cache_path(cache_dir, source: str) -> str:
    return os.path.join(str(cache_dir), CACHE_NAMES[source])


def write_json_atomic(path: str, data: dict) -> None:
    """Replace `path` with `data` (whole file, via a PID-tagged temp file)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, path)


def load_backend_cache(cache_path: str, source: str) -> dict:
    """Cached records of `source` from its cache file.

    The other backend's cache is assumed to sit in the same directory under its default name
    (CACHE_NAMES). Records found in the wrong file are moved to their backend's file first
    (a record already in its own file wins), so switching backends never mixes or overwrites
    annotations. Only records of `source` are returned.
    """
    other = next(s for s in SOURCES if s != source)
    paths = {source: os.path.abspath(cache_path),
             other: os.path.abspath(default_cache_path(os.path.dirname(os.path.abspath(cache_path)), other))}
    if paths[source] == paths[other]:                       # both backends configured to one file
        return {k: v for k, v in load_json_cache(cache_path).items() if record_source(v) == source}

    data = {s: load_json_cache(p) for s, p in paths.items()}
    final = {s: {k: v for k, v in data[s].items() if record_source(v) == s} for s in paths}
    moved = {s: 0 for s in paths}
    for holder, records in data.items():
        for rsid, record in records.items():
            owner = record_source(record)
            if owner != holder and rsid not in final[owner]:
                final[owner][rsid] = record
                moved[owner] += 1
    if any(moved.values()) or any(len(final[s]) != len(data[s]) for s in paths):
        # Write the files that gain records before the ones that lose them: an interruption
        # can leave a duplicate (filtered on the next load), never a lost record
        for s in sorted(paths, key=lambda s: -moved[s]):
            if final[s] or os.path.exists(paths[s]):
                write_json_atomic(paths[s], final[s])
        for s, n in moved.items():
            if n:
                print(f"{LOG} Moved {n:,} cached {s} annotations into {paths[s]} "
                      f"(caches are kept separate per backend).", flush=True)
    return final[source]
