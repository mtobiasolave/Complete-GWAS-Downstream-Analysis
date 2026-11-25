#!/usr/bin/env python3
"""
Step B2: Positional Gene Annotation of BH-Significant RNA m6A SNPs via Ensembl VEP.

1. Collects every SNP that is BH-significant in at least one trait (Step B1 outputs).
2. Queries the Ensembl VEP REST API (`POST /vep/homo_sapiens/id`, GRCh37 server by
   default to match the hg19 coordinates of the DNA track) in batches of rsIDs.
3. Keeps the genes whose transcripts the variant physically falls inside: any transcript
   consequence other than upstream/downstream_gene_variant (VEP reports those up to 5 kb
   away). Only loci on primary chromosomes (1-22, X, Y, MT) are used.
4. Uses a persistent atomic JSON cache (`results/.cache/ensembl_grch37_genes_cache.json`),
   flushed after every batch so partial progress survives a failed run. Only successfully
   queried SNPs are cached; failed requests make the job fail so the remaining SNPs are
   retried on the next run.

Resilience for the GRCh37 archive server (unmirrored; transient 500/502/503/504, timeouts,
rate limits): the server is pinged first (an outage fails fast with the HTTP status), batches
are small (25 rsIDs) with long read timeouts, retryable failures (429, 500, 502, 503, 504,
connection errors, timeouts) back off exponentially with full jitter and honor Retry-After,
a batch that still fails is retried one SNP at a time, and every failure logs its HTTP
status and message.

Output columns:
    SNP, GENES, ENSEMBL_GENE_IDS, STATUS, CHR, POS, CONSEQUENCES, TRANSCRIPT_BIOTYPES,
    ENSEMBL_RELEASE

STATUS: genic | intergenic | not_in_ensembl (rsID unknown to Ensembl) |
        no_primary_locus (maps only to patches/alternative haplotypes).
GENES holds the gene symbols (Ensembl gene ID when no symbol), or the status label.

Usage:
    python scripts/annotate_rna_genes.py --inputs results/rna/*_overlap_bh.tsv \
        --output results/rna/rna_sig_genes.tsv
"""

import argparse
import html
import os
import random
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from email.utils import parsedate_to_datetime

import pandas as pd
import requests
from requests.adapters import HTTPAdapter

from common_utils import PROJECT_ROOT, load_json_cache, save_json_cache
from pleiotropy import is_true

DEFAULT_SERVER = "https://grch37.rest.ensembl.org"
DEFAULT_CACHE_FILE = str(PROJECT_ROOT / "results" / ".cache" / "ensembl_grch37_genes_cache.json")
DEFAULT_BATCH_SIZE = 25
DEFAULT_WORKERS = 2

# (connect, read) timeouts in seconds: VEP on the archive server can take minutes per batch
TIMEOUT = (15, 300)
RETRYABLE_STATUS = {429, 500, 502, 503, 504}
RETRYABLE_EXCEPTIONS = (
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
    requests.exceptions.ChunkedEncodingError,
)
BATCH_ATTEMPTS = 6
SINGLE_ATTEMPTS = 3
BACKOFF_BASE = 2.0
BACKOFF_CAP = 120.0
LOG = "[annotate_rna_genes]"

PRIMARY_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
# Transcript consequences assigned to nearby (not overlapping) genes
NEARBY_TERMS = {"upstream_gene_variant", "downstream_gene_variant"}
STATUS_LABELS = {
    "intergenic": "intergenic",
    "not_in_ensembl": "not in Ensembl",
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


def significant_snps(input_paths: list) -> list:
    """Unique SNPs BH-significant in at least one Step B1 overlap table."""
    snps = set()
    for path in input_paths:
        df = pd.read_csv(path, sep="\t", usecols=["SNP", "BH_SIGNIFICANT"], dtype=str)
        snps.update(df.loc[is_true(df["BH_SIGNIFICANT"]), "SNP"].str.strip())
    return sorted(snps)


def create_vep_session(pool_size: int) -> requests.Session:
    """Pooled session WITHOUT urllib3 retries, so every HTTP status reaches the retry loop
    below (urllib3's own retries hide it behind 'too many 500 error responses')."""
    session = requests.Session()
    adapter = HTTPAdapter(pool_connections=pool_size, pool_maxsize=pool_size, max_retries=0)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def ascii_safe(text: str) -> str:
    """Log-safe text: Windows consoles (cp1252) crash on characters such as U+2019 in EBI pages."""
    return text.encode("ascii", "replace").decode("ascii")


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
    message = ascii_safe(re.sub(r"\s+", " ", message).strip()[:200])
    return f"HTTP {resp.status_code}" + (f": {message}" if message else "")


def retry_delay(attempt: int, resp: requests.Response = None) -> float:
    """Server Retry-After (seconds or HTTP date) if given, else capped exponential backoff
    with full jitter."""
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


def request_with_retries(session, method: str, url: str, label: str, attempts: int, **kwargs):
    """Send a request, retrying retryable statuses and connection errors.
    Returns (response or None, last error description)."""
    last_error = "no attempt made"
    for attempt in range(attempts):
        started = time.time()
        resp = None
        try:
            resp = session.request(method, url, timeout=TIMEOUT, **kwargs)
            if resp.status_code == 200:
                return resp, ""
            last_error = describe_response(resp)
            retryable = resp.status_code in RETRYABLE_STATUS
        except RETRYABLE_EXCEPTIONS as exc:
            last_error = ascii_safe(f"{type(exc).__name__}: {exc}")
            retryable = True
        except requests.exceptions.RequestException as exc:
            last_error = ascii_safe(f"{type(exc).__name__}: {exc}")
            retryable = False
        elapsed = time.time() - started
        if not retryable or attempt == attempts - 1:
            print(f"{LOG} {label} failed after {attempt + 1} attempt(s) ({elapsed:.1f}s): {last_error}", flush=True)
            break
        delay = retry_delay(attempt, resp)
        print(
            f"{LOG} {label} attempt {attempt + 1}/{attempts} ({elapsed:.1f}s): {last_error}; "
            f"retrying in {delay:.1f}s",
            flush=True,
        )
        time.sleep(delay)
    return None, last_error


def check_server(server: str, session) -> None:
    """Fail fast with the HTTP status when the REST server does not answer /info/ping."""
    resp, error = request_with_retries(
        session, "GET", f"{server}/info/ping", "ping", attempts=3, headers={"Accept": "application/json"}
    )
    try:
        ok = resp is not None and resp.json().get("ping") == 1
    except ValueError:
        ok, error = False, ascii_safe(f"unexpected ping response: {resp.text[:200]!r}")
    if not ok:
        raise RuntimeError(
            f"{LOG} Ensembl REST server {server} is unavailable ({error}). "
            "No SNPs were queried; re-run when the server is back."
        )


def fetch_release(server: str, session) -> str:
    """Ensembl release served by `server` (recorded for reproducibility)."""
    resp, error = request_with_retries(
        session, "GET", f"{server}/info/data", "release lookup", attempts=3, headers={"Accept": "application/json"}
    )
    if resp is not None:
        try:
            return ", ".join(str(r) for r in resp.json().get("releases", [])) or "NA"
        except ValueError as exc:
            error = f"invalid JSON: {exc}"
    print(f"{LOG} Could not read the Ensembl release ({error}); recording 'NA'.")
    return "NA"


def post_vep(rsids: list, server: str, session, attempts: int) -> tuple:
    """POST rsIDs to VEP. Returns ({input rsID: [VEP entries]} or None, error description)."""
    label = f"VEP {rsids[0]}" if len(rsids) == 1 else f"VEP batch of {len(rsids)} ({rsids[0]}..{rsids[-1]})"
    resp, error = request_with_retries(
        session, "POST", f"{server}/vep/homo_sapiens/id", label, attempts,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        json={"ids": rsids},
    )
    if resp is None:
        return None, error
    try:
        entries = {}
        for entry in resp.json():
            entries.setdefault(entry.get("input") or entry.get("id"), []).append(entry)
        return entries, ""
    except (ValueError, AttributeError, TypeError) as exc:
        return None, f"invalid VEP JSON ({type(exc).__name__}: {exc})"


def fetch_vep_batch(rsids: list, server: str, session, abort: threading.Event) -> tuple:
    """Query one batch; if it still fails, retry its SNPs one at a time.
    Returns ({rsID: [VEP entries]} for resolved SNPs, [failed rsIDs])."""
    if abort.is_set():
        return {}, list(rsids)
    entries, error = post_vep(rsids, server, session, BATCH_ATTEMPTS)
    if entries is not None:
        return {rsid: entries.get(rsid, []) for rsid in rsids}, []
    if len(rsids) == 1:
        return {}, list(rsids)

    # A failing batch on a healthy server is usually one problematic rsID: isolate it
    try:
        check_server(server, session)
    except RuntimeError as exc:
        print(exc, flush=True)
        abort.set()
        return {}, list(rsids)
    print(f"{LOG} Batch {rsids[0]}..{rsids[-1]} failed ({error}); retrying its {len(rsids)} SNPs one at a time.", flush=True)
    resolved, failed = {}, []
    for rsid in rsids:
        if abort.is_set():
            failed.append(rsid)
            continue
        single, _ = post_vep([rsid], server, session, SINGLE_ATTEMPTS)
        if single is None:
            failed.append(rsid)
        else:
            resolved[rsid] = single.get(rsid, [])
    return resolved, failed


def summarize_vep_entries(entries: list, release: str) -> dict:
    """Collapse the VEP entries of one rsID into its overlapping genes and status."""
    loci = [e for e in entries if str(e.get("seq_region_name", "")) in PRIMARY_CHROMS]
    genes, consequences, biotypes = {}, set(), set()
    for entry in loci:
        for tc in entry.get("transcript_consequences", []):
            terms = set(tc.get("consequence_terms") or []) - NEARBY_TERMS
            gene_id = tc.get("gene_id")
            if not terms or not gene_id:
                continue
            genes[gene_id] = tc.get("gene_symbol") or gene_id
            consequences |= terms
            if tc.get("biotype"):
                biotypes.add(tc["biotype"])

    if not entries:
        status = "not_in_ensembl"
    elif not loci:
        status = "no_primary_locus"
    else:
        status = "genic" if genes else "intergenic"
    return {
        "status": status,
        "genes": sorted(set(genes.values())),
        "gene_ids": sorted(genes),
        "chr": ", ".join(str(e["seq_region_name"]) for e in loci),
        "pos": ", ".join(str(e.get("start")) for e in loci),
        "consequences": sorted(consequences),
        "biotypes": sorted(biotypes),
        "release": release,
    }


def annotation_row(snp: str, record: dict) -> dict:
    return {
        "SNP": snp,
        "GENES": ", ".join(record["genes"]) or STATUS_LABELS[record["status"]],
        "ENSEMBL_GENE_IDS": ", ".join(record["gene_ids"]) or "NA",
        "STATUS": record["status"],
        "CHR": record["chr"] or "NA",
        "POS": record["pos"] or "NA",
        "CONSEQUENCES": ", ".join(record["consequences"]) or "NA",
        "TRANSCRIPT_BIOTYPES": ", ".join(record["biotypes"]) or "NA",
        "ENSEMBL_RELEASE": record["release"],
    }


def run_annotate_rna_genes(
    input_paths: list,
    output_path: str,
    server: str = DEFAULT_SERVER,
    cache_path: str = DEFAULT_CACHE_FILE,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_workers: int = DEFAULT_WORKERS,
) -> pd.DataFrame:
    server = server.rstrip("/")
    snps = significant_snps(input_paths)
    print(f"{LOG} {len(input_paths)} traits: {len(snps):,} unique BH-significant SNPs to annotate.")

    cache = load_json_cache(cache_path)
    uncached = [s for s in snps if s not in cache]
    print(f"{LOG} Cache hits: {len(snps) - len(uncached):,} | Uncached SNPs to query: {len(uncached):,}")

    if uncached:
        session = create_vep_session(pool_size=max_workers + 2)
        check_server(server, session)
        release = fetch_release(server, session)
        batches = [uncached[i : i + batch_size] for i in range(0, len(uncached), batch_size)]
        print(
            f"{LOG} Querying Ensembl VEP at {server} (release {release}): {len(batches)} batches of "
            f"<= {batch_size} rsIDs, {max_workers} concurrent request(s)",
            flush=True,
        )
        n_new, failed_snps = 0, []
        abort = threading.Event()  # set when the server stops answering /info/ping mid-run

        with ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as executor:
            futures = {executor.submit(fetch_vep_batch, b, server, session, abort): b for b in batches}
            for done, future in enumerate(as_completed(futures), start=1):
                resolved, failed = future.result()
                failed_snps.extend(failed)
                records = {rsid: summarize_vep_entries(entries, release) for rsid, entries in resolved.items()}
                if records:
                    # Flush after every batch (main thread only) so progress survives a failed run
                    cache.update(records)
                    save_json_cache(cache_path, records)
                    n_new += len(records)
                print(
                    f"{LOG} {done}/{len(batches)} batches done | resolved {n_new:,} | failed {len(failed_snps):,}",
                    flush=True,
                )

        # Only fully retrieved SNPs are cached; failed SNPs are retried on the next run
        print(f"{LOG} Cached {n_new:,} newly queried SNPs at {cache_path}")

        if failed_snps:
            reason = " The server stopped responding mid-run." if abort.is_set() else ""
            raise RuntimeError(
                f"{LOG} Ensembl VEP requests failed for {len(failed_snps):,} SNPs "
                f"(e.g. {sorted(failed_snps)[:5]}).{reason} Successful results were cached; "
                f"re-run the pipeline to retry the remaining SNPs."
            )

    out = pd.DataFrame([annotation_row(s, cache[s]) for s in snps], columns=OUTPUT_COLUMNS)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    counts = out["STATUS"].value_counts().to_dict()
    print(f"[annotate_rna_genes] {len(out):,} SNPs {counts} -> {output_path}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate BH-significant RNA m6A SNPs with the Ensembl genes they fall inside."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Step B1 results/rna/{trait}_overlap_bh.tsv files")
    parser.add_argument("--output", "-o", required=True, help="Output results/rna/rna_sig_genes.tsv")
    parser.add_argument("--server", default=DEFAULT_SERVER, help=f"Ensembl REST server (default: {DEFAULT_SERVER})")
    parser.add_argument("--cache-file", default=DEFAULT_CACHE_FILE, help="Persistent JSON cache path")
    parser.add_argument(
        "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help=f"rsIDs per VEP POST request (default: {DEFAULT_BATCH_SIZE})",
    )
    parser.add_argument(
        "--workers", type=int, default=DEFAULT_WORKERS,
        help=f"Max concurrent VEP requests (default: {DEFAULT_WORKERS})",
    )
    args = parser.parse_args()

    run_annotate_rna_genes(
        args.inputs,
        args.output,
        server=args.server,
        cache_path=args.cache_file,
        batch_size=args.batch_size,
        max_workers=args.workers,
    )


if __name__ == "__main__":
    main()
