#!/usr/bin/env python3
"""
Step B2 backend "ensembl": genes of BH-significant RNA m6A SNPs from the Ensembl VEP REST API
(GRCh37 server by default, matching the hg19 coordinates of the DNA track). Called by
annotate_rna_genes.py (--source ensembl).

1. Queries `POST /vep/homo_sapiens/id` in batches of rsIDs (2 concurrent requests).
2. Keeps the genes whose transcripts the variant physically falls inside: any transcript
   consequence other than upstream/downstream_gene_variant (VEP reports those up to 5 kb away).
   Only loci on primary chromosomes (1-22, X, Y, MT) are used.
3. Successfully queried SNPs are cached after every batch; failed requests are never cached.

Resilience for the GRCh37 archive server (unmirrored; transient 500/502/503/504, timeouts, rate
limits), on top of rna_genes_common.ApiClient (proxy-independent session, timeouts, backoff with
jitter and Retry-After): the server is pinged first (an outage fails fast with the HTTP status),
batches are small (25 rsIDs) with long read timeouts, a batch that still fails is retried one
SNP at a time after a health check, and the run stops querying if the server stops answering.

STATUS: genic | intergenic | not_in_ensembl (rsID unknown to Ensembl) |
        no_primary_locus (maps only to patches/alternative haplotypes).
"""

import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

from common_utils import save_json_cache
from rna_genes_common import LOG, PRIMARY_CHROMS, ApiClient

SOURCE_NAME = "Ensembl VEP"
DEFAULT_SERVER = "https://grch37.rest.ensembl.org"
DEFAULT_WORKERS = 2
# (connect, read) timeouts in seconds: VEP on the archive server can take minutes per batch
TIMEOUT = (15, 300)
BATCH_ATTEMPTS = 6
SINGLE_ATTEMPTS = 3
JSON_HEADERS = {"Accept": "application/json"}
# Transcript consequences assigned to nearby (not overlapping) genes
NEARBY_TERMS = {"upstream_gene_variant", "downstream_gene_variant"}


class VepClient(ApiClient):
    """Client for the Ensembl REST API (concurrent requests allowed, no pacing)."""

    name = "Ensembl REST server"

    def __init__(self, server: str = DEFAULT_SERVER, session=None, timeout: tuple = TIMEOUT,
                 pool_size: int = DEFAULT_WORKERS + 2):
        super().__init__(server, timeout, BATCH_ATTEMPTS, min_interval=0.0, session=session, pool_size=pool_size)

    def check_server(self, allow_proxy_fallback: bool = True) -> None:
        """Fail fast with the HTTP status when the server does not answer /info/ping."""
        self.check("/info/ping", None, lambda payload: isinstance(payload, dict) and payload.get("ping") == 1,
                   headers=JSON_HEADERS, allow_proxy_fallback=allow_proxy_fallback)

    def fetch_release(self) -> str:
        """Ensembl release served by the server (recorded for reproducibility), or 'NA'."""
        payload, error = self.get_json("/info/data", None, "release lookup", attempts=3, headers=JSON_HEADERS)
        if isinstance(payload, dict):
            return ", ".join(str(r) for r in payload.get("releases", [])) or "NA"
        print(f"{LOG} Could not read the Ensembl release ({error or 'unexpected response'}); recording 'NA'.")
        return "NA"

    def post_vep(self, rsids: list, attempts: int) -> tuple:
        """POST rsIDs to VEP. Returns ({input rsID: [VEP entries]} or None, error description)."""
        label = f"VEP {rsids[0]}" if len(rsids) == 1 else f"VEP batch of {len(rsids)} ({rsids[0]}..{rsids[-1]})"
        resp, error = self.request(
            "POST", f"{self.server}/vep/homo_sapiens/id", label, attempts,
            headers={"Content-Type": "application/json", **JSON_HEADERS}, json={"ids": rsids},
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


def fetch_vep_batch(client: VepClient, rsids: list, abort: threading.Event) -> tuple:
    """Query one batch; if it still fails, retry its SNPs one at a time.
    Returns ({rsID: [VEP entries]} for resolved SNPs, [failed rsIDs])."""
    if abort.is_set():
        return {}, list(rsids)
    entries, error = client.post_vep(rsids, BATCH_ATTEMPTS)
    if entries is not None:
        return {rsid: entries.get(rsid, []) for rsid in rsids}, []
    if len(rsids) == 1:
        return {}, list(rsids)

    # A failing batch on a healthy server is usually one problematic rsID: isolate it
    try:
        client.check_server(allow_proxy_fallback=False)
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
        single, _ = client.post_vep([rsid], SINGLE_ATTEMPTS)
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


def annotate(snps: list, cache: dict, cache_path: str, server: str = DEFAULT_SERVER, batch_size: int = 25,
             max_workers: int = None, client: VepClient = None) -> list:
    """Annotate uncached `snps` into `cache` (+ its file, flushed every batch). Returns failed rsIDs."""
    max_workers = max(1, max_workers or DEFAULT_WORKERS)
    client = client or VepClient(server, pool_size=max_workers + 2)
    client.check_server()
    release = client.fetch_release()
    batches = [snps[i: i + batch_size] for i in range(0, len(snps), batch_size)]
    print(f"{LOG} Querying Ensembl VEP at {client.server} (release {release}): {len(batches)} batches of "
          f"<= {batch_size} rsIDs, {max_workers} concurrent request(s)", flush=True)
    n_new, failed_snps = 0, []
    abort = threading.Event()  # set when the server stops answering /info/ping mid-run

    with ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as executor:
        futures = [executor.submit(fetch_vep_batch, client, batch, abort) for batch in batches]
        for done, future in enumerate(as_completed(futures), start=1):
            resolved, failed = future.result()
            failed_snps.extend(failed)
            records = {rsid: summarize_vep_entries(entries, release) for rsid, entries in resolved.items()}
            if records:
                # Flush after every batch (main thread only) so progress survives a failed run
                cache.update(records)
                save_json_cache(cache_path, records)
                n_new += len(records)
            print(f"{LOG} {done}/{len(batches)} batches done | resolved {n_new:,} | failed {len(failed_snps):,}",
                  flush=True)

    print(f"{LOG} Cached {n_new:,} newly annotated SNPs at {cache_path} ({client.requests_sent:,} requests)")
    if abort.is_set():
        print(f"{LOG} The Ensembl server stopped responding mid-run.", flush=True)
    return sorted(set(failed_snps))
