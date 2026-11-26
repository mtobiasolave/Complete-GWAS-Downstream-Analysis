#!/usr/bin/env python3
"""
Step B2: Positional Gene Annotation of BH-Significant RNA m6A SNPs (UCSC hg19 / GENCODE V19).

1. Collects every SNP that is BH-significant in at least one trait (Step B1 outputs, rsIDs only).
2. Resolves each uncached rsID to its hg19 (= GRCh37) position(s) with the UCSC Genome Browser
   REST API search endpoint (`/search?search=<rsID>&genome=hg19`), using the newest dbSNP track
   that contains the rsID (dbSnp155, then older builds). Only loci on primary chromosomes
   (1-22, X, Y, MT) are used.
3. Fetches the GENCODE V19 transcripts around those SNPs (`/getData/track`, tracks
   wgEncodeGencodeCompV19 + wgEncodeGencodePseudoGeneV19; nearby SNPs share one region query)
   and keeps the genes whose transcripts the variant physically falls inside (exon or intron;
   no upstream/downstream window), all biotypes: the rule of the former Ensembl VEP step.
   Gene IDs and transcript biotypes come from UCSC's GENCODE V19 attributes table
   (wgEncodeGencodeAttrsV19, downloaded once from hgdownload.soe.ucsc.edu).
4. Persistent atomic JSON caches, flushed after every batch and on exit: the annotation cache
   (--cache-file; it also holds the SNPs annotated earlier by Ensembl VEP, which are reused
   as-is) and the resolved rsID positions (ucsc_hg19_rsid_positions.json next to it).
   Failed lookups are never cached; they make the job fail so a re-run retries only them.

Coordinates: output CHR uses Ensembl-style names (1..22, X, Y, MT) and POS is 1-based, as before.
UCSC uses chr-prefixed names and 0-based half-open intervals; conversions are explicit below.

Network: UCSC asks for sequential requests and about one request per second, so requests are
paced, retried with exponential backoff and jitter (honoring Retry-After) on 429/500/502/503/504,
connection errors and timeouts. Proxy environment variables are ignored (session.trust_env =
False: a Windows proxy broke TLS); only if the direct connection fails and proxy variables are
set (e.g. an HPC that requires one) is the proxy from the environment used, and this is logged.

Output columns (unchanged):
    SNP, GENES, ENSEMBL_GENE_IDS, STATUS, CHR, POS, CONSEQUENCES, TRANSCRIPT_BIOTYPES,
    ENSEMBL_RELEASE (annotation source: "UCSC hg19 GENCODE V19 + <dbSNP track>" for new SNPs,
    the Ensembl release for SNPs annotated earlier by VEP)

STATUS: genic | intergenic | not_in_dbsnp (rsID not found in any UCSC hg19 dbSNP track) |
        no_primary_locus (maps only to alternative haplotypes) | not_in_ensembl (legacy VEP records).
CONSEQUENCES are positional Sequence Ontology terms (coding_sequence_variant,
5_prime_UTR_variant, 3_prime_UTR_variant, intron_variant, non_coding_transcript_exon_variant,
non_coding_transcript_variant); allele-dependent terms (e.g. missense_variant) need the
alleles and are not derived.

Usage:
    python scripts/annotate_rna_genes.py --inputs results/rna/*_overlap_bh.tsv \
        --output results/rna/rna_sig_genes.tsv
"""

import argparse
import gzip
import html
import os
import random
import re
import time
from email.utils import parsedate_to_datetime
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
import requests
from requests.adapters import HTTPAdapter

from common_utils import PROJECT_ROOT, load_json_cache, save_json_cache
from pleiotropy import is_true

DEFAULT_SERVER = "https://api.genome.ucsc.edu"
DEFAULT_CACHE_FILE = str(PROJECT_ROOT / "results" / ".cache" / "ensembl_grch37_genes_cache.json")
POSITIONS_CACHE_NAME = "ucsc_hg19_rsid_positions.json"
GENOME = "hg19"
# GENCODE V19 comprehensive transcripts + pseudogenes (pseudogenes are a separate UCSC table)
GENE_TRACKS = ("wgEncodeGencodeCompV19", "wgEncodeGencodePseudoGeneV19")
# Small tracks (16,175 pseudogenes genome-wide) are fetched once per chromosome, not per region
WHOLE_CHROMOSOME_TRACKS = {"wgEncodeGencodePseudoGeneV19"}
ATTRS_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeAttrsV19.txt.gz"
ATTRS_CACHE_NAME = "ucsc_hg19_wgEncodeGencodeAttrsV19.txt.gz"
SOURCE_LABEL = "UCSC hg19 GENCODE V19"
# dbSNP tracks in UCSC hg19 searches: dbSnp155, dbSnp153, snp151, snp150, ... (subsets such as
# dbSnp155Common don't match this pattern); the highest build that contains the rsID wins
SNP_TRACK_RE = re.compile(r"^(?:dbSnp|snp)(\d+)$")
HIGHLIGHT_RE = re.compile(r"^(?:[^:]*\.)?(chr[^:.]+):(\d+)-(\d+)")

DEFAULT_BATCH_SIZE = 25          # SNPs per cache checkpoint
CLUSTER_MAX_GAP = 500_000        # SNPs closer than this share one region query ...
CLUSTER_MAX_SPAN = 5_000_000     # ... as long as the region stays below this size
MIN_INTERVAL = 1.0               # seconds between request starts (UCSC: ~1 request per second)
TIMEOUT = (15, 180)              # (connect, read) seconds
ATTEMPTS = 6
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
    "not_in_dbsnp": "not in dbSNP",
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
LOG = "[annotate_rna_genes]"


# =====================================================================
# 1. Coordinates & Chromosome Names
# =====================================================================

def to_ucsc_chrom(chrom) -> str:
    """Ensembl/bare chromosome name -> UCSC ('1' -> 'chr1', 'MT' -> 'chrM')."""
    body = str(chrom).strip()
    if body.lower().startswith("chr"):
        body = body[3:]
    if body.upper() in ("M", "MT"):
        return "chrM"
    if body.upper() in ("X", "Y"):
        body = body.upper()
    return f"chr{body}"


def from_ucsc_chrom(chrom) -> str:
    """UCSC chromosome name -> Ensembl style ('chr1' -> '1', 'chrM' -> 'MT')."""
    body = str(chrom).strip()
    if body.lower().startswith("chr"):
        body = body[3:]
    if body.upper() in ("M", "MT"):
        return "MT"
    if body.upper() in ("X", "Y"):
        return body.upper()
    return body


def is_primary(chrom) -> bool:
    """True for 1-22, X, Y, MT (in either naming); False for haplotypes/patches/unplaced."""
    return from_ucsc_chrom(chrom) in PRIMARY_CHROMS


def ensembl_to_ucsc(start: int, end: int) -> tuple:
    """1-based closed [start, end] -> UCSC 0-based half-open [start - 1, end)."""
    return int(start) - 1, int(end)


def ucsc_to_ensembl(start: int, end: int) -> tuple:
    """UCSC 0-based half-open [start, end) -> 1-based closed [start + 1, end]."""
    return int(start) + 1, int(end)


def locus_query_interval(start1: int, end1: int) -> tuple:
    """0-based half-open interval of a 1-based locus; insertions (end < start) cover base `start`."""
    start0, end0 = ensembl_to_ucsc(start1, end1)
    return start0, max(end0, start0 + 1)


def intervals_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    """Overlap of two 0-based half-open intervals."""
    return start_a < end_b and start_b < end_a


# =====================================================================
# 2. UCSC Response Parsing
# =====================================================================

def parse_highlight(highlight: str):
    """'hg19.chr22:30771458-30771458#fcfcac' -> ('chr22', 30771458, 30771458) (1-based), or None."""
    match = HIGHLIGHT_RE.match(str(highlight or "").strip())
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def parse_search_loci(payload: dict, rsid: str) -> tuple:
    """Exact rsID loci from a UCSC /search response.

    Returns (dbSNP track used, [(ucsc_chrom, start1, end1), ...]) from the highest dbSNP build
    that contains the rsID, or (None, []) when no dbSNP track matches it.
    """
    target = rsid.strip().lower()
    by_track = {}
    for group in payload.get("positionMatches") or []:
        track = str(group.get("trackName") or group.get("name") or "")
        build = SNP_TRACK_RE.match(track)
        if not build:
            continue
        for match in group.get("matches") or []:
            if str(match.get("posName", "")).strip().lower() != target:
                continue
            locus = parse_highlight(match.get("highlight"))
            if locus:
                by_track.setdefault((int(build.group(1)), track), set()).add(locus)
    if not by_track:
        return None, []
    (_, track), loci = max(by_track.items(), key=lambda item: item[0][0])
    return track, sorted(loci, key=lambda l: (l[0], l[1], l[2]))


def parse_gene_track(payload: dict, track: str) -> list:
    """Transcripts of a UCSC genePred track response (0-based half-open coordinates)."""
    transcripts = []
    for item in payload.get(track) or []:
        starts = [int(x) for x in str(item.get("exonStarts", "")).split(",") if x.strip()]
        ends = [int(x) for x in str(item.get("exonEnds", "")).split(",") if x.strip()]
        transcripts.append({
            "transcript_id": str(item.get("name", "")),
            "gene_symbol": str(item.get("name2", "") or ""),
            "chrom": str(item.get("chrom", "")),
            "strand": str(item.get("strand", "+")),
            "tx_start": int(item["txStart"]),
            "tx_end": int(item["txEnd"]),
            "cds_start": int(item.get("cdsStart", item["txEnd"])),
            "cds_end": int(item.get("cdsEnd", item["txEnd"])),
            "exons": list(zip(starts, ends)),
        })
    return transcripts


def transcript_consequences(tx: dict, start0: int, end0: int) -> list:
    """Positional SO terms of a variant [start0, end0) in a transcript; [] if it is outside."""
    if not intervals_overlap(start0, end0, tx["tx_start"], tx["tx_end"]):
        return []
    in_exon = any(intervals_overlap(start0, end0, s, e) for s, e in tx["exons"])
    coding = tx["cds_start"] < tx["cds_end"]
    terms = set()
    if not coding:
        if in_exon:
            terms.add("non_coding_transcript_exon_variant")
        else:
            terms.update(("intron_variant", "non_coding_transcript_variant"))
    elif not in_exon:
        terms.add("intron_variant")
    else:
        plus = tx["strand"] != "-"
        if intervals_overlap(start0, end0, tx["cds_start"], tx["cds_end"]):
            terms.add("coding_sequence_variant")
        if start0 < tx["cds_start"]:
            terms.add("5_prime_UTR_variant" if plus else "3_prime_UTR_variant")
        if end0 > tx["cds_end"]:
            terms.add("3_prime_UTR_variant" if plus else "5_prime_UTR_variant")
    return sorted(terms)


def strip_version(identifier: str) -> str:
    """'ENSG00000141510.11' -> 'ENSG00000141510' (Ensembl VEP reported unversioned IDs)."""
    return str(identifier).split(".", 1)[0]


def parse_gencode_attrs(lines) -> dict:
    """wgEncodeGencodeAttrsV19 rows -> {transcript ID (versioned and unversioned): attributes}."""
    attrs = {}
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        gene_id, gene_name, gene_type, _, transcript_id, _, transcript_type = fields[:7]
        record = {"gene_id": strip_version(gene_id), "gene_name": gene_name,
                  "gene_type": gene_type, "transcript_type": transcript_type}
        attrs[transcript_id] = record
        attrs.setdefault(strip_version(transcript_id), record)
    return attrs


def build_record(snp_track, loci: list, hits: dict, attrs: dict) -> dict:
    """Annotation record of one rsID (same keys as the former VEP records).

    loci: [(ucsc_chrom, start1, end1)] from dbSNP (None track = rsID not found);
    hits: {locus: [(transcript, [SO terms])]} for the primary loci.
    """
    primary = [l for l in loci if is_primary(l[0])]
    genes, consequences, biotypes = {}, set(), set()
    for locus in primary:
        for tx, terms in hits.get(tuple(locus), []):
            info = attrs.get(tx["transcript_id"]) or attrs.get(strip_version(tx["transcript_id"])) or {}
            symbol = tx["gene_symbol"] or info.get("gene_name") or ""
            gene_key = info.get("gene_id") or symbol or tx["transcript_id"]
            genes[gene_key] = symbol or gene_key
            consequences.update(terms)
            if info.get("transcript_type"):
                biotypes.add(info["transcript_type"])

    if snp_track is None or not loci:
        status = "not_in_dbsnp"
    elif not primary:
        status = "no_primary_locus"
    else:
        status = "genic" if genes else "intergenic"
    return {
        "status": status,
        "genes": sorted(set(genes.values())),
        "gene_ids": sorted(g for g in genes if g.startswith("ENSG")),
        "chr": ", ".join(from_ucsc_chrom(l[0]) for l in primary),
        "pos": ", ".join(str(l[1]) for l in primary),
        "consequences": sorted(consequences),
        "biotypes": sorted(biotypes),
        "release": f"{SOURCE_LABEL} + {snp_track}" if snp_track else SOURCE_LABEL,
    }


def cluster_loci(loci, max_gap: int = CLUSTER_MAX_GAP, max_span: int = CLUSTER_MAX_SPAN) -> list:
    """Group 1-based loci (ucsc_chrom, start1, end1) into region queries.

    Returns [(ucsc_chrom, start0, end0, [loci])] with 0-based half-open regions.
    """
    clusters = []
    for locus in sorted(set(loci), key=lambda l: (l[0], l[1], l[2])):
        start0, end0 = locus_query_interval(locus[1], locus[2])
        if clusters:
            chrom, c_start, c_end, members = clusters[-1]
            if chrom == locus[0] and start0 - c_end <= max_gap and max(end0, c_end) - c_start <= max_span:
                clusters[-1] = (chrom, c_start, max(end0, c_end), members + [locus])
                continue
        clusters.append((locus[0], start0, end0, [locus]))
    return clusters


# =====================================================================
# 3. Network Layer (paced, retried, proxy-independent)
# =====================================================================

def ascii_safe(text: str) -> str:
    """Log-safe text: Windows consoles (cp1252) crash on characters such as U+2019."""
    return str(text).encode("ascii", "replace").decode("ascii")


def describe_response(resp: requests.Response) -> str:
    """'HTTP <status>' plus a one-line message (JSON error, HTML page title, or text)."""
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
    """Session that ignores proxy environment variables and has no hidden urllib3 retries."""
    session = requests.Session()
    session.trust_env = False
    session.proxies = {}
    adapter = HTTPAdapter(pool_connections=pool_size, pool_maxsize=pool_size, max_retries=0)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


class UcscClient:
    """Sequential, paced client for the UCSC Genome Browser REST API."""

    def __init__(self, server: str = DEFAULT_SERVER, session: requests.Session = None,
                 min_interval: float = MIN_INTERVAL):
        self.server = server.rstrip("/")
        self.session = session or create_session()
        self.min_interval = min_interval
        self.requests_sent = 0
        self._last_start = 0.0
        self.last_failure_was_connection = False

    def _pace(self) -> None:
        wait = self._last_start + self.min_interval - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        self._last_start = time.monotonic()

    def request(self, url: str, label: str, params: dict = None, attempts: int = ATTEMPTS):
        """GET with retries. Returns (response or None, error description)."""
        last_error = "no attempt made"
        for attempt in range(attempts):
            self._pace()
            started = time.time()
            resp = None
            self.requests_sent += 1
            try:
                resp = self.session.get(url, params=params, timeout=TIMEOUT)
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

    def get_json(self, path: str, params: dict, label: str, attempts: int = ATTEMPTS):
        """JSON payload of an API endpoint, or (None, error) on failure or an API 'error' field."""
        resp, error = self.request(f"{self.server}{path}", label, params, attempts)
        if resp is None:
            return None, error
        try:
            payload = resp.json()
        except ValueError as exc:
            return None, f"invalid JSON ({exc})"
        if isinstance(payload, dict) and payload.get("error"):
            return None, ascii_safe(f"API error: {payload['error']}")
        return payload, ""

    def check_server(self) -> None:
        """Fail fast when the API (and the GENCODE V19 track) is unreachable; fall back to the
        environment proxy only when the direct connection fails and one is configured."""
        params = {"genome": GENOME, "track": GENE_TRACKS[0]}
        payload, error = self.get_json("/list/schema", params, "server check", attempts=3)
        if payload is None and self.last_failure_was_connection:
            env_proxies = requests.utils.get_environ_proxies(self.server)
            if env_proxies:
                print(f"{LOG} Direct connection to {self.server} failed ({error}); "
                      f"retrying through the proxy set in the environment.", flush=True)
                self.session.proxies = dict(env_proxies)
                payload, error = self.get_json("/list/schema", params, "server check (proxy)", attempts=3)
                if payload is None:
                    self.session.proxies = {}
        if payload is None:
            raise RuntimeError(
                f"{LOG} UCSC API {self.server} is unavailable ({error}). "
                "No SNPs were queried; re-run when the server is reachable."
            )

    def search_rsid(self, rsid: str):
        """(dbSNP track, loci) of an rsID, or None when the request failed."""
        payload, _ = self.get_json("/search", {"search": rsid, "genome": GENOME}, f"search {rsid}")
        return None if payload is None else parse_search_loci(payload, rsid)

    def transcripts(self, track: str, chrom: str, start0: int = None, end0: int = None):
        """Transcripts of `track` overlapping [start0, end0) on `chrom` (whole chromosome when no
        region is given), or None on failure."""
        params = {"genome": GENOME, "track": track, "chrom": chrom}
        label = f"{track} {chrom}"
        if start0 is not None:
            params.update(start=start0, end=end0)
            label += f":{start0 + 1}-{end0}"
        payload, _ = self.get_json("/getData/track", params, label)
        return None if payload is None else parse_gene_track(payload, track)


# =====================================================================
# 4. Pipeline Steps
# =====================================================================

def significant_snps(input_paths: list) -> list:
    """Unique SNPs BH-significant in at least one Step B1 overlap table."""
    snps = set()
    for path in input_paths:
        df = pd.read_csv(path, sep="\t", usecols=["SNP", "BH_SIGNIFICANT"], dtype=str)
        snps.update(df.loc[is_true(df["BH_SIGNIFICANT"]), "SNP"].str.strip())
    return sorted(snps)


def resolve_server(server: str) -> str:
    """Use the UCSC API; an Ensembl URL (legacy `rna.ensembl_server` value) maps to it."""
    server = (server or DEFAULT_SERVER).rstrip("/")
    if "ensembl.org" in (urlparse(server).hostname or ""):
        print(f"{LOG} --server {server} is the legacy Ensembl setting; using {DEFAULT_SERVER}.", flush=True)
        return DEFAULT_SERVER
    return server


def resolve_positions(client: UcscClient, rsids: list, positions: dict, positions_path: str,
                      batch_size: int) -> list:
    """Look up hg19 loci of `rsids` (stored in `positions`, flushed every batch). Returns failures."""
    failed, pending = [], {}
    started = time.time()
    try:
        for done, rsid in enumerate(rsids, start=1):
            result = client.search_rsid(rsid)
            if result is None:
                failed.append(rsid)
            else:
                track, loci = result
                pending[rsid] = {"track": track, "loci": [list(l) for l in loci]}
            if len(pending) >= batch_size:
                positions.update(pending)
                save_json_cache(positions_path, pending)
                pending = {}
            if done % batch_size == 0 or done == len(rsids):
                rate = (time.time() - started) / done
                print(f"{LOG} rsID positions: {done:,}/{len(rsids):,} looked up | failed {len(failed):,} | "
                      f"~{rate * (len(rsids) - done) / 60:.0f} min left", flush=True)
    finally:
        if pending:
            positions.update(pending)
            save_json_cache(positions_path, pending)
    return failed


def load_gencode_attrs(cache_dir: str, client: UcscClient, attrs_url: str = ATTRS_URL) -> dict:
    """GENCODE V19 transcript attributes, downloaded once into `cache_dir`."""
    path = os.path.join(cache_dir, ATTRS_CACHE_NAME)
    if not os.path.exists(path):
        resp, error = client.request(attrs_url, "GENCODE V19 attributes download")
        if resp is None:
            raise RuntimeError(f"{LOG} Could not download {attrs_url} ({error}); re-run later.")
        os.makedirs(cache_dir, exist_ok=True)
        tmp = f"{path}.{os.getpid()}.tmp"
        with open(tmp, "wb") as f:
            f.write(resp.content)
        os.replace(tmp, path)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        attrs = parse_gencode_attrs(f)
    if not attrs:
        raise RuntimeError(f"{LOG} {path} holds no GENCODE attributes; delete it and re-run.")
    return attrs


def annotate_genes(client: UcscClient, rsids: list, positions: dict, attrs: dict, cache: dict,
                   cache_path: str, batch_size: int) -> list:
    """Gene records of resolved `rsids` into `cache` (flushed every batch). Returns failures."""
    pending, failed = {}, set()

    def flush():
        if pending:
            cache.update(pending)
            save_json_cache(cache_path, dict(pending))
            pending.clear()

    loci_of = {r: [tuple(l) for l in positions[r]["loci"]] for r in rsids}
    remaining = {r: {l for l in loci if is_primary(l[0])} for r, loci in loci_of.items()}
    hits = {r: {} for r in rsids}
    whole_chromosome = {}                      # (track, chrom) -> transcripts, fetched once
    owners = {}
    for rsid, loci in remaining.items():
        for locus in loci:
            owners.setdefault(locus, []).append(rsid)

    try:
        # SNPs without a primary locus need no gene query
        for rsid in rsids:
            if not remaining[rsid]:
                pending[rsid] = build_record(positions[rsid]["track"], loci_of[rsid], {}, attrs)
        clusters = cluster_loci(owners)
        print(f"{LOG} Gene lookup: {len(owners):,} loci in {len(clusters):,} regions on "
              f"{len({c[0] for c in clusters})} chromosomes (GENCODE V19: {', '.join(GENE_TRACKS)})", flush=True)
        for n, (chrom, start0, end0, members) in enumerate(clusters, start=1):
            transcripts = []
            for track in GENE_TRACKS:
                if track in WHOLE_CHROMOSOME_TRACKS:
                    if (track, chrom) not in whole_chromosome:
                        whole_chromosome[(track, chrom)] = client.transcripts(track, chrom)
                    found = whole_chromosome[(track, chrom)]
                    if found is not None:
                        found = [t for t in found if intervals_overlap(start0, end0, t["tx_start"], t["tx_end"])]
                    else:
                        del whole_chromosome[(track, chrom)]   # retry on this chromosome's next region
                else:
                    found = client.transcripts(track, chrom, start0, end0)
                if found is None:
                    transcripts = None
                    break
                transcripts.extend(found)
            for locus in members:
                for rsid in owners[locus]:
                    if transcripts is None:
                        failed.add(rsid)
                        continue
                    q_start, q_end = locus_query_interval(locus[1], locus[2])
                    hits[rsid][locus] = [
                        (tx, terms) for tx in transcripts
                        if tx["chrom"] == chrom and (terms := transcript_consequences(tx, q_start, q_end))
                    ]
                    remaining[rsid].discard(locus)
                    if not remaining[rsid] and rsid not in failed:
                        pending[rsid] = build_record(positions[rsid]["track"], loci_of[rsid], hits[rsid], attrs)
            if len(pending) >= batch_size:
                flush()
            if n % 50 == 0 or n == len(clusters):
                print(f"{LOG} Gene lookup: {n:,}/{len(clusters):,} regions done", flush=True)
    finally:
        flush()
    return sorted(failed)


def annotation_row(snp: str, record: dict) -> dict:
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


def run_annotate_rna_genes(
    input_paths: list,
    output_path: str,
    server: str = DEFAULT_SERVER,
    cache_path: str = DEFAULT_CACHE_FILE,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_workers: int = 1,
    attrs_url: str = ATTRS_URL,
    client: UcscClient = None,
) -> pd.DataFrame:
    server = resolve_server(server)
    snps = significant_snps(input_paths)
    print(f"{LOG} {len(input_paths)} traits: {len(snps):,} unique BH-significant SNPs to annotate.")

    cache = load_json_cache(cache_path)
    uncached = [s for s in snps if s not in cache]
    print(f"{LOG} Cache hits: {len(snps) - len(uncached):,} | Uncached SNPs to query: {len(uncached):,}")

    if uncached:
        if max_workers != 1:
            print(f"{LOG} UCSC asks for sequential requests; --workers {max_workers} is ignored.")
        client = client or UcscClient(server)
        client.check_server()
        cache_dir = os.path.dirname(os.path.abspath(cache_path))
        positions_path = os.path.join(cache_dir, POSITIONS_CACHE_NAME)
        positions = load_json_cache(positions_path)

        to_search = [s for s in uncached if s not in positions]
        print(f"{LOG} {SOURCE_LABEL} via {server}: {len(uncached) - len(to_search):,} rsID positions "
              f"cached, {len(to_search):,} to look up (~1 request/s)", flush=True)
        failed = resolve_positions(client, to_search, positions, positions_path, batch_size)

        attrs = load_gencode_attrs(cache_dir, client, attrs_url)
        resolved = [s for s in uncached if s in positions]
        failed += annotate_genes(client, resolved, positions, attrs, cache, cache_path, batch_size)
        n_new = sum(1 for s in uncached if s in cache)
        print(f"{LOG} Cached {n_new:,} newly annotated SNPs at {cache_path} ({client.requests_sent:,} requests)")

        if failed:
            raise RuntimeError(
                f"{LOG} UCSC lookups failed for {len(failed):,} SNPs (e.g. {sorted(failed)[:5]}). "
                f"Successful results were cached; re-run the pipeline to retry the remaining SNPs."
            )

    out = pd.DataFrame([annotation_row(s, cache[s]) for s in snps], columns=OUTPUT_COLUMNS)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    counts = out["STATUS"].value_counts().to_dict()
    print(f"{LOG} {len(out):,} SNPs {counts} -> {output_path}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate BH-significant RNA m6A SNPs with the GENCODE V19 genes they fall inside (UCSC hg19)."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Step B1 results/rna/{trait}_overlap_bh.tsv files")
    parser.add_argument("--output", "-o", required=True, help="Output results/rna/rna_sig_genes.tsv")
    parser.add_argument(
        "--server", default=DEFAULT_SERVER,
        help=f"UCSC Genome Browser REST API (default: {DEFAULT_SERVER}; an Ensembl URL maps to it)",
    )
    parser.add_argument("--cache-file", default=DEFAULT_CACHE_FILE, help="Persistent JSON annotation cache path")
    parser.add_argument(
        "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help=f"SNPs per cache checkpoint (default: {DEFAULT_BATCH_SIZE})",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Kept for compatibility; UCSC requests are sent sequentially (~1 per second)",
    )
    args = parser.parse_args()

    run_annotate_rna_genes(
        args.inputs,
        args.output,
        server=args.server,
        cache_path=args.cache_file,
        batch_size=max(1, args.batch_size),
        max_workers=args.workers,
    )


if __name__ == "__main__":
    main()
