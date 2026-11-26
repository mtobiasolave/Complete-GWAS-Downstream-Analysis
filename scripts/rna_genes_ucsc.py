#!/usr/bin/env python3
"""
Step B2 backend "ucsc": genes of BH-significant RNA m6A SNPs from the UCSC Genome Browser REST API
(hg19 = GRCh37, GENCODE V19). Called by annotate_rna_genes.py (--source ucsc).

1. Resolves each uncached rsID to its hg19 position(s) with the search endpoint
   (`/search?search=<rsID>&genome=hg19`), using the newest dbSNP track that contains the rsID
   (dbSnp155, then older builds). Only loci on primary chromosomes (1-22, X, Y, MT) are used.
2. Fetches the GENCODE V19 transcripts around those SNPs (`/getData/track`, tracks
   wgEncodeGencodeCompV19 + wgEncodeGencodePseudoGeneV19; nearby SNPs share one region query)
   and keeps the genes whose transcripts the variant physically falls inside (exon or intron;
   no upstream/downstream window), all biotypes: the same rule as the Ensembl backend.
   Gene IDs and transcript biotypes come from UCSC's GENCODE V19 attributes table
   (wgEncodeGencodeAttrsV19, downloaded once from hgdownload.soe.ucsc.edu).
3. rsID positions are cached in ucsc_hg19_rsid_positions.json next to the annotation cache;
   both are flushed every batch and on exit.

Coordinates: records use Ensembl-style chromosome names (1..22, X, Y, MT) and 1-based POS, like
the Ensembl backend. UCSC uses chr-prefixed names and 0-based half-open intervals; conversions
are explicit below.

Network: UCSC asks for sequential requests and about one request per second, so requests are
paced (rna_genes_common.ApiClient provides retries, backoff and proxy handling).

STATUS: genic | intergenic | not_in_dbsnp (rsID in no UCSC hg19 dbSNP track) |
        no_primary_locus (maps only to alternative haplotypes).
CONSEQUENCES are positional Sequence Ontology terms (coding_sequence_variant,
5_prime_UTR_variant, 3_prime_UTR_variant, intron_variant, non_coding_transcript_exon_variant,
non_coding_transcript_variant); allele-dependent terms (e.g. missense_variant) need the
alleles and are not derived.
"""

import gzip
import os
import re
import time

from common_utils import load_json_cache, save_json_cache
from rna_genes_common import LOG, PRIMARY_CHROMS, ApiClient

SOURCE_NAME = "UCSC"
DEFAULT_SERVER = "https://api.genome.ucsc.edu"
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

CLUSTER_MAX_GAP = 500_000        # SNPs closer than this share one region query ...
CLUSTER_MAX_SPAN = 5_000_000     # ... as long as the region stays below this size
MIN_INTERVAL = 1.0               # seconds between request starts (UCSC: ~1 request per second)
TIMEOUT = (15, 180)              # (connect, read) seconds
ATTEMPTS = 6


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
    """'ENSG00000141510.11' -> 'ENSG00000141510' (the Ensembl backend reports unversioned IDs)."""
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
    """Annotation record of one rsID (same keys as the Ensembl backend's records).

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
# 3. UCSC Client (paced, retried, proxy-independent)
# =====================================================================

class UcscClient(ApiClient):
    """Sequential, paced client for the UCSC Genome Browser REST API."""

    name = "UCSC API"

    def __init__(self, server: str = DEFAULT_SERVER, session=None, min_interval: float = MIN_INTERVAL,
                 timeout: tuple = TIMEOUT):
        super().__init__(server, timeout, ATTEMPTS, min_interval=min_interval, session=session)

    def check_server(self, allow_proxy_fallback: bool = True) -> None:
        """Fail fast when the API (and the GENCODE V19 track) is unreachable."""
        self.check("/list/schema", {"genome": GENOME, "track": GENE_TRACKS[0]},
                   lambda payload: isinstance(payload, dict), allow_proxy_fallback=allow_proxy_fallback)

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
        resp, error = client.request("GET", attrs_url, "GENCODE V19 attributes download")
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


def annotate(snps: list, cache: dict, cache_path: str, server: str = DEFAULT_SERVER, batch_size: int = 25,
             max_workers: int = None, attrs_url: str = None, client: UcscClient = None) -> list:
    """Annotate uncached `snps` into `cache` (+ its file, flushed every batch). Returns failed rsIDs."""
    if max_workers not in (None, 1):
        print(f"{LOG} UCSC asks for sequential requests; --workers {max_workers} is ignored.")
    client = client or UcscClient(server)
    client.check_server()
    cache_dir = os.path.dirname(os.path.abspath(cache_path))
    positions_path = os.path.join(cache_dir, POSITIONS_CACHE_NAME)
    positions = load_json_cache(positions_path)

    to_search = [s for s in snps if s not in positions]
    print(f"{LOG} {SOURCE_LABEL} via {client.server}: {len(snps) - len(to_search):,} rsID positions "
          f"cached, {len(to_search):,} to look up (~1 request/s)", flush=True)
    failed = resolve_positions(client, to_search, positions, positions_path, batch_size)

    attrs = load_gencode_attrs(cache_dir, client, attrs_url or ATTRS_URL)
    resolved = [s for s in snps if s in positions]
    failed += annotate_genes(client, resolved, positions, attrs, cache, cache_path, batch_size)
    n_new = sum(1 for s in snps if s in cache)
    print(f"{LOG} Cached {n_new:,} newly annotated SNPs at {cache_path} ({client.requests_sent:,} requests)")
    return sorted(set(failed))
