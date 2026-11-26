"""Tests for Step B2 (scripts/annotate_rna_genes.py and its backends rna_genes_ucsc.py,
rna_genes_ensembl.py, rna_genes_common.py): RNA m6A SNP gene annotation from UCSC hg19 /
GENCODE V19 or Ensembl GRCh37 VEP.

Run from the repository root with snakemake9_env activated:
    python -m unittest discover -s tests -v
Live smoke tests are skipped unless RUN_LIVE_UCSC=1 (api.genome.ucsc.edu) and/or
RUN_LIVE_ENSEMBL=1 (grch37.rest.ensembl.org) are set:
    RUN_LIVE_UCSC=1 RUN_LIVE_ENSEMBL=1 python -m unittest discover -s tests -k Live -v
"""

import gzip
import json
import os
import socket
import sys
import tempfile
import threading
import time
import unittest
from collections import deque
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlparse

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import annotate_rna_genes as A  # noqa: E402
import rna_genes_common as C  # noqa: E402
import rna_genes_ensembl as E  # noqa: E402
import rna_genes_ucsc as U  # noqa: E402

FAST = {"BACKOFF_BASE": 0.01, "BACKOFF_CAP": 0.05}


def dump_json(path, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# =====================================================================
# Fixtures: synthetic GENCODE V19 tracks, dbSNP searches, attributes, VEP answers
# =====================================================================

def gene_item(name, symbol, chrom, strand, tx, cds, exons):
    return {
        "bin": 585, "name": name, "chrom": chrom, "strand": strand,
        "txStart": tx[0], "txEnd": tx[1], "cdsStart": cds[0], "cdsEnd": cds[1],
        "exonCount": len(exons),
        "exonStarts": "".join(f"{s}," for s, _ in exons),
        "exonEnds": "".join(f"{e}," for _, e in exons),
        "score": 0, "name2": symbol, "cdsStartStat": "cmpl", "cdsEndStat": "cmpl", "exonFrames": "",
    }


TRACK_ITEMS = {
    "wgEncodeGencodeCompV19": [
        # protein coding, + strand: exons [1000,1600) and [3000,5000), CDS [1500,4000)
        gene_item("ENST00000000001.1", "GENEA", "chr1", "+", (1000, 5000), (1500, 4000), [(1000, 1600), (3000, 5000)]),
        # non-coding (cdsStart == cdsEnd), exons [20000,20100) and [20900,21000)
        gene_item("ENST00000000002.1", "GENEB", "chr1", "-", (20000, 21000), (21000, 21000), [(20000, 20100), (20900, 21000)]),
    ],
    "wgEncodeGencodePseudoGeneV19": [
        gene_item("ENST00000000003.1", "PSEUDO1", "chr2", "+", (500, 800), (800, 800), [(500, 800)]),
    ],
}
ATTRS_LINES = [
    "ENSG00000000001.3\tGENEA\tprotein_coding\tKNOWN\tENST00000000001.1\tGENEA-001\tprotein_coding\tKNOWN\tH1\tH1\t\t2\tcoding\tENSP1",
    "ENSG00000000002.1\tGENEB\tlincRNA\tKNOWN\tENST00000000002.1\tGENEB-001\tlincRNA\tKNOWN\tH2\tH2\t\t2\tnonCoding\t",
    "ENSG00000000003.2\tPSEUDO1\tpseudogene\tKNOWN\tENST00000000003.1\tPSEUDO1-001\tprocessed_pseudogene\tKNOWN\tH3\tH3\t\t1\tpseudo\t",
]


def snp_match(track, rsid, chrom, pos):
    return {"name": track, "trackName": track, "matches": [{
        "position": f"{chrom}:{pos - 100}-{pos + 100}", "hgFindMatches": rsid, "posName": rsid,
        "highlight": f"hg19.{chrom}:{pos}-{pos}#fcfcac", "canonical": False}]}


SEARCH = {
    "rs1": [snp_match("dbSnp155", "rs1", "chr1", 1001)],     # 1st exon base, before CDS on + -> 5' UTR
    "rs2": [snp_match("dbSnp155", "rs2", "chr1", 2000)],     # intron of GENEA
    "rs3": [snp_match("dbSnp155", "rs3", "chr1", 20950)],    # non-coding exon of GENEB
    "rs4": [snp_match("snp151", "rs4", "chr1", 10000)],      # intergenic; only in an older dbSNP build
    "rs5": [],                                               # unknown rsID
    "rs6": [snp_match("dbSnp155", "rs6", "chr6_cox_hap2", 5000)],  # alternative haplotype only
    "rs7": [snp_match("dbSnp155", "rs7", "chr2", 800)],      # last base of PSEUDO1 (0-based 799)
    "rs8": [snp_match("dbSnp155", "rs8", "chr2", 801)],      # first base after PSEUDO1 -> intergenic
    "rs9": [snp_match("dbSnp155", "rs9", "chr1", 1001), snp_match("dbSnp155", "rs9", "chr2", 800)],
    "rs10": [snp_match("dbSnp155", "rs10", "chr1", 900000)],  # > 500 kb away: 2nd chr1 region
}


def vep_entry(rsid, chrom, pos, consequences):
    return {"input": rsid, "id": rsid, "seq_region_name": chrom, "start": pos, "end": pos,
            "transcript_consequences": [{"gene_id": g, "gene_symbol": s, "consequence_terms": t, "biotype": b}
                                        for g, s, t, b in consequences]}


GENEA_UTR = ("ENSG00000000001", "GENEA", ["5_prime_UTR_variant"], "protein_coding")
PSEUDO1_EXON = ("ENSG00000000003", "PSEUDO1", ["non_coding_transcript_exon_variant"], "processed_pseudogene")
# The same SNPs as SEARCH, as Ensembl VEP reports them (rs5 is unknown: VEP omits it)
VEP = {
    "rs1": [vep_entry("rs1", "1", 1001, [GENEA_UTR, ("ENSG00000000008", "NEARBY", ["upstream_gene_variant"], "protein_coding")])],
    "rs2": [vep_entry("rs2", "1", 2000, [("ENSG00000000001", "GENEA", ["intron_variant"], "protein_coding")])],
    "rs3": [vep_entry("rs3", "1", 20950, [("ENSG00000000002", "GENEB", ["non_coding_transcript_exon_variant"], "lincRNA")])],
    "rs4": [vep_entry("rs4", "1", 10000, [])],
    "rs6": [vep_entry("rs6", "HSCHR6_MHC_COX", 5000, [("ENSG00000000007", "HAPGENE", ["intron_variant"], "protein_coding")])],
    "rs7": [vep_entry("rs7", "2", 800, [PSEUDO1_EXON])],
    "rs8": [vep_entry("rs8", "2", 801, [("ENSG00000000003", "PSEUDO1", ["downstream_gene_variant"], "processed_pseudogene")])],
    "rs9": [vep_entry("rs9", "1", 1001, [GENEA_UTR]), vep_entry("rs9", "2", 800, [PSEUDO1_EXON])],
    "rs10": [vep_entry("rs10", "1", 900000, [])],
}
VEP_RELEASE = 75
LEGACY_VEP_RECORD = {
    "status": "genic", "genes": ["OLDGENE"], "gene_ids": ["ENSG00000999999"], "chr": "3", "pos": "12345",
    "consequences": ["intron_variant"], "biotypes": ["protein_coding"], "release": "116",
}
UCSC_RECORD = {
    "status": "genic", "genes": ["UGENE"], "gene_ids": ["ENSG00000888888"], "chr": "5", "pos": "777",
    "consequences": ["intron_variant"], "biotypes": ["lincRNA"], "release": "UCSC hg19 GENCODE V19 + dbSnp155",
}
# Expected (STATUS, GENES, ENSEMBL_GENE_IDS, CHR, POS, CONSEQUENCES, TRANSCRIPT_BIOTYPES), both backends
EXPECTED_SHARED = {
    "rs1": ("genic", "GENEA", "ENSG00000000001", "1", "1001", "5_prime_UTR_variant", "protein_coding"),
    "rs2": ("genic", "GENEA", "ENSG00000000001", "1", "2000", "intron_variant", "protein_coding"),
    "rs3": ("genic", "GENEB", "ENSG00000000002", "1", "20950", "non_coding_transcript_exon_variant", "lincRNA"),
    "rs4": ("intergenic", "intergenic", "NA", "1", "10000", "NA", "NA"),
    "rs6": ("no_primary_locus", "no primary-assembly locus", "NA", "NA", "NA", "NA", "NA"),
    "rs7": ("genic", "PSEUDO1", "ENSG00000000003", "2", "800", "non_coding_transcript_exon_variant", "processed_pseudogene"),
    "rs8": ("intergenic", "intergenic", "NA", "2", "801", "NA", "NA"),
    "rs9": ("genic", "GENEA, PSEUDO1", "ENSG00000000001, ENSG00000000003", "1, 2", "1001, 800",
            "5_prime_UTR_variant, non_coding_transcript_exon_variant", "processed_pseudogene, protein_coding"),
    "rs10": ("intergenic", "intergenic", "NA", "1", "900000", "NA", "NA"),
}
VALUE_COLUMNS = ["STATUS", "GENES", "ENSEMBL_GENE_IDS", "CHR", "POS", "CONSEQUENCES", "TRANSCRIPT_BIOTYPES"]


class MockApi(BaseHTTPRequestHandler):
    """Minimal UCSC API (/list/schema, /search, /getData/track, /attrs.txt.gz) and Ensembl REST
    (/info/ping, /info/data, POST /vep/homo_sapiens/id), plus a scripted /flaky endpoint."""
    hits = []
    fail_chroms = set()
    script = deque()
    vep_script = deque()
    vep_fail_ids = set()
    vep_batches = []
    ping_down = False

    @classmethod
    def reset(cls):
        cls.hits.clear()
        cls.fail_chroms.clear()
        cls.script.clear()
        cls.vep_script.clear()
        cls.vep_fail_ids.clear()
        cls.vep_batches.clear()
        cls.ping_down = False

    def log_message(self, *args):
        pass

    def reply(self, code, body, ctype="application/json", headers=None):
        data = body if isinstance(body, bytes) else body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        for key, value in (headers or {}).items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        url = urlparse(self.path)            # absolute URI when the request came through a proxy
        q = {k: v[0] for k, v in parse_qs(url.query).items()}
        MockApi.hits.append(url.path)
        if url.path == "/flaky":
            code, headers, body, delay = MockApi.script.popleft()
            time.sleep(delay)
            return self.reply(code, body, headers=headers)
        if url.path == "/list/schema":
            return self.reply(200, json.dumps({"track": q.get("track"), "columnTypes": []}))
        if url.path == "/search":
            return self.reply(200, json.dumps({"genome": "hg19", "positionMatches": SEARCH.get(q["search"], [])}))
        if url.path == "/getData/track":
            if q["chrom"] in MockApi.fail_chroms:
                return self.reply(503, json.dumps({"error": "busy"}), headers={"Retry-After": "0"})
            MockApi.hits.append(f"getData:{q['track']}:{q['chrom']}:{'region' if 'start' in q else 'chromosome'}")
            start, end = int(q.get("start", 0)), int(q.get("end", 10**12))
            items = [i for i in TRACK_ITEMS.get(q["track"], [])
                     if i["chrom"] == q["chrom"] and i["txStart"] < end and i["txEnd"] > start]
            return self.reply(200, json.dumps({"track": q["track"], q["track"]: items}))
        if url.path == "/attrs.txt.gz":
            return self.reply(200, gzip.compress(("\n".join(ATTRS_LINES) + "\n").encode()), "application/gzip")
        if url.path == "/info/ping":
            if MockApi.ping_down:
                return self.reply(500, "<html><title>Error: 500 – service unavailable</title></html>", "text/html")
            return self.reply(200, json.dumps({"ping": 1}))
        if url.path == "/info/data":
            return self.reply(200, json.dumps({"releases": [VEP_RELEASE]}))
        return self.reply(400, json.dumps({"error": f"unknown endpoint {url.path}"}))

    def do_POST(self):
        url = urlparse(self.path)
        body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        MockApi.hits.append(url.path)
        if url.path != "/vep/homo_sapiens/id":
            return self.reply(400, json.dumps({"error": f"unknown endpoint {url.path}"}))
        ids = list(body.get("ids", []))
        MockApi.vep_batches.append(ids)
        if MockApi.vep_script:
            code, headers = MockApi.vep_script.popleft()
            return self.reply(code, json.dumps({"error": "busy"}), headers=headers)
        if set(ids) & MockApi.vep_fail_ids:
            return self.reply(500, "<html><title>500 Internal Server Error</title></html>", "text/html")
        return self.reply(200, json.dumps([e for i in ids for e in VEP.get(i, [])]))


class QuietServer(ThreadingHTTPServer):
    def handle_error(self, request, client_address):
        # A client that timed out on purpose closes its socket before the late reply
        if not isinstance(sys.exc_info()[1], (ConnectionAbortedError, ConnectionResetError, BrokenPipeError)):
            super().handle_error(request, client_address)


def free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class MockServerCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = QuietServer(("127.0.0.1", 0), MockApi)
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()
        cls.url = f"http://127.0.0.1:{cls.httpd.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def setUp(self):
        MockApi.reset()
        patcher = mock.patch.multiple(C, **FAST)
        patcher.start()
        self.addCleanup(patcher.stop)

    def ucsc_client(self, url=None):
        return U.UcscClient(url or self.url, min_interval=0)

    def vep_client(self, url=None):
        return E.VepClient(url or self.url)

    def clients(self, url=None):
        return {"ucsc": self.ucsc_client(url), "ensembl": self.vep_client(url)}


# =====================================================================
# 1. Chromosome names
# =====================================================================

class ChromosomeTests(unittest.TestCase):
    def test_to_ucsc(self):
        for bare, ucsc in [("1", "chr1"), ("22", "chr22"), ("X", "chrX"), ("x", "chrX"), ("Y", "chrY"),
                           ("MT", "chrM"), ("M", "chrM"), ("chr7", "chr7"), ("chrM", "chrM"), (17, "chr17")]:
            self.assertEqual(U.to_ucsc_chrom(bare), ucsc, bare)

    def test_from_ucsc(self):
        for ucsc, bare in [("chr1", "1"), ("chr22", "22"), ("chrX", "X"), ("chrM", "MT"), ("chrMT", "MT"), ("5", "5")]:
            self.assertEqual(U.from_ucsc_chrom(ucsc), bare, ucsc)

    def test_round_trip(self):
        for bare in sorted(C.PRIMARY_CHROMS):
            self.assertEqual(U.from_ucsc_chrom(U.to_ucsc_chrom(bare)), bare)

    def test_primary(self):
        self.assertTrue(U.is_primary("chr1") and U.is_primary("chrM") and U.is_primary("MT"))
        for name in ("chr6_cox_hap2", "chr17_gl000205_random", "chrUn_gl000220", "GL000192.1"):
            self.assertFalse(U.is_primary(name), name)


# =====================================================================
# 2. Coordinates
# =====================================================================

class CoordinateTests(unittest.TestCase):
    def test_conversion(self):
        self.assertEqual(U.ensembl_to_ucsc(7565096, 7590856), (7565095, 7590856))
        self.assertEqual(U.ucsc_to_ensembl(7565095, 7590856), (7565096, 7590856))
        self.assertEqual(U.ucsc_to_ensembl(*U.ensembl_to_ucsc(100, 100)), (100, 100))
        self.assertEqual(U.ensembl_to_ucsc(30771458, 30771458), (30771457, 30771458))  # rs1002189

    def test_single_base_interval(self):
        self.assertEqual(U.locus_query_interval(1001, 1001), (1000, 1001))
        self.assertEqual(U.locus_query_interval(1001, 1000), (1000, 1001))   # insertion: end < start

    def test_transcript_boundaries(self):
        # UCSC transcript [100, 200) covers 1-based bases 101..200
        tx = U.parse_gene_track({"t": [gene_item("ENST1.1", "G", "chr1", "+", (100, 200), (200, 200), [(100, 200)])]}, "t")[0]
        inside = {pos: bool(U.transcript_consequences(tx, *U.locus_query_interval(pos, pos)))
                  for pos in (100, 101, 200, 201)}
        self.assertEqual(inside, {100: False, 101: True, 200: True, 201: False})

    def test_clusters(self):
        loci = [("chr1", 1001, 1001), ("chr1", 50000, 50000), ("chr1", 400000, 400000), ("chr2", 10, 10)]
        clusters = U.cluster_loci(loci, max_gap=100_000, max_span=2_000_000)
        self.assertEqual([(c[0], c[1], c[2], len(c[3])) for c in clusters],
                         [("chr1", 1000, 50000, 2), ("chr1", 399999, 400000, 1), ("chr2", 9, 10, 1)])
        self.assertEqual(len(U.cluster_loci(loci[:2], max_gap=100_000, max_span=10_000)), 2)


# =====================================================================
# 3. JSON parsing & record mapping (UCSC and Ensembl)
# =====================================================================

class ParsingTests(unittest.TestCase):
    def test_highlight(self):
        self.assertEqual(U.parse_highlight("hg19.chr22:30771458-30771458#fcfcac"), ("chr22", 30771458, 30771458))
        self.assertEqual(U.parse_highlight("chr6_cox_hap2:2377304-2377305"), ("chr6_cox_hap2", 2377304, 2377305))
        self.assertIsNone(U.parse_highlight(""))
        self.assertIsNone(U.parse_highlight(None))

    def test_search_prefers_newest_dbsnp_and_exact_name(self):
        payload = {"positionMatches": [
            {"trackName": "gtexEqtlCluster", "matches": [{"posName": "rs1002189", "highlight": "hg19.chr22:1-1"}]},
            snp_match("snp151", "rs1002189", "chr22", 30771458),
            snp_match("dbSnp155Common", "rs1002189", "chr22", 999),        # subset track: ignored
            snp_match("dbSnp155", "rs1002189", "chr22", 30771458),
            snp_match("dbSnp155", "rs10021890", "chr4", 5),               # other rsID: ignored
        ]}
        self.assertEqual(U.parse_search_loci(payload, "rs1002189"), ("dbSnp155", [("chr22", 30771458, 30771458)]))
        self.assertEqual(U.parse_search_loci({"positionMatches": [snp_match("snp151", "rsX", "chr1", 5)]}, "RSX"),
                         ("snp151", [("chr1", 5, 5)]))
        self.assertEqual(U.parse_search_loci({"genome": "hg19"}, "rs1"), (None, []))

    def test_search_keeps_all_loci(self):
        payload = {"positionMatches": [snp_match("dbSnp155", "rs9", "chr1", 7), snp_match("dbSnp155", "rs9", "chr6_cox_hap2", 3)]}
        self.assertEqual(U.parse_search_loci(payload, "rs9")[1], [("chr1", 7, 7), ("chr6_cox_hap2", 3, 3)])

    def test_gene_track_fields(self):
        txs = U.parse_gene_track({"wgEncodeGencodeCompV19": TRACK_ITEMS["wgEncodeGencodeCompV19"]}, "wgEncodeGencodeCompV19")
        self.assertEqual(txs[0]["transcript_id"], "ENST00000000001.1")
        self.assertEqual(txs[0]["gene_symbol"], "GENEA")
        self.assertEqual((txs[0]["tx_start"], txs[0]["tx_end"], txs[0]["cds_start"], txs[0]["cds_end"]), (1000, 5000, 1500, 4000))
        self.assertEqual(txs[0]["exons"], [(1000, 1600), (3000, 5000)])
        self.assertEqual(U.parse_gene_track({"other": []}, "wgEncodeGencodeCompV19"), [])

    def test_consequences(self):
        coding, noncoding = U.parse_gene_track({"t": TRACK_ITEMS["wgEncodeGencodeCompV19"]}, "t")
        cons = lambda tx, pos: U.transcript_consequences(tx, *U.locus_query_interval(pos, pos))
        self.assertEqual(cons(coding, 1001), ["5_prime_UTR_variant"])
        self.assertEqual(cons(coding, 1501), ["coding_sequence_variant"])
        self.assertEqual(cons(coding, 2000), ["intron_variant"])
        self.assertEqual(cons(coding, 4500), ["3_prime_UTR_variant"])
        self.assertEqual(cons(coding, 6000), [])
        self.assertEqual(cons(noncoding, 20950), ["non_coding_transcript_exon_variant"])
        self.assertEqual(cons(noncoding, 20500), ["intron_variant", "non_coding_transcript_variant"])
        minus = dict(coding, strand="-")                        # UTR sides swap on the - strand
        self.assertEqual(cons(minus, 1001), ["3_prime_UTR_variant"])
        self.assertEqual(cons(minus, 4500), ["5_prime_UTR_variant"])

    def test_attrs_and_record(self):
        attrs = U.parse_gencode_attrs(line + "\n" for line in ATTRS_LINES)
        self.assertEqual(attrs["ENST00000000001.1"]["gene_id"], "ENSG00000000001")
        self.assertIs(attrs["ENST00000000001"], attrs["ENST00000000001.1"])
        tx = U.parse_gene_track({"t": TRACK_ITEMS["wgEncodeGencodeCompV19"]}, "t")[0]
        locus = ("chr1", 1001, 1001)
        record = U.build_record("dbSnp155", [locus, ("chr6_cox_hap2", 1, 1)], {locus: [(tx, ["5_prime_UTR_variant"])]}, attrs)
        self.assertEqual(record, {
            "status": "genic", "genes": ["GENEA"], "gene_ids": ["ENSG00000000001"], "chr": "1", "pos": "1001",
            "consequences": ["5_prime_UTR_variant"], "biotypes": ["protein_coding"],
            "release": "UCSC hg19 GENCODE V19 + dbSnp155",
        })
        self.assertEqual(U.build_record(None, [], {}, attrs)["status"], "not_in_dbsnp")
        self.assertEqual(U.build_record("dbSnp155", [("chr6_cox_hap2", 1, 1)], {}, attrs)["status"], "no_primary_locus")
        self.assertEqual(U.build_record("dbSnp155", [locus], {locus: []}, attrs)["status"], "intergenic")

    def test_vep_entries(self):
        record = E.summarize_vep_entries(VEP["rs1"], "75")
        self.assertEqual(record, {
            "status": "genic", "genes": ["GENEA"], "gene_ids": ["ENSG00000000001"], "chr": "1", "pos": "1001",
            "consequences": ["5_prime_UTR_variant"], "biotypes": ["protein_coding"], "release": "75",
        })                                                       # upstream-only NEARBY gene excluded
        self.assertEqual(E.summarize_vep_entries([], "75")["status"], "not_in_ensembl")
        self.assertEqual(E.summarize_vep_entries(VEP["rs6"], "75")["status"], "no_primary_locus")
        self.assertEqual(E.summarize_vep_entries(VEP["rs8"], "75")["status"], "intergenic")  # downstream only

    def test_rows_share_one_schema(self):
        ucsc_row = C.annotation_row("rs5", U.build_record(None, [], {}, {}))
        vep_row = C.annotation_row("rs5", E.summarize_vep_entries([], "75"))
        legacy_row = C.annotation_row("rsOLD", LEGACY_VEP_RECORD)
        for row in (ucsc_row, vep_row, legacy_row):
            self.assertEqual(list(row), C.OUTPUT_COLUMNS)
        self.assertEqual((ucsc_row["GENES"], ucsc_row["CHR"], ucsc_row["POS"]), ("not in dbSNP", "NA", "NA"))
        self.assertEqual((vep_row["GENES"], vep_row["ENSEMBL_RELEASE"]), ("not in Ensembl", "75"))
        self.assertEqual((legacy_row["GENES"], legacy_row["CHR"], legacy_row["ENSEMBL_RELEASE"]), ("OLDGENE", "3", "116"))


# =====================================================================
# 4. CLI, server selection
# =====================================================================

class CliTests(unittest.TestCase):
    def test_source_argument(self):
        base = ["--inputs", "a.tsv", "--output", "o.tsv"]
        args = A.build_parser().parse_args(base)
        self.assertEqual((args.source, args.server, args.cache_file, args.workers), ("ucsc", None, None, None))
        self.assertEqual(A.build_parser().parse_args(base + ["--source", "ensembl"]).source, "ensembl")
        with mock.patch("sys.stderr"), self.assertRaises(SystemExit):
            A.build_parser().parse_args(base + ["--source", "ncbi"])

    def test_server_defaults_and_mismatch(self):
        self.assertEqual(A.resolve_server("ucsc"), U.DEFAULT_SERVER)
        self.assertEqual(A.resolve_server("ensembl"), E.DEFAULT_SERVER)
        self.assertEqual(A.resolve_server("ucsc", "https://api.genome.ucsc.edu/"), "https://api.genome.ucsc.edu")
        self.assertEqual(A.resolve_server("ensembl", "http://127.0.0.1:9"), "http://127.0.0.1:9")
        with self.assertRaisesRegex(ValueError, "ensembl API"):
            A.resolve_server("ucsc", "https://grch37.rest.ensembl.org")
        with self.assertRaisesRegex(ValueError, "ucsc API"):
            A.resolve_server("ensembl", "https://api.genome.ucsc.edu")
        with self.assertRaisesRegex(ValueError, "Unknown --source"):
            A.run_annotate_rna_genes([], "o.tsv", source="ncbi")

    def test_default_caches_are_separate(self):
        ucsc = C.default_cache_path(A.DEFAULT_CACHE_DIR, "ucsc")
        ensembl = C.default_cache_path(A.DEFAULT_CACHE_DIR, "ensembl")
        self.assertNotEqual(ucsc, ensembl)
        self.assertEqual((Path(ucsc).name, Path(ensembl).name),
                         ("ucsc_hg19_genes_cache.json", "ensembl_grch37_genes_cache.json"))


# =====================================================================
# 5. Per-backend caches (no mixing, no overwriting)
# =====================================================================

class CacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.ucsc = os.path.join(self.tmp.name, C.CACHE_NAMES["ucsc"])
        self.ensembl = os.path.join(self.tmp.name, C.CACHE_NAMES["ensembl"])

    def write(self, path, data):
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f)

    def read(self, path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)

    def test_record_source(self):
        self.assertEqual(C.record_source(UCSC_RECORD), "ucsc")
        self.assertEqual(C.record_source(LEGACY_VEP_RECORD), "ensembl")
        self.assertEqual(C.record_source(dict(LEGACY_VEP_RECORD, release="NA")), "ensembl")

    def test_mixed_legacy_cache_is_split(self):
        # The single cache used before both backends existed held VEP and UCSC records
        self.write(self.ensembl, {"rsOLD": LEGACY_VEP_RECORD, "rsU": UCSC_RECORD})
        with mock.patch("sys.stdout"):
            self.assertEqual(C.load_backend_cache(self.ucsc, "ucsc"), {"rsU": UCSC_RECORD})
        self.assertEqual(self.read(self.ucsc), {"rsU": UCSC_RECORD})
        self.assertEqual(self.read(self.ensembl), {"rsOLD": LEGACY_VEP_RECORD})
        before = (os.path.getmtime(self.ucsc), os.path.getmtime(self.ensembl))
        self.assertEqual(C.load_backend_cache(self.ensembl, "ensembl"), {"rsOLD": LEGACY_VEP_RECORD})
        self.assertEqual((os.path.getmtime(self.ucsc), os.path.getmtime(self.ensembl)), before)  # nothing rewritten

    def test_record_in_its_own_file_wins(self):
        stale = dict(UCSC_RECORD, genes=["STALE"])
        self.write(self.ucsc, {"rsX": UCSC_RECORD})
        self.write(self.ensembl, {"rsX": stale, "rsOLD": LEGACY_VEP_RECORD})
        with mock.patch("sys.stdout"):
            self.assertEqual(C.load_backend_cache(self.ensembl, "ensembl"), {"rsOLD": LEGACY_VEP_RECORD})
        self.assertEqual(self.read(self.ucsc), {"rsX": UCSC_RECORD})
        self.assertEqual(self.read(self.ensembl), {"rsOLD": LEGACY_VEP_RECORD})

    def test_no_files_created_when_nothing_to_move(self):
        self.assertEqual(C.load_backend_cache(self.ucsc, "ucsc"), {})
        self.assertFalse(os.path.exists(self.ucsc) or os.path.exists(self.ensembl))

    def test_one_shared_file_is_filtered_not_rewritten(self):
        shared = os.path.join(self.tmp.name, "shared.json")
        self.write(shared, {"rsOLD": LEGACY_VEP_RECORD, "rsU": UCSC_RECORD})
        with mock.patch.object(C, "CACHE_NAMES", {"ucsc": "shared.json", "ensembl": "shared.json"}):
            self.assertEqual(C.load_backend_cache(shared, "ucsc"), {"rsU": UCSC_RECORD})
            self.assertEqual(C.load_backend_cache(shared, "ensembl"), {"rsOLD": LEGACY_VEP_RECORD})
        self.assertEqual(len(self.read(shared)), 2)


# =====================================================================
# 6. Network: proxy bypass, retries/backoff, fail-fast checks, pacing (both clients)
# =====================================================================

class NetworkTests(MockServerCase):
    def test_sessions_ignore_proxy_environment(self):
        dead_proxy = f"http://127.0.0.1:{free_port()}"
        env = {"HTTP_PROXY": dead_proxy, "HTTPS_PROXY": dead_proxy, "http_proxy": dead_proxy,
               "https_proxy": dead_proxy, "NO_PROXY": "", "no_proxy": ""}
        with mock.patch.dict(os.environ, env):
            for source, client in self.clients().items():
                with self.subTest(source):
                    self.assertFalse(client.session.trust_env)
                    self.assertEqual(client.session.proxies, {})
                    client.check_server()                          # direct connection, proxy ignored
            entries, error = self.vep_client().post_vep(["rs2"], attempts=1)
            self.assertEqual((sorted(entries), error), (["rs2"], ""))
            # Control: a default session follows the dead proxy and fails
            with self.assertRaises(requests.exceptions.ConnectionError):
                requests.get(f"{self.url}/info/ping", timeout=5)

    def test_proxy_fallback_only_when_direct_connection_fails(self):
        unreachable = f"http://127.0.0.1:{free_port()}"            # nothing listens here
        env = {"HTTP_PROXY": self.url, "http_proxy": self.url, "NO_PROXY": "", "no_proxy": ""}
        with mock.patch.dict(os.environ, env):
            for source, client in self.clients(unreachable).items():
                with self.subTest(source):
                    client.check_server()                          # direct fails -> mock acts as proxy
                    self.assertEqual(client.session.proxies.get("http"), self.url)
                    self.assertFalse(client.session.trust_env)
        for source, direct in self.clients().items():
            with self.subTest(f"{source} direct"):
                direct.check_server()
                self.assertEqual(direct.session.proxies, {})       # direct works -> no proxy

    def test_unreachable_server_fails_fast(self):
        no_proxy = {"HTTP_PROXY": "", "http_proxy": "", "HTTPS_PROXY": "", "https_proxy": ""}
        with mock.patch.dict(os.environ, no_proxy):
            for source, client in self.clients(f"http://127.0.0.1:{free_port()}").items():
                with self.subTest(source):
                    with self.assertRaisesRegex(RuntimeError, "unavailable"):
                        client.check_server()
                    self.assertEqual(client.requests_sent, 3)

    def test_ensembl_outage_page_fails_fast(self):
        MockApi.ping_down = True
        client = self.vep_client()
        with self.assertRaisesRegex(RuntimeError, r"Ensembl REST server .* unavailable \(HTTP 500: Error: 500"):
            client.check_server()
        self.assertEqual(client.requests_sent, 3)

    def test_retry_after_and_backoff_on_retryable_statuses(self):
        MockApi.script.extend([
            (503, {"Retry-After": "0"}, '{"error": "busy"}', 0),
            (429, {"Retry-After": "0"}, '{"error": "slow down"}', 0),
            (502, {}, "<html><title>Bad gateway ’</title></html>", 0),
            (200, {}, '{"ok": 1}', 0),
        ])
        client = self.ucsc_client()
        payload, error = client.get_json("/flaky", {}, "flaky")
        self.assertEqual((payload, error, client.requests_sent), ({"ok": 1}, "", 4))

    def test_vep_post_is_retried(self):
        MockApi.vep_script.extend([(503, {"Retry-After": "0"}), (429, {"Retry-After": "0"})])
        client = self.vep_client()
        entries, error = client.post_vep(["rs1", "rs2"], attempts=6)
        self.assertEqual((sorted(entries), error, client.requests_sent), (["rs1", "rs2"], "", 3))

    def test_read_timeout_is_retried(self):
        MockApi.script.extend([(200, {}, '{"late": 1}', 1.0), (200, {}, '{"ok": 2}', 0)])
        client = U.UcscClient(self.url, min_interval=0, timeout=(2, 0.3))
        payload, _ = client.get_json("/flaky", {}, "flaky")
        self.assertEqual((payload, client.requests_sent), ({"ok": 2}, 2))

    def test_client_errors_are_not_retried(self):
        MockApi.script.extend([(400, {}, '{"error": "bad chrom"}', 0)])
        client = self.ucsc_client()
        payload, error = client.get_json("/flaky", {}, "flaky")
        self.assertIsNone(payload)
        self.assertIn("bad chrom", error)
        self.assertEqual(client.requests_sent, 1)

    def test_api_error_field_is_a_failure(self):
        MockApi.script.extend([(200, {}, '{"error": "track not found"}', 0)])
        payload, error = self.ucsc_client().get_json("/flaky", {}, "flaky")
        self.assertIsNone(payload)
        self.assertIn("track not found", error)

    def test_backoff_delay(self):
        with mock.patch.multiple(C, BACKOFF_BASE=2.0, BACKOFF_CAP=120.0):
            for attempt in range(8):
                self.assertLessEqual(C.retry_delay(attempt), min(120.0, 2.0 * 2 ** attempt))
            resp = requests.Response()
            resp.headers["Retry-After"] = "7"
            self.assertEqual(C.retry_delay(0, resp), 7.0)
            resp.headers["Retry-After"] = "999"
            self.assertEqual(C.retry_delay(0, resp), 120.0)

    def test_ucsc_requests_are_paced(self):
        client = U.UcscClient(self.url, min_interval=0.2)
        started = time.monotonic()
        for _ in range(3):
            client.get_json("/list/schema", {"track": "t"}, "schema")
        self.assertGreaterEqual(time.monotonic() - started, 0.38)


# =====================================================================
# 7. End to end against the mock APIs (both sources, caches, schema, failures)
# =====================================================================

class EndToEndTests(MockServerCase):
    def setUp(self):
        super().setUp()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = Path(self.tmp.name)
        self.caches = {s: str(d / C.CACHE_NAMES[s]) for s in C.SOURCES}
        self.output = str(d / "rna_sig_genes.tsv")
        self.positions = str(d / U.POSITIONS_CACHE_NAME)
        # A record annotated earlier by Ensembl VEP
        dump_json(self.caches["ensembl"], {"rsOLD": LEGACY_VEP_RECORD})
        rows = [("T1", s, "True") for s in ["rs1", "rs2", "rs3", "rs4", "rs5", "rsOLD"]] + \
               [("T2", s, "True") for s in ["rs6", "rs7", "rs8", "rs9", "rs10", "rs1"]] + [("T2", "rsNOTSIG", "False")]
        self.inputs = []
        for trait in ("T1", "T2"):
            path = d / f"{trait}_overlap_bh.tsv"
            pd.DataFrame([r for r in rows if r[0] == trait], columns=["TRAIT", "SNP", "BH_SIGNIFICANT"]).to_csv(path, sep="\t", index=False)
            self.inputs.append(str(path))

    def run_step(self, source="ucsc", client=None):
        client = client or (self.ucsc_client() if source == "ucsc" else self.vep_client())
        out = A.run_annotate_rna_genes(self.inputs, self.output, source=source, server=self.url,
                                       cache_path=self.caches[source], batch_size=2,
                                       attrs_url=f"{self.url}/attrs.txt.gz", client=client)
        return out.set_index("SNP"), client

    def read_cache(self, source):
        with open(self.caches[source], encoding="utf-8") as f:
            return json.load(f)

    def assert_rows(self, out, expected):
        for snp, values in expected.items():
            self.assertEqual(tuple(str(out.loc[snp, c]) for c in VALUE_COLUMNS), values, snp)

    def test_ucsc_full_annotation(self):
        out, client = self.run_step("ucsc")
        self.assertEqual(list(pd.read_csv(self.output, sep="\t").columns), C.OUTPUT_COLUMNS)
        self.assertNotIn("rsNOTSIG", out.index)
        self.assert_rows(out, EXPECTED_SHARED)
        # rsOLD's Ensembl record is not used: UCSC annotates it itself (unknown to the mock dbSNP)
        self.assert_rows(out, {"rs5": ("not_in_dbsnp", "not in dbSNP", "NA", "NA", "NA", "NA", "NA"),
                               "rsOLD": ("not_in_dbsnp", "not in dbSNP", "NA", "NA", "NA", "NA", "NA")})
        self.assertEqual(out.loc["rs4", "ENSEMBL_RELEASE"], "UCSC hg19 GENCODE V19 + snp151")
        self.assertTrue(out["ENSEMBL_RELEASE"].str.startswith("UCSC hg19 GENCODE V19").all())
        # 11 searches; comprehensive track per region (chr1 x 2, chr2), pseudogene track once per
        # chromosome (chr1, chr2); + server check + attributes download
        self.assertEqual(MockApi.hits.count("/search"), 11)
        gene_queries = sorted(h for h in MockApi.hits if h.startswith("getData:"))
        self.assertEqual(gene_queries, [
            "getData:wgEncodeGencodeCompV19:chr1:region", "getData:wgEncodeGencodeCompV19:chr1:region",
            "getData:wgEncodeGencodeCompV19:chr2:region",
            "getData:wgEncodeGencodePseudoGeneV19:chr1:chromosome", "getData:wgEncodeGencodePseudoGeneV19:chr2:chromosome",
        ])
        self.assertEqual(client.requests_sent, 18)
        self.assertEqual(self.read_cache("ensembl"), {"rsOLD": LEGACY_VEP_RECORD})   # untouched
        self.assertEqual({C.record_source(r) for r in self.read_cache("ucsc").values()}, {"ucsc"})

        MockApi.hits.clear()
        again, client2 = self.run_step("ucsc")
        self.assertEqual(client2.requests_sent, 0)                 # everything cached: no requests
        pd.testing.assert_frame_equal(again, out)

    def test_ensembl_full_annotation(self):
        out, client = self.run_step("ensembl")
        self.assertEqual(list(pd.read_csv(self.output, sep="\t").columns), C.OUTPUT_COLUMNS)
        self.assert_rows(out, EXPECTED_SHARED)
        self.assert_rows(out, {"rs5": ("not_in_ensembl", "not in Ensembl", "NA", "NA", "NA", "NA", "NA"),
                               "rsOLD": ("genic", "OLDGENE", "ENSG00000999999", "3", "12345", "intron_variant", "protein_coding")})
        self.assertEqual(out.loc["rsOLD", "ENSEMBL_RELEASE"], "116")               # from the cache
        self.assertEqual(out.loc["rs1", "ENSEMBL_RELEASE"], str(VEP_RELEASE))      # from /info/data
        # ping + release + 10 uncached SNPs in 5 batches of 2; rsOLD was not queried
        self.assertEqual(client.requests_sent, 7)
        self.assertNotIn("rsOLD", [i for batch in MockApi.vep_batches for i in batch])
        self.assertFalse(os.path.exists(self.caches["ucsc"]))
        self.assertEqual({C.record_source(r) for r in self.read_cache("ensembl").values()}, {"ensembl"})

        again, client2 = self.run_step("ensembl")
        self.assertEqual(client2.requests_sent, 0)
        pd.testing.assert_frame_equal(again, out)

    def test_schema_parity_and_switching_sources(self):
        ucsc_out, _ = self.run_step("ucsc")
        ucsc_cache = self.read_cache("ucsc")
        ensembl_out, _ = self.run_step("ensembl")
        self.assertEqual(list(ucsc_out.reset_index().columns), list(ensembl_out.reset_index().columns))
        self.assertEqual(list(ucsc_out.index), list(ensembl_out.index))
        self.assertEqual(dict(ucsc_out.dtypes), dict(ensembl_out.dtypes))
        shared = list(EXPECTED_SHARED)
        pd.testing.assert_frame_equal(ucsc_out.loc[shared, VALUE_COLUMNS], ensembl_out.loc[shared, VALUE_COLUMNS])
        # Switching never touches the other backend's cache, and switching back needs no requests
        self.assertEqual(self.read_cache("ucsc"), ucsc_cache)
        self.assertEqual({C.record_source(r) for r in self.read_cache("ensembl").values()}, {"ensembl"})
        back, client = self.run_step("ucsc")
        self.assertEqual(client.requests_sent, 0)
        pd.testing.assert_frame_equal(back, ucsc_out)

    def test_mixed_cache_is_split_before_a_run(self):
        # The pre-split single cache (Ensembl file) also held UCSC annotations of rs1 and rs2
        ucsc_rs1 = U.build_record("dbSnp155", [("chr1", 1001, 1001)], {}, {})
        dump_json(self.caches["ensembl"], {"rsOLD": LEGACY_VEP_RECORD, "rs1": ucsc_rs1, "rs2": UCSC_RECORD})
        out, _ = self.run_step("ucsc")
        self.assertEqual(MockApi.hits.count("/search"), 9)          # rs1, rs2 reused, not searched
        self.assertEqual(out.loc["rs2", "GENES"], "UGENE")
        self.assertEqual(self.read_cache("ensembl"), {"rsOLD": LEGACY_VEP_RECORD})
        ensembl_out, _ = self.run_step("ensembl")                   # rs1, rs2 queried from VEP
        self.assertEqual(ensembl_out.loc["rs2", "GENES"], "GENEA")

    def test_ucsc_failed_regions_are_not_cached_and_retried(self):
        MockApi.fail_chroms.add("chr2")
        with self.assertRaisesRegex(RuntimeError, r"UCSC lookups failed for 3 SNPs .*'rs7', 'rs8', 'rs9'"):
            self.run_step("ucsc")
        cache = self.read_cache("ucsc")
        self.assertTrue({"rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rsOLD"} <= set(cache))
        self.assertFalse({"rs7", "rs8", "rs9"} & set(cache))
        self.assertEqual(len(load_json(self.positions)), 11)    # positions kept for the retry

        MockApi.fail_chroms.clear()
        MockApi.hits.clear()
        out, _ = self.run_step("ucsc")
        self.assertEqual(MockApi.hits.count("/search"), 0)          # no repeated searches
        self.assertEqual(sorted({h for h in MockApi.hits if h.startswith("/")} - {"/list/schema"}), ["/getData/track"])
        self.assertEqual(out.loc["rs9", "GENES"], "GENEA, PSEUDO1")

    def test_ensembl_failing_snp_is_isolated_and_retried(self):
        MockApi.vep_fail_ids.add("rs3")
        with self.assertRaisesRegex(RuntimeError, r"Ensembl VEP lookups failed for 1 SNPs .*'rs3'"):
            self.run_step("ensembl")
        cache = self.read_cache("ensembl")
        self.assertNotIn("rs3", cache)
        self.assertTrue({"rs1", "rs2", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9", "rs10", "rsOLD"} <= set(cache))

        MockApi.reset()
        out, _ = self.run_step("ensembl")
        self.assertEqual(MockApi.vep_batches, [["rs3"]])            # only the failed SNP
        self.assertEqual(out.loc["rs3", "GENES"], "GENEB")

    def test_ensembl_outage_queries_nothing(self):
        MockApi.ping_down = True
        with self.assertRaisesRegex(RuntimeError, "unavailable"):
            self.run_step("ensembl")
        self.assertEqual(MockApi.vep_batches, [])
        self.assertEqual(self.read_cache("ensembl"), {"rsOLD": LEGACY_VEP_RECORD})
        self.assertFalse(os.path.exists(self.output))

    def test_ucsc_progress_is_flushed_when_interrupted(self):
        client = self.ucsc_client()
        original = client.search_rsid
        calls = {"n": 0}

        def crash_after_three(rsid):
            calls["n"] += 1
            if calls["n"] > 3:
                raise KeyboardInterrupt
            return original(rsid)

        client.search_rsid = crash_after_three
        with self.assertRaises(KeyboardInterrupt):
            self.run_step("ucsc", client)
        self.assertEqual(len(load_json(self.positions)), 3)

    def test_cli_runs_the_ensembl_backend(self):
        A.main(["--inputs", *self.inputs, "--output", self.output, "--source", "ensembl",
                "--server", self.url, "--cache-file", self.caches["ensembl"], "--workers", "1"])
        out = pd.read_csv(self.output, sep="\t", dtype=str, keep_default_na=False).set_index("SNP")
        self.assert_rows(out, EXPECTED_SHARED)
        self.assertEqual(MockApi.hits.count("/info/ping"), 1)


# =====================================================================
# 8. Live smoke tests (RUN_LIVE_UCSC=1 / RUN_LIVE_ENSEMBL=1)
# =====================================================================

@unittest.skipUnless(os.environ.get("RUN_LIVE_UCSC") == "1", "set RUN_LIVE_UCSC=1 to query api.genome.ucsc.edu")
class LiveUcscTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.client = U.UcscClient(U.DEFAULT_SERVER)
        cls.client.check_server()

    def test_tp53_region(self):
        start0, end0 = U.ensembl_to_ucsc(7565096, 7590856)
        txs = self.client.transcripts("wgEncodeGencodeCompV19", "chr17", start0, end0)
        self.assertIsNotNone(txs)
        self.assertTrue({"TP53", "WRAP53"} <= {t["gene_symbol"] for t in txs})
        self.assertTrue(all(U.intervals_overlap(start0, end0, t["tx_start"], t["tx_end"]) for t in txs))
        tp53 = [t for t in txs if t["gene_symbol"] == "TP53"]
        self.assertEqual(min(t["tx_start"] for t in tp53), 7565096)  # GENCODE V19 TP53 starts at 0-based 7565096

    def test_rs1002189(self):
        track, loci = self.client.search_rsid("rs1002189")
        self.assertEqual((track, loci), ("dbSnp155", [("chr22", 30771458, 30771458)]))
        with tempfile.TemporaryDirectory() as d:
            attrs = U.load_gencode_attrs(d, self.client)
            cache = {}
            failed = U.annotate_genes(self.client, ["rs1002189"], {"rs1002189": {"track": track, "loci": loci}},
                                      attrs, cache, os.path.join(d, "cache.json"), 25)
        self.assertEqual(failed, [])
        record = cache["rs1002189"]
        print(f"\n    UCSC rs1002189 -> {C.annotation_row('rs1002189', record)}")
        self.assertEqual((record["chr"], record["pos"]), ("22", "30771458"))
        self.assertIn(record["status"], {"genic", "intergenic"})


@unittest.skipUnless(os.environ.get("RUN_LIVE_ENSEMBL") == "1", "set RUN_LIVE_ENSEMBL=1 to query grch37.rest.ensembl.org")
class LiveEnsemblTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.client = E.VepClient(E.DEFAULT_SERVER)
        try:
            cls.client.check_server()
        except RuntimeError as exc:
            raise unittest.SkipTest(f"Ensembl GRCh37 server unavailable: {exc}")

    def test_rs1002189(self):
        release = self.client.fetch_release()
        entries, error = self.client.post_vep(["rs1002189"], attempts=3)
        self.assertEqual(error, "")
        record = E.summarize_vep_entries(entries.get("rs1002189", []), release)
        print(f"\n    Ensembl rs1002189 -> {C.annotation_row('rs1002189', record)}")
        self.assertEqual((record["chr"], record["pos"]), ("22", "30771458"))   # same hg19 position as UCSC
        self.assertIn(record["status"], {"genic", "intergenic"})


if __name__ == "__main__":
    unittest.main()
