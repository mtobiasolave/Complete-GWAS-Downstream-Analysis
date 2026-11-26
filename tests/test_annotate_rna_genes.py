"""Tests for scripts/annotate_rna_genes.py (UCSC hg19 / GENCODE V19 gene annotation, Step B2).

Run from the repository root with snakemake9_env activated:
    python -m unittest discover -s tests -v
Live smoke tests against https://api.genome.ucsc.edu are skipped unless RUN_LIVE_UCSC=1.
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

FAST = {"BACKOFF_BASE": 0.01, "BACKOFF_CAP": 0.05}


# =====================================================================
# Fixtures: synthetic GENCODE V19 tracks, dbSNP search results, attributes
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
LEGACY_VEP_RECORD = {
    "status": "genic", "genes": ["OLDGENE"], "gene_ids": ["ENSG00000999999"], "chr": "3", "pos": "12345",
    "consequences": ["intron_variant"], "biotypes": ["protein_coding"], "release": "116",
}


class MockUcsc(BaseHTTPRequestHandler):
    """Minimal UCSC API: /list/schema, /search, /getData/track, /attrs.txt.gz, scripted /flaky."""
    hits = []
    fail_chroms = set()
    script = deque()

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
        MockUcsc.hits.append(url.path)
        if url.path == "/flaky":
            code, headers, body, delay = MockUcsc.script.popleft()
            time.sleep(delay)
            return self.reply(code, body, headers=headers)
        if url.path == "/list/schema":
            return self.reply(200, json.dumps({"track": q.get("track"), "columnTypes": []}))
        if url.path == "/search":
            return self.reply(200, json.dumps({"genome": "hg19", "positionMatches": SEARCH.get(q["search"], [])}))
        if url.path == "/getData/track":
            if q["chrom"] in MockUcsc.fail_chroms:
                return self.reply(503, json.dumps({"error": "busy"}), headers={"Retry-After": "0"})
            MockUcsc.hits.append(f"getData:{q['track']}:{q['chrom']}:{'region' if 'start' in q else 'chromosome'}")
            start, end = int(q.get("start", 0)), int(q.get("end", 10**12))
            items = [i for i in TRACK_ITEMS.get(q["track"], [])
                     if i["chrom"] == q["chrom"] and i["txStart"] < end and i["txEnd"] > start]
            return self.reply(200, json.dumps({"track": q["track"], q["track"]: items}))
        if url.path == "/attrs.txt.gz":
            return self.reply(200, gzip.compress(("\n".join(ATTRS_LINES) + "\n").encode()), "application/gzip")
        return self.reply(400, json.dumps({"error": f"unknown endpoint {url.path}"}))


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
        cls.httpd = QuietServer(("127.0.0.1", 0), MockUcsc)
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()
        cls.url = f"http://127.0.0.1:{cls.httpd.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def setUp(self):
        MockUcsc.hits.clear()
        MockUcsc.fail_chroms.clear()
        MockUcsc.script.clear()
        patcher = mock.patch.multiple(A, **FAST)
        patcher.start()
        self.addCleanup(patcher.stop)

    def client(self, url=None):
        return A.UcscClient(url or self.url, min_interval=0)


# =====================================================================
# 1. Chromosome names
# =====================================================================

class ChromosomeTests(unittest.TestCase):
    def test_to_ucsc(self):
        for bare, ucsc in [("1", "chr1"), ("22", "chr22"), ("X", "chrX"), ("x", "chrX"), ("Y", "chrY"),
                           ("MT", "chrM"), ("M", "chrM"), ("chr7", "chr7"), ("chrM", "chrM"), (17, "chr17")]:
            self.assertEqual(A.to_ucsc_chrom(bare), ucsc, bare)

    def test_from_ucsc(self):
        for ucsc, bare in [("chr1", "1"), ("chr22", "22"), ("chrX", "X"), ("chrM", "MT"), ("chrMT", "MT"), ("5", "5")]:
            self.assertEqual(A.from_ucsc_chrom(ucsc), bare, ucsc)

    def test_round_trip(self):
        for bare in sorted(A.PRIMARY_CHROMS):
            self.assertEqual(A.from_ucsc_chrom(A.to_ucsc_chrom(bare)), bare)

    def test_primary(self):
        self.assertTrue(A.is_primary("chr1") and A.is_primary("chrM") and A.is_primary("MT"))
        for name in ("chr6_cox_hap2", "chr17_gl000205_random", "chrUn_gl000220", "GL000192.1"):
            self.assertFalse(A.is_primary(name), name)


# =====================================================================
# 2. Coordinates
# =====================================================================

class CoordinateTests(unittest.TestCase):
    def test_conversion(self):
        self.assertEqual(A.ensembl_to_ucsc(7565096, 7590856), (7565095, 7590856))
        self.assertEqual(A.ucsc_to_ensembl(7565095, 7590856), (7565096, 7590856))
        self.assertEqual(A.ucsc_to_ensembl(*A.ensembl_to_ucsc(100, 100)), (100, 100))
        self.assertEqual(A.ensembl_to_ucsc(30771458, 30771458), (30771457, 30771458))  # rs1002189

    def test_single_base_interval(self):
        self.assertEqual(A.locus_query_interval(1001, 1001), (1000, 1001))
        self.assertEqual(A.locus_query_interval(1001, 1000), (1000, 1001))   # insertion: end < start

    def test_transcript_boundaries(self):
        # UCSC transcript [100, 200) covers 1-based bases 101..200
        tx = A.parse_gene_track({"t": [gene_item("ENST1.1", "G", "chr1", "+", (100, 200), (200, 200), [(100, 200)])]}, "t")[0]
        inside = {pos: bool(A.transcript_consequences(tx, *A.locus_query_interval(pos, pos)))
                  for pos in (100, 101, 200, 201)}
        self.assertEqual(inside, {100: False, 101: True, 200: True, 201: False})

    def test_clusters(self):
        loci = [("chr1", 1001, 1001), ("chr1", 50000, 50000), ("chr1", 400000, 400000), ("chr2", 10, 10)]
        clusters = A.cluster_loci(loci, max_gap=100_000, max_span=2_000_000)
        self.assertEqual([(c[0], c[1], c[2], len(c[3])) for c in clusters],
                         [("chr1", 1000, 50000, 2), ("chr1", 399999, 400000, 1), ("chr2", 9, 10, 1)])
        self.assertEqual(len(A.cluster_loci(loci[:2], max_gap=100_000, max_span=10_000)), 2)


# =====================================================================
# 3. UCSC JSON parsing & record mapping
# =====================================================================

class ParsingTests(unittest.TestCase):
    def test_highlight(self):
        self.assertEqual(A.parse_highlight("hg19.chr22:30771458-30771458#fcfcac"), ("chr22", 30771458, 30771458))
        self.assertEqual(A.parse_highlight("chr6_cox_hap2:2377304-2377305"), ("chr6_cox_hap2", 2377304, 2377305))
        self.assertIsNone(A.parse_highlight(""))
        self.assertIsNone(A.parse_highlight(None))

    def test_search_prefers_newest_dbsnp_and_exact_name(self):
        payload = {"positionMatches": [
            {"trackName": "gtexEqtlCluster", "matches": [{"posName": "rs1002189", "highlight": "hg19.chr22:1-1"}]},
            snp_match("snp151", "rs1002189", "chr22", 30771458),
            snp_match("dbSnp155Common", "rs1002189", "chr22", 999),        # subset track: ignored
            snp_match("dbSnp155", "rs1002189", "chr22", 30771458),
            snp_match("dbSnp155", "rs10021890", "chr4", 5),               # other rsID: ignored
        ]}
        self.assertEqual(A.parse_search_loci(payload, "rs1002189"), ("dbSnp155", [("chr22", 30771458, 30771458)]))
        self.assertEqual(A.parse_search_loci({"positionMatches": [snp_match("snp151", "rsX", "chr1", 5)]}, "RSX"),
                         ("snp151", [("chr1", 5, 5)]))
        self.assertEqual(A.parse_search_loci({"genome": "hg19"}, "rs1"), (None, []))

    def test_search_keeps_all_loci(self):
        payload = {"positionMatches": [snp_match("dbSnp155", "rs9", "chr1", 7), snp_match("dbSnp155", "rs9", "chr6_cox_hap2", 3)]}
        self.assertEqual(A.parse_search_loci(payload, "rs9")[1], [("chr1", 7, 7), ("chr6_cox_hap2", 3, 3)])

    def test_gene_track_fields(self):
        txs = A.parse_gene_track({"wgEncodeGencodeCompV19": TRACK_ITEMS["wgEncodeGencodeCompV19"]}, "wgEncodeGencodeCompV19")
        self.assertEqual(txs[0]["transcript_id"], "ENST00000000001.1")
        self.assertEqual(txs[0]["gene_symbol"], "GENEA")
        self.assertEqual((txs[0]["tx_start"], txs[0]["tx_end"], txs[0]["cds_start"], txs[0]["cds_end"]), (1000, 5000, 1500, 4000))
        self.assertEqual(txs[0]["exons"], [(1000, 1600), (3000, 5000)])
        self.assertEqual(A.parse_gene_track({"other": []}, "wgEncodeGencodeCompV19"), [])

    def test_consequences(self):
        coding, noncoding = A.parse_gene_track({"t": TRACK_ITEMS["wgEncodeGencodeCompV19"]}, "t")
        cons = lambda tx, pos: A.transcript_consequences(tx, *A.locus_query_interval(pos, pos))
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
        attrs = A.parse_gencode_attrs(line + "\n" for line in ATTRS_LINES)
        self.assertEqual(attrs["ENST00000000001.1"]["gene_id"], "ENSG00000000001")
        self.assertIs(attrs["ENST00000000001"], attrs["ENST00000000001.1"])
        tx = A.parse_gene_track({"t": TRACK_ITEMS["wgEncodeGencodeCompV19"]}, "t")[0]
        locus = ("chr1", 1001, 1001)
        record = A.build_record("dbSnp155", [locus, ("chr6_cox_hap2", 1, 1)], {locus: [(tx, ["5_prime_UTR_variant"])]}, attrs)
        self.assertEqual(record, {
            "status": "genic", "genes": ["GENEA"], "gene_ids": ["ENSG00000000001"], "chr": "1", "pos": "1001",
            "consequences": ["5_prime_UTR_variant"], "biotypes": ["protein_coding"],
            "release": "UCSC hg19 GENCODE V19 + dbSnp155",
        })
        self.assertEqual(A.build_record(None, [], {}, attrs)["status"], "not_in_dbsnp")
        self.assertEqual(A.build_record("dbSnp155", [("chr6_cox_hap2", 1, 1)], {}, attrs)["status"], "no_primary_locus")
        self.assertEqual(A.build_record("dbSnp155", [locus], {locus: []}, attrs)["status"], "intergenic")

    def test_rows_keep_schema_for_new_and_legacy_records(self):
        legacy = A.annotation_row("rsOLD", LEGACY_VEP_RECORD)
        self.assertEqual(list(legacy), A.OUTPUT_COLUMNS)
        self.assertEqual((legacy["GENES"], legacy["CHR"], legacy["ENSEMBL_RELEASE"]), ("OLDGENE", "3", "116"))
        missing = A.annotation_row("rs5", A.build_record(None, [], {}, {}))
        self.assertEqual((missing["GENES"], missing["CHR"], missing["POS"]), ("not in dbSNP", "NA", "NA"))
        legacy_missing = A.annotation_row("rsZ", dict(LEGACY_VEP_RECORD, status="not_in_ensembl", genes=[], gene_ids=[]))
        self.assertEqual(legacy_missing["GENES"], "not in Ensembl")

    def test_legacy_server_maps_to_ucsc(self):
        self.assertEqual(A.resolve_server("https://grch37.rest.ensembl.org"), A.DEFAULT_SERVER)
        self.assertEqual(A.resolve_server("https://api.genome.ucsc.edu/"), "https://api.genome.ucsc.edu")


# =====================================================================
# 4. Network: proxy bypass, retries/backoff, pacing
# =====================================================================

class NetworkTests(MockServerCase):
    def test_session_ignores_proxy_environment(self):
        dead_proxy = f"http://127.0.0.1:{free_port()}"
        env = {"HTTP_PROXY": dead_proxy, "HTTPS_PROXY": dead_proxy, "http_proxy": dead_proxy,
               "https_proxy": dead_proxy, "NO_PROXY": "", "no_proxy": ""}
        with mock.patch.dict(os.environ, env):
            client = self.client()
            self.assertFalse(client.session.trust_env)
            self.assertEqual(client.session.proxies, {})
            payload, error = client.get_json("/list/schema", {"genome": "hg19", "track": "t"}, "schema")
            self.assertEqual((payload["track"], error), ("t", ""))
            # Control: a default session follows the dead proxy and fails
            with self.assertRaises(requests.exceptions.ConnectionError):
                requests.get(f"{self.url}/list/schema", timeout=5)

    def test_proxy_fallback_only_when_direct_connection_fails(self):
        unreachable = f"http://127.0.0.1:{free_port()}"            # nothing listens here
        env = {"HTTP_PROXY": self.url, "http_proxy": self.url, "NO_PROXY": "", "no_proxy": ""}
        with mock.patch.dict(os.environ, env):
            client = self.client(unreachable)
            client.check_server()                                  # direct fails -> mock acts as proxy
            self.assertEqual(client.session.proxies.get("http"), self.url)
            self.assertFalse(client.session.trust_env)
        direct = self.client()
        direct.check_server()
        self.assertEqual(direct.session.proxies, {})               # direct works -> no proxy

    def test_unreachable_server_fails_fast(self):
        with mock.patch.dict(os.environ, {"HTTP_PROXY": "", "http_proxy": "", "HTTPS_PROXY": "", "https_proxy": ""}):
            client = self.client(f"http://127.0.0.1:{free_port()}")
            with self.assertRaisesRegex(RuntimeError, "unavailable"):
                client.check_server()
            self.assertEqual(client.requests_sent, 3)

    def test_retry_after_and_backoff_on_retryable_statuses(self):
        MockUcsc.script.extend([
            (503, {"Retry-After": "0"}, '{"error": "busy"}', 0),
            (429, {"Retry-After": "0"}, '{"error": "slow down"}', 0),
            (502, {}, "<html><title>Bad gateway ’</title></html>", 0),
            (200, {}, '{"ok": 1}', 0),
        ])
        client = self.client()
        payload, error = client.get_json("/flaky", {}, "flaky")
        self.assertEqual((payload, error, client.requests_sent), ({"ok": 1}, "", 4))

    def test_read_timeout_is_retried(self):
        MockUcsc.script.extend([(200, {}, '{"late": 1}', 1.0), (200, {}, '{"ok": 2}', 0)])
        with mock.patch.object(A, "TIMEOUT", (2, 0.3)):
            client = self.client()
            payload, _ = client.get_json("/flaky", {}, "flaky")
        self.assertEqual((payload, client.requests_sent), ({"ok": 2}, 2))

    def test_client_errors_are_not_retried(self):
        MockUcsc.script.extend([(400, {}, '{"error": "bad chrom"}', 0)])
        client = self.client()
        payload, error = client.get_json("/flaky", {}, "flaky")
        self.assertIsNone(payload)
        self.assertIn("bad chrom", error)
        self.assertEqual(client.requests_sent, 1)

    def test_api_error_field_is_a_failure(self):
        MockUcsc.script.extend([(200, {}, '{"error": "track not found"}', 0)])
        payload, error = self.client().get_json("/flaky", {}, "flaky")
        self.assertIsNone(payload)
        self.assertIn("track not found", error)

    def test_backoff_delay(self):
        with mock.patch.multiple(A, BACKOFF_BASE=2.0, BACKOFF_CAP=120.0):
            for attempt in range(8):
                self.assertLessEqual(A.retry_delay(attempt), min(120.0, 2.0 * 2 ** attempt))
            resp = requests.Response()
            resp.headers["Retry-After"] = "7"
            self.assertEqual(A.retry_delay(0, resp), 7.0)
            resp.headers["Retry-After"] = "999"
            self.assertEqual(A.retry_delay(0, resp), 120.0)

    def test_requests_are_paced(self):
        client = A.UcscClient(self.url, min_interval=0.2)
        started = time.monotonic()
        for _ in range(3):
            client.get_json("/list/schema", {"track": "t"}, "schema")
        self.assertGreaterEqual(time.monotonic() - started, 0.38)


# =====================================================================
# 5. End to end against the mock API (caches, schema, failures)
# =====================================================================

class EndToEndTests(MockServerCase):
    def setUp(self):
        super().setUp()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = Path(self.tmp.name)
        self.cache = str(d / "genes_cache.json")
        self.output = str(d / "rna_sig_genes.tsv")
        self.positions = str(d / A.POSITIONS_CACHE_NAME)
        json.dump({"rsOLD": LEGACY_VEP_RECORD}, open(self.cache, "w"))
        rows = [("T1", s, "True") for s in ["rs1", "rs2", "rs3", "rs4", "rs5", "rsOLD"]] + \
               [("T2", s, "True") for s in ["rs6", "rs7", "rs8", "rs9", "rs10", "rs1"]] + [("T2", "rsNOTSIG", "False")]
        self.inputs = []
        for trait in ("T1", "T2"):
            path = d / f"{trait}_overlap_bh.tsv"
            pd.DataFrame([r for r in rows if r[0] == trait], columns=["TRAIT", "SNP", "BH_SIGNIFICANT"]).to_csv(path, sep="\t", index=False)
            self.inputs.append(str(path))

    def run_step(self, client=None):
        client = client or self.client()
        out = A.run_annotate_rna_genes(self.inputs, self.output, server=self.url, cache_path=self.cache,
                                       batch_size=2, attrs_url=f"{self.url}/attrs.txt.gz", client=client)
        return out.set_index("SNP"), client

    def test_full_annotation(self):
        out, client = self.run_step()
        self.assertEqual(list(pd.read_csv(self.output, sep="\t").columns), A.OUTPUT_COLUMNS)
        self.assertNotIn("rsNOTSIG", out.index)
        expect = {
            "rs1": ("genic", "GENEA", "ENSG00000000001", "1", "1001", "5_prime_UTR_variant", "protein_coding"),
            "rs2": ("genic", "GENEA", "ENSG00000000001", "1", "2000", "intron_variant", "protein_coding"),
            "rs3": ("genic", "GENEB", "ENSG00000000002", "1", "20950", "non_coding_transcript_exon_variant", "lincRNA"),
            "rs4": ("intergenic", "intergenic", "NA", "1", "10000", "NA", "NA"),
            "rs5": ("not_in_dbsnp", "not in dbSNP", "NA", "NA", "NA", "NA", "NA"),
            "rs6": ("no_primary_locus", "no primary-assembly locus", "NA", "NA", "NA", "NA", "NA"),
            "rs7": ("genic", "PSEUDO1", "ENSG00000000003", "2", "800", "non_coding_transcript_exon_variant", "processed_pseudogene"),
            "rs8": ("intergenic", "intergenic", "NA", "2", "801", "NA", "NA"),
            "rs9": ("genic", "GENEA, PSEUDO1", "ENSG00000000001, ENSG00000000003", "1, 2", "1001, 800",
                    "5_prime_UTR_variant, non_coding_transcript_exon_variant", "processed_pseudogene, protein_coding"),
            "rs10": ("intergenic", "intergenic", "NA", "1", "900000", "NA", "NA"),
            "rsOLD": ("genic", "OLDGENE", "ENSG00000999999", "3", "12345", "intron_variant", "protein_coding"),
        }
        cols = ["STATUS", "GENES", "ENSEMBL_GENE_IDS", "CHR", "POS", "CONSEQUENCES", "TRANSCRIPT_BIOTYPES"]
        for snp, values in expect.items():
            self.assertEqual(tuple(str(out.loc[snp, c]) for c in cols), values, snp)
        self.assertEqual(out.loc["rs4", "ENSEMBL_RELEASE"], "UCSC hg19 GENCODE V19 + snp151")
        self.assertEqual(out.loc["rsOLD", "ENSEMBL_RELEASE"], "116")
        # rsOLD came from the cache. 10 searches; comprehensive track per region (chr1 x 2, chr2),
        # pseudogene track once per chromosome (chr1, chr2); + schema check + attributes download
        self.assertNotIn("rsOLD", json.load(open(self.positions)))
        self.assertEqual(MockUcsc.hits.count("/search"), 10)
        gene_queries = sorted(h for h in MockUcsc.hits if h.startswith("getData:"))
        self.assertEqual(gene_queries, [
            "getData:wgEncodeGencodeCompV19:chr1:region", "getData:wgEncodeGencodeCompV19:chr1:region",
            "getData:wgEncodeGencodeCompV19:chr2:region",
            "getData:wgEncodeGencodePseudoGeneV19:chr1:chromosome", "getData:wgEncodeGencodePseudoGeneV19:chr2:chromosome",
        ])
        self.assertEqual(client.requests_sent, 17)

        MockUcsc.hits.clear()
        again, client2 = self.run_step()
        self.assertEqual(client2.requests_sent, 0)                 # everything cached: no requests
        pd.testing.assert_frame_equal(again, out)

    def test_failed_regions_are_not_cached_and_retried(self):
        MockUcsc.fail_chroms.add("chr2")
        with self.assertRaisesRegex(RuntimeError, r"failed for 3 SNPs .*'rs7', 'rs8', 'rs9'"):
            self.run_step()
        cache = json.load(open(self.cache))
        self.assertTrue({"rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rsOLD"} <= set(cache))
        self.assertFalse({"rs7", "rs8", "rs9"} & set(cache))
        self.assertEqual(len(json.load(open(self.positions))), 10)    # positions kept for the retry

        MockUcsc.fail_chroms.clear()
        MockUcsc.hits.clear()
        out, client = self.run_step()
        self.assertEqual(MockUcsc.hits.count("/search"), 0)          # no repeated searches
        self.assertEqual(sorted({h for h in MockUcsc.hits if h.startswith("/")} - {"/list/schema"}), ["/getData/track"])
        self.assertEqual(out.loc["rs9", "GENES"], "GENEA, PSEUDO1")

    def test_progress_is_flushed_when_interrupted(self):
        client = self.client()
        original = client.search_rsid
        calls = {"n": 0}

        def crash_after_three(rsid):
            calls["n"] += 1
            if calls["n"] > 3:
                raise KeyboardInterrupt
            return original(rsid)

        client.search_rsid = crash_after_three
        with self.assertRaises(KeyboardInterrupt):
            self.run_step(client)
        self.assertEqual(len(json.load(open(self.positions))), 3)


# =====================================================================
# 6. Live smoke tests (RUN_LIVE_UCSC=1)
# =====================================================================

@unittest.skipUnless(os.environ.get("RUN_LIVE_UCSC") == "1", "set RUN_LIVE_UCSC=1 to query api.genome.ucsc.edu")
class LiveUcscTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.client = A.UcscClient(A.DEFAULT_SERVER)
        cls.client.check_server()

    def test_tp53_region(self):
        start0, end0 = A.ensembl_to_ucsc(7565096, 7590856)
        txs = self.client.transcripts("wgEncodeGencodeCompV19", "chr17", start0, end0)
        self.assertIsNotNone(txs)
        self.assertTrue({"TP53", "WRAP53"} <= {t["gene_symbol"] for t in txs})
        self.assertTrue(all(A.intervals_overlap(start0, end0, t["tx_start"], t["tx_end"]) for t in txs))
        tp53 = [t for t in txs if t["gene_symbol"] == "TP53"]
        self.assertEqual(min(t["tx_start"] for t in tp53), 7565096)  # GENCODE V19 TP53 starts at 0-based 7565096

    def test_rs1002189(self):
        track, loci = self.client.search_rsid("rs1002189")
        self.assertEqual((track, loci), ("dbSnp155", [("chr22", 30771458, 30771458)]))
        with tempfile.TemporaryDirectory() as d:
            attrs = A.load_gencode_attrs(d, self.client)
            cache = {}
            failed = A.annotate_genes(self.client, ["rs1002189"], {"rs1002189": {"track": track, "loci": loci}},
                                      attrs, cache, os.path.join(d, "cache.json"), 25)
        self.assertEqual(failed, [])
        record = cache["rs1002189"]
        print(f"\n    rs1002189 -> {A.annotation_row('rs1002189', record)}")
        self.assertEqual((record["chr"], record["pos"]), ("22", "30771458"))
        self.assertIn(record["status"], {"genic", "intergenic"})


if __name__ == "__main__":
    unittest.main()
