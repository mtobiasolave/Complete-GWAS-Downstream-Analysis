#!/usr/bin/env python3
"""
Step B2: Positional Gene Annotation of BH-Significant RNA m6A SNPs (hg19 / GRCh37).

1. Collects every SNP that is BH-significant in at least one trait (Step B1 outputs, rsIDs only).
2. Annotates each uncached SNP with the genes whose transcripts it physically falls inside
   (exon or intron, all biotypes, no upstream/downstream window), using one of two backends:
     --source ucsc     UCSC Genome Browser REST API, hg19, GENCODE V19 tracks, rsID positions
                       from UCSC's dbSNP tracks (rna_genes_ucsc.py; default)
     --source ensembl  Ensembl VEP REST API on the GRCh37 server (rna_genes_ensembl.py)
   Both share the hardened HTTP client of rna_genes_common.py (proxy environment ignored,
   explicit timeouts, backoff with jitter + Retry-After, fail-fast server check).
3. Each backend has its own persistent cache (results/.cache/ucsc_hg19_genes_cache.json,
   results/.cache/ensembl_grch37_genes_cache.json), flushed during the run. Records of the
   other backend found in a cache are moved to that backend's file, so annotations from the two
   databases are never mixed. Failed lookups are never cached; they make the job fail so a
   re-run retries only them.

Output columns (identical for both backends):
    SNP, GENES, ENSEMBL_GENE_IDS, STATUS, CHR, POS (1-based, chromosomes 1..22, X, Y, MT),
    CONSEQUENCES, TRANSCRIPT_BIOTYPES,
    ENSEMBL_RELEASE (annotation source: "UCSC hg19 GENCODE V19 + <dbSNP track>" or the
    Ensembl release number)
STATUS: genic | intergenic | no_primary_locus | not_in_dbsnp (UCSC) | not_in_ensembl (Ensembl).

Usage:
    python scripts/annotate_rna_genes.py --inputs results/rna/*_overlap_bh.tsv \
        --output results/rna/rna_sig_genes.tsv [--source ucsc|ensembl]
"""

import argparse
import os
from urllib.parse import urlparse

import pandas as pd

import rna_genes_ensembl
import rna_genes_ucsc
from common_utils import PROJECT_ROOT
from pleiotropy import is_true
from rna_genes_common import (
    DEFAULT_SOURCE,
    LOG,
    OUTPUT_COLUMNS,
    SOURCES,
    annotation_row,
    default_cache_path,
    load_backend_cache,
)

BACKENDS = {"ucsc": rna_genes_ucsc, "ensembl": rna_genes_ensembl}
DEFAULT_CACHE_DIR = PROJECT_ROOT / "results" / ".cache"
DEFAULT_BATCH_SIZE = 25
# Host suffix of each backend's API, to catch a server configured for the other backend
SERVER_HOSTS = {"ucsc": "genome.ucsc.edu", "ensembl": "ensembl.org"}


def significant_snps(input_paths: list) -> list:
    """Unique SNPs BH-significant in at least one Step B1 overlap table."""
    snps = set()
    for path in input_paths:
        df = pd.read_csv(path, sep="\t", usecols=["SNP", "BH_SIGNIFICANT"], dtype=str)
        snps.update(df.loc[is_true(df["BH_SIGNIFICANT"]), "SNP"].str.strip())
    return sorted(snps)


def resolve_server(source: str, server: str = None) -> str:
    """Server of `source` (its default when empty); refuse the other backend's API host."""
    server = (server or BACKENDS[source].DEFAULT_SERVER).rstrip("/")
    host = urlparse(server).hostname or ""
    for other, suffix in SERVER_HOSTS.items():
        if other != source and host.endswith(suffix):
            raise ValueError(
                f"{LOG} --server {server} is the {other} API, but --source is {source}. "
                f"Set the {source} server (default {BACKENDS[source].DEFAULT_SERVER}) or use --source {other}."
            )
    return server


def run_annotate_rna_genes(
    input_paths: list,
    output_path: str,
    source: str = DEFAULT_SOURCE,
    server: str = None,
    cache_path: str = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
    max_workers: int = None,
    attrs_url: str = None,
    client=None,
) -> pd.DataFrame:
    if source not in BACKENDS:
        raise ValueError(f"{LOG} Unknown --source {source!r}; choose one of {', '.join(SOURCES)}.")
    backend = BACKENDS[source]
    server = resolve_server(source, server)
    cache_path = str(cache_path or default_cache_path(DEFAULT_CACHE_DIR, source))
    snps = significant_snps(input_paths)
    print(f"{LOG} {len(input_paths)} traits: {len(snps):,} unique BH-significant SNPs to annotate "
          f"(source: {source}).")

    cache = load_backend_cache(cache_path, source)
    uncached = [s for s in snps if s not in cache]
    print(f"{LOG} Cache hits: {len(snps) - len(uncached):,} | Uncached SNPs to query: {len(uncached):,}")

    if uncached:
        options = {"attrs_url": attrs_url} if source == "ucsc" else {}
        failed = backend.annotate(uncached, cache, cache_path, server=server, batch_size=max(1, batch_size),
                                  max_workers=max_workers, client=client, **options)
        if failed:
            raise RuntimeError(
                f"{LOG} {backend.SOURCE_NAME} lookups failed for {len(failed):,} SNPs (e.g. {failed[:5]}). "
                f"Successful results were cached; re-run the pipeline to retry the remaining SNPs."
            )

    out = pd.DataFrame([annotation_row(s, cache[s]) for s in snps], columns=OUTPUT_COLUMNS)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    counts = out["STATUS"].value_counts().to_dict()
    print(f"{LOG} {len(out):,} SNPs {counts} -> {output_path}")
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Annotate BH-significant RNA m6A SNPs with the hg19/GRCh37 genes they fall inside."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Step B1 results/rna/{trait}_overlap_bh.tsv files")
    parser.add_argument("--output", "-o", required=True, help="Output results/rna/rna_sig_genes.tsv")
    parser.add_argument(
        "--source", choices=SOURCES, default=DEFAULT_SOURCE,
        help="Annotation backend: ucsc = UCSC Genome Browser hg19 / GENCODE V19 (default); "
             "ensembl = Ensembl VEP REST, GRCh37",
    )
    parser.add_argument(
        "--server", default=None,
        help=f"REST server of the chosen source (defaults: ucsc {rna_genes_ucsc.DEFAULT_SERVER}, "
             f"ensembl {rna_genes_ensembl.DEFAULT_SERVER})",
    )
    parser.add_argument(
        "--cache-file", default=None,
        help="Persistent JSON annotation cache of the chosen source "
             "(default: results/.cache/ucsc_hg19_genes_cache.json or ensembl_grch37_genes_cache.json)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help=f"SNPs per cache checkpoint / VEP request (default: {DEFAULT_BATCH_SIZE})",
    )
    parser.add_argument(
        "--workers", type=int, default=None,
        help=f"Concurrent Ensembl VEP requests (default: {rna_genes_ensembl.DEFAULT_WORKERS}); "
             "UCSC requests are always sequential (~1 per second)",
    )
    return parser


def main(argv: list = None) -> None:
    args = build_parser().parse_args(argv)
    run_annotate_rna_genes(
        args.inputs,
        args.output,
        source=args.source,
        server=args.server,
        cache_path=args.cache_file,
        batch_size=args.batch_size,
        max_workers=args.workers,
    )


if __name__ == "__main__":
    main()
