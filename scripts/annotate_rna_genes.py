#!/usr/bin/env python3
"""
Step B2: Positional Gene Annotation of BH-Significant RNA m6A SNPs via Ensembl VEP.

1. Collects every SNP that is BH-significant in at least one trait (Step B1 outputs).
2. Queries the Ensembl VEP REST API (`POST /vep/homo_sapiens/id`, GRCh37 server by
   default to match the hg19 coordinates of the DNA track) in batches of rsIDs.
3. Keeps the genes whose transcripts the variant physically falls inside: any transcript
   consequence other than upstream/downstream_gene_variant (VEP reports those up to 5 kb
   away). Only loci on primary chromosomes (1-22, X, Y, MT) are used.
4. Uses a persistent atomic JSON cache (`results/.cache/ensembl_grch37_genes_cache.json`).
   Only successfully queried SNPs are cached; failed requests make the job fail so the
   remaining SNPs are retried on the next run.

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
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from common_utils import PROJECT_ROOT, create_http_session, load_json_cache, save_json_cache
from pleiotropy import is_true

DEFAULT_SERVER = "https://grch37.rest.ensembl.org"
DEFAULT_CACHE_FILE = str(PROJECT_ROOT / "results" / ".cache" / "ensembl_grch37_genes_cache.json")

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


def fetch_release(server: str, session) -> str:
    """Ensembl release served by `server` (recorded for reproducibility)."""
    try:
        resp = session.get(f"{server}/info/data", headers={"Accept": "application/json"}, timeout=30)
        if resp.status_code == 200:
            return ", ".join(str(r) for r in resp.json().get("releases", [])) or "NA"
    except Exception:
        pass
    return "NA"


def fetch_vep_batch(rsids: list, server: str, session) -> tuple:
    """POST one batch to VEP. Returns ({input rsID: [VEP entries]}, complete)."""
    url = f"{server}/vep/homo_sapiens/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    for attempt in range(3):
        try:
            resp = session.post(url, headers=headers, json={"ids": rsids}, timeout=120)
            if resp.status_code == 200:
                entries = {}
                for entry in resp.json():
                    entries.setdefault(entry.get("input") or entry.get("id"), []).append(entry)
                return entries, True
            if resp.status_code == 429:
                time.sleep(float(resp.headers.get("Retry-After", 2 * (attempt + 1))))
                continue
            print(f"[annotate_rna_genes] VEP HTTP {resp.status_code}: {resp.text[:200]}")
        except Exception as exc:
            print(f"[annotate_rna_genes] VEP request error: {exc}")
        time.sleep(2 * (attempt + 1))
    return {}, False


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
    batch_size: int = 100,
    max_workers: int = 4,
) -> pd.DataFrame:
    server = server.rstrip("/")
    snps = significant_snps(input_paths)
    print(f"[annotate_rna_genes] {len(input_paths)} traits: {len(snps):,} unique BH-significant SNPs to annotate.")

    cache = load_json_cache(cache_path)
    uncached = [s for s in snps if s not in cache]
    print(f"[annotate_rna_genes] Cache hits: {len(snps) - len(uncached):,} | Uncached SNPs to query: {len(uncached):,}")

    if uncached:
        session = create_http_session(pool_size=max_workers + 2)
        release = fetch_release(server, session)
        print(f"[annotate_rna_genes] Querying Ensembl VEP at {server} (release {release})")
        batches = [uncached[i : i + batch_size] for i in range(0, len(uncached), batch_size)]
        new_records, failed_snps = {}, []

        with ThreadPoolExecutor(max_workers=min(max_workers, len(batches))) as executor:
            futures = {executor.submit(fetch_vep_batch, b, server, session): b for b in batches}
            for done, future in enumerate(as_completed(futures), start=1):
                entries, complete = future.result()
                batch = futures[future]
                if not complete:
                    failed_snps.extend(batch)
                    continue
                for rsid in batch:
                    new_records[rsid] = summarize_vep_entries(entries.get(rsid, []), release)
                print(f"[annotate_rna_genes] {done}/{len(batches)} batches complete", flush=True)

        # Only fully retrieved SNPs are cached; failed batches are retried on the next run
        cache.update(new_records)
        save_json_cache(cache_path, new_records)
        print(f"[annotate_rna_genes] Cached {len(new_records):,} newly queried SNPs at {cache_path}")

        if failed_snps:
            raise RuntimeError(
                f"[annotate_rna_genes] Ensembl VEP requests failed for {len(failed_snps):,} SNPs "
                f"(e.g. {sorted(failed_snps)[:5]}). Successful results were cached; "
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
    parser.add_argument("--batch-size", type=int, default=100, help="rsIDs per VEP POST request (default: 100)")
    parser.add_argument("--workers", type=int, default=4, help="Max concurrent VEP requests (default: 4)")
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
