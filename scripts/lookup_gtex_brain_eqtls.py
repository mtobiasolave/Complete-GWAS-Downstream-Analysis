#!/usr/bin/env python3
"""
Step A3: Query GTEx v10 Across 13 Brain Tissues for BH-Significant DNA ASM Variants.

Key Features:
1. Queries every unique QUERY_SNP in `{trait}_overlap_sig.tsv` (BH FDR <= 5% at Step A1).
2. Uses a persistent atomic JSON cache (`results/.cache/gtex_eqtl_cache.json`)
   shared across traits so overlapping SNPs are never queried twice.
3. Resolves rsID -> GTEx variantId in parallel threads and queries the GTEx v2
   `singleTissueEqtl` endpoint in batches of 40 variantIds per HTTP request across
   all 13 brain regions simultaneously.
4. Collapses multi-tissue eQTL hits per (QUERY_SNP, GTEX_GENE_SYMBOL) into:
   QUERY_SNP, GTEX_GENE_SYMBOL, GTEX_TISSUES, N_BRAIN_TISSUES, MIN_GTEX_PVALUE, MEAN_GTEX_SLOPE.
"""

import argparse
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests

from common_utils import (
    DEFAULT_BRAIN_TISSUES,
    TISSUE_TYPO_FIXES,
    PROJECT_ROOT,
    create_http_session,
    load_json_cache,
    normalize_tissues,
    save_json_cache,
)

DEFAULT_CACHE_FILE = str(PROJECT_ROOT / "results" / ".cache" / "gtex_eqtl_cache.json")


def fetch_variant_id(rsid: str, session: requests.Session) -> tuple:
    """Resolve rsID to GTEx v10 variantId (e.g. chr12_2203945_A_G_b38).

    Returns (rsid, variant_id_or_None, resolved); resolved is False when the request
    itself failed, so the SNP is not cached as "absent from GTEx".
    """
    url = "https://gtexportal.org/api/v2/dataset/variant"
    params = {"snpId": rsid, "datasetId": "gtex_v10"}
    for attempt in range(3):
        try:
            resp = session.get(url, params=params, timeout=20)
            if resp.status_code == 200:
                data = resp.json().get("data", [])
                if data:
                    vid = data[0].get("variantId")
                    if vid:
                        return rsid, vid, True
                return rsid, None, True
            if resp.status_code == 429:
                time.sleep(1.5 * (attempt + 1))
                continue
        except Exception:
            time.sleep(0.5 * (attempt + 1))
    return rsid, None, False


def fetch_eqtl_batch(
    variant_ids: list, tissues: list, session: requests.Session
) -> list:
    """Query singleTissueEqtl for a batch of GTEx variantIds across configured brain tissues.

    Returns (items, complete); complete is False if any page could not be retrieved.
    """
    url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
    all_items = []
    page = 0
    while True:
        params = {
            "variantId": variant_ids,
            "datasetId": "gtex_v10",
            "tissueSiteDetailId": tissues,
            "itemsPerPage": 100000,
            "page": page,
        }
        success = False
        for attempt in range(3):
            try:
                resp = session.get(url, params=params, timeout=30)
                if resp.status_code == 200:
                    payload = resp.json()
                    data = payload.get("data", [])
                    all_items.extend(data)
                    paging = payload.get("paging_info", {})
                    max_pages = paging.get("numberOfPages", 1)
                    page += 1
                    success = True
                    if page >= max_pages or not data:
                        return all_items, True
                    break
                if resp.status_code == 429:
                    time.sleep(1.5 * (attempt + 1))
            except Exception:
                time.sleep(0.5 * (attempt + 1))
        if not success:
            return all_items, False


def run_lookup_gtex(
    overlap_path: str,
    output_path: str,
    brain_tissues: list = None,
    p_threshold: float = 0.05,
    cache_path: str = DEFAULT_CACHE_FILE,
    max_workers: int = 16,
) -> pd.DataFrame:
    tissues = normalize_tissues(brain_tissues or DEFAULT_BRAIN_TISSUES)
    tissue_set = set(tissues)

    overlap_df = pd.read_csv(overlap_path, sep="\t", low_memory=False)
    if "QUERY_SNP" not in overlap_df.columns:
        raise ValueError(f"Missing QUERY_SNP column in {overlap_path}")

    overlap_df["QUERY_SNP"] = overlap_df["QUERY_SNP"].astype(str).str.strip()
    candidate_snps = sorted(overlap_df["QUERY_SNP"].dropna().unique())
    print(
        f"[lookup_gtex] {os.path.basename(overlap_path)}: {len(candidate_snps):,} "
        f"unique BH-significant QUERY_SNPs to annotate."
    )

    cache = load_json_cache(cache_path)
    uncached_snps = [s for s in candidate_snps if s not in cache]
    print(
        f"[lookup_gtex] Cache hits: {len(candidate_snps) - len(uncached_snps):,} | "
        f"Uncached SNPs to query: {len(uncached_snps):,}"
    )

    if uncached_snps:
        session = create_http_session(pool_size=max_workers + 4)
        vid_to_rsids = {}
        failed_snps = set()

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(fetch_variant_id, rsid, session): rsid
                for rsid in uncached_snps
            }
            for future in as_completed(futures):
                rsid, vid, resolved = future.result()
                if not resolved:
                    failed_snps.add(rsid)
                elif vid:
                    vid_to_rsids.setdefault(vid, []).append(rsid)

        new_hits = {rsid: [] for rsid in uncached_snps if rsid not in failed_snps}

        all_vids = list(vid_to_rsids.keys())
        batch_size = 40
        batches = [
            all_vids[i : i + batch_size] for i in range(0, len(all_vids), batch_size)
        ]

        with ThreadPoolExecutor(max_workers=min(8, max(1, len(batches)))) as executor:
            batch_futures = {
                executor.submit(fetch_eqtl_batch, batch, tissues, session): batch
                for batch in batches
            }
            for future in as_completed(batch_futures):
                items, complete = future.result()
                if not complete:
                    for vid in batch_futures[future]:
                        failed_snps.update(vid_to_rsids[vid])
                    continue
                for item in items:
                    vid = item.get("variantId")
                    snp_id = item.get("snpId")
                    target_rsids = vid_to_rsids.get(vid, [snp_id] if snp_id else [])
                    hit_record = {
                        "geneSymbol": item.get("geneSymbol"),
                        "tissueSiteDetailId": item.get("tissueSiteDetailId"),
                        "pValue": item.get("pValue"),
                        "nes": item.get("nes"),
                    }
                    for rsid in target_rsids:
                        if rsid in new_hits:
                            new_hits[rsid].append(hit_record)

        # Only fully retrieved SNPs are cached; failed requests are retried on the next run
        for rsid in failed_snps:
            new_hits.pop(rsid, None)
        cache.update(new_hits)
        save_json_cache(cache_path, new_hits)
        print(f"[lookup_gtex] Cached {len(new_hits):,} newly queried SNPs at {cache_path}")

        if failed_snps:
            raise RuntimeError(
                f"[lookup_gtex] GTEx requests failed for {len(failed_snps):,} SNPs "
                f"(e.g. {sorted(failed_snps)[:5]}). Successful results were cached; "
                f"re-run the pipeline to retry the remaining SNPs."
            )

    raw_rows = []
    all_overlap_snps = overlap_df["QUERY_SNP"].dropna().unique()
    for rsid in all_overlap_snps:
        hits = cache.get(rsid, [])
        if not isinstance(hits, list):
            continue
        for h in hits:
            raw_tissue = str(h.get("tissueSiteDetailId", ""))
            tissue = TISSUE_TYPO_FIXES.get(raw_tissue, raw_tissue)
            gene = str(h.get("geneSymbol", "")).strip()
            pval = pd.to_numeric(h.get("pValue"), errors="coerce")
            slope = pd.to_numeric(h.get("nes"), errors="coerce")
            if (
                gene
                and tissue in tissue_set
                and pd.notna(pval)
                and float(pval) <= p_threshold
            ):
                raw_rows.append(
                    {
                        "QUERY_SNP": rsid,
                        "GTEX_GENE_SYMBOL": gene,
                        "TISSUE": tissue,
                        "PVALUE": float(pval),
                        "SLOPE": float(slope) if pd.notna(slope) else np.nan,
                    }
                )

    out_cols = [
        "QUERY_SNP",
        "GTEX_GENE_SYMBOL",
        "GTEX_TISSUES",
        "N_BRAIN_TISSUES",
        "MIN_GTEX_PVALUE",
        "MEAN_GTEX_SLOPE",
    ]

    if not raw_rows:
        agg_df = pd.DataFrame(columns=out_cols)
    else:
        hits_df = pd.DataFrame(raw_rows)
        hits_df = hits_df.sort_values("PVALUE", ascending=True).drop_duplicates(
            subset=["QUERY_SNP", "GTEX_GENE_SYMBOL", "TISSUE"]
        )

        grouped = (
            hits_df.groupby(["QUERY_SNP", "GTEX_GENE_SYMBOL"], as_index=False)
            .agg(
                GTEX_TISSUES=("TISSUE", lambda s: ", ".join(sorted(set(s)))),
                N_BRAIN_TISSUES=("TISSUE", "nunique"),
                MIN_GTEX_PVALUE=("PVALUE", "min"),
                MEAN_GTEX_SLOPE=("SLOPE", "mean"),
            )
            .sort_values(
                by=["MIN_GTEX_PVALUE", "N_BRAIN_TISSUES"], ascending=[True, False]
            )
        )
        agg_df = grouped[out_cols]

    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    agg_df.to_csv(output_path, sep="\t", index=False, na_rep="NA")
    print(
        f"[lookup_gtex] Wrote {len(agg_df):,} aggregated (QUERY_SNP, GTEX_GENE_SYMBOL) pairs to {output_path}"
    )
    return agg_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Query GTEx v10 across 13 brain tissues for BH-significant DNA ASM variants."
    )
    parser.add_argument("--overlap", required=True, help="Input results/dna/{trait}_overlap_sig.tsv")
    parser.add_argument("--output", "-o", required=True, help="Output results/dna/{trait}_gtex_brain_eqtls.tsv")
    parser.add_argument("--tissues", nargs="+", default=DEFAULT_BRAIN_TISSUES, help="GTEx brain tissue IDs")
    parser.add_argument("--p-threshold", type=float, default=0.05, help="GTEx eQTL p-value cutoff")
    parser.add_argument("--cache-file", default=DEFAULT_CACHE_FILE, help="Persistent JSON cache path")
    parser.add_argument("--workers", type=int, default=16, help="Max concurrent HTTP workers")
    args = parser.parse_args()

    run_lookup_gtex(
        args.overlap,
        args.output,
        brain_tissues=args.tissues,
        p_threshold=args.p_threshold,
        cache_path=args.cache_file,
        max_workers=args.workers,
    )


if __name__ == "__main__":
    main()
