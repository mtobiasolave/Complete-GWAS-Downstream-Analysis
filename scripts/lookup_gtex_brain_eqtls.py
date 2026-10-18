#!/usr/bin/env python3
"""
Step 5c: Query GTEx v10 Across 13 Brain Tissues for DNA ASM Overlapping Variants.

Key Improvements:
1. Prioritizes candidate variants reaching either candidate-space Storey q <= 0.05
   or GWAS_P <= gwas_p_threshold (default 1e-3), avoiding tens of thousands of
   redundant HTTP calls on null variants while supporting full queries if desired.
2. Uses a persistent atomic JSON cache (`results/methylation/.gtex_eqtl_cache.json`)
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
    create_http_session,
    load_json_cache,
    normalize_tissues,
    save_json_cache,
    storey_qvalues,
)

DEFAULT_CACHE_FILE = "results/methylation/.gtex_eqtl_cache.json"


def fetch_variant_id(rsid: str, session: requests.Session) -> tuple:
    """Resolve rsID to GTEx v10 variantId (e.g. chr12_2203945_A_G_b38)."""
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
                        return rsid, vid
                return rsid, None
            if resp.status_code == 429:
                time.sleep(1.5 * (attempt + 1))
                continue
        except Exception:
            time.sleep(0.5 * (attempt + 1))
    return rsid, None


def fetch_eqtl_batch(
    variant_ids: list, tissues: list, session: requests.Session
) -> list:
    """Query singleTissueEqtl for a batch of GTEx variantIds across configured brain tissues."""
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
                        return all_items
                    break
                if resp.status_code == 429:
                    time.sleep(1.5 * (attempt + 1))
            except Exception:
                time.sleep(0.5 * (attempt + 1))
        if not success:
            break
    return all_items


def select_candidate_snps(
    overlap_df: pd.DataFrame, gwas_p_threshold: float, fdr_threshold: float = 0.05
) -> list:
    """
    Select unique QUERY_SNPs to query in GTEx:
    Retains variants with GWAS_P <= gwas_p_threshold OR candidate-space Storey q <= fdr_threshold.
    """
    p_col = "GWAS_P" if "GWAS_P" in overlap_df.columns else None
    if p_col is None:
        p_candidates = [c for c in overlap_df.columns if c.upper().startswith("P")]
        p_col = p_candidates[0] if p_candidates else None

    unique_df = overlap_df.drop_duplicates(subset=["QUERY_SNP"]).copy()
    if not p_col or gwas_p_threshold >= 1.0:
        return unique_df["QUERY_SNP"].dropna().astype(str).str.strip().tolist()

    unique_df[p_col] = pd.to_numeric(unique_df[p_col], errors="coerce")
    unique_df = unique_df.dropna(subset=["QUERY_SNP", p_col])

    mask_p = unique_df[p_col] <= gwas_p_threshold
    if len(unique_df) > 10:
        qvals = storey_qvalues(unique_df[p_col].to_numpy(dtype=float))
        mask_q = qvals <= fdr_threshold
        selected = unique_df[mask_p | mask_q]
    else:
        selected = unique_df[mask_p]

    return selected["QUERY_SNP"].astype(str).str.strip().tolist()


def run_lookup_gtex(
    overlap_path: str,
    output_path: str,
    brain_tissues: list = None,
    p_threshold: float = 0.05,
    gwas_p_threshold: float = 0.001,
    cache_path: str = DEFAULT_CACHE_FILE,
    max_workers: int = 16,
) -> pd.DataFrame:
    tissues = normalize_tissues(brain_tissues or DEFAULT_BRAIN_TISSUES)
    tissue_set = set(tissues)

    overlap_df = pd.read_csv(overlap_path, sep="\t", low_memory=False)
    if "QUERY_SNP" not in overlap_df.columns:
        raise ValueError(f"Missing QUERY_SNP column in {overlap_path}")

    overlap_df["QUERY_SNP"] = overlap_df["QUERY_SNP"].astype(str).str.strip()
    total_unique = overlap_df["QUERY_SNP"].nunique()

    candidate_snps = select_candidate_snps(overlap_df, gwas_p_threshold=gwas_p_threshold)
    print(
        f"[lookup_gtex] {os.path.basename(overlap_path)}: {len(candidate_snps):,} / {total_unique:,} "
        f"unique QUERY_SNPs selected (GWAS P <= {gwas_p_threshold} or Storey q <= 0.05)."
    )

    cache = load_json_cache(cache_path)
    uncached_snps = [s for s in candidate_snps if s not in cache]
    print(
        f"[lookup_gtex] Cache hits: {len(candidate_snps) - len(uncached_snps):,} | "
        f"Uncached SNPs to query: {len(uncached_snps):,}"
    )

    if uncached_snps:
        session = create_http_session(pool_size=max_workers + 4)
        rsid_to_vid = {}
        vid_to_rsids = {}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(fetch_variant_id, rsid, session): rsid
                for rsid in uncached_snps
            }
            for future in as_completed(futures):
                rsid, vid = future.result()
                rsid_to_vid[rsid] = vid
                if vid:
                    vid_to_rsids.setdefault(vid, []).append(rsid)

        for rsid in uncached_snps:
            cache[rsid] = []

        all_vids = list(vid_to_rsids.keys())
        batch_size = 40
        batches = [
            all_vids[i : i + batch_size] for i in range(0, len(all_vids), batch_size)
        ]

        with ThreadPoolExecutor(max_workers=min(8, max(1, len(batches)))) as executor:
            batch_futures = [
                executor.submit(fetch_eqtl_batch, batch, tissues, session)
                for batch in batches
            ]
            for future in as_completed(batch_futures):
                items = future.result()
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
                        if rsid in cache:
                            cache[rsid].append(hit_record)

        save_json_cache(cache_path, cache)
        print(f"[lookup_gtex] Updated persistent GTEx cache at {cache_path}")

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
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        overlap_path = smk.input.overlap
        output_path = smk.output[0]
        tissues = getattr(smk.params, "brain_tissues", DEFAULT_BRAIN_TISSUES)
        p_thresh = float(getattr(smk.params, "p_threshold", 0.05))
        gwas_p_thresh = float(getattr(smk.params, "gwas_p_threshold", 0.001))
        cache_file = getattr(smk.params, "cache_file", DEFAULT_CACHE_FILE)
        run_lookup_gtex(
            overlap_path,
            output_path,
            brain_tissues=tissues,
            p_threshold=p_thresh,
            gwas_p_threshold=gwas_p_thresh,
            cache_path=cache_file,
        )
        return

    parser = argparse.ArgumentParser(
        description="Query GTEx v10 across 13 brain tissues for DNA ASM overlapping variants."
    )
    parser.add_argument("--overlap", required=True, help="Input dna_asm_{trait}_gwas_overlap.tsv")
    parser.add_argument("--output", "-o", required=True, help="Output dna_asm_{trait}_gtex_brain_eqtls.tsv")
    parser.add_argument("--tissues", nargs="+", default=DEFAULT_BRAIN_TISSUES, help="GTEx brain tissue IDs")
    parser.add_argument("--p-threshold", type=float, default=0.05, help="GTEx eQTL p-value cutoff")
    parser.add_argument(
        "--gwas-p-threshold",
        type=float,
        default=0.001,
        help="GWAS p-value pre-filter (variants with Storey q <= 0.05 are also included)",
    )
    parser.add_argument("--cache-file", default=DEFAULT_CACHE_FILE, help="Persistent JSON cache path")
    parser.add_argument("--workers", type=int, default=16, help="Max concurrent HTTP workers")
    args = parser.parse_args()

    run_lookup_gtex(
        args.overlap,
        args.output,
        brain_tissues=args.tissues,
        p_threshold=args.p_threshold,
        gwas_p_threshold=args.gwas_p_threshold,
        cache_path=args.cache_file,
        max_workers=args.workers,
    )


if __name__ == "__main__":
    main()
