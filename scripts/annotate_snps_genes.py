#!/usr/bin/env python3
"""
Annotate Cross-Trait Methylation SNPs with Ensembl VEP Host/Consequence Genes and GTEx v10 Brain eQTLs.

Uses shared HTTP connection pooling, persistent JSON caching, and vectorized/list-comprehension
gene concordance matching via `common_utils`.
"""

import argparse
import glob
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from common_utils import (
    DEFAULT_BRAIN_TISSUES,
    create_http_session,
    detect_column,
    find_common_genes as _find_common_genes_base,
    load_json_cache,
    save_json_cache,
)

COMPARISONS_DIR = "results/methylation/comparisons"
DEFAULT_CACHE_FILE = "results/methylation/.api_cache.json"
GTEX_BRAIN_TISSUES = DEFAULT_BRAIN_TISSUES
cache_lock = threading.Lock()


def clean_rsid(val):
    """Sanitizes variant strings into a standard rsID format."""
    val_str = str(val).strip().lower()
    if not val_str or val_str in ("nan", "none", ""):
        return None
    if "." in val_str:
        val_str = val_str.split(".")[0]
    if not val_str.startswith("rs"):
        val_str = "rs" + val_str
    return val_str


def get_ensembl_genes(rsid, session, retries=3):
    """Fetch structurally associated/consequence genes via Ensembl VEP REST API."""
    if not rsid:
        return "None"

    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?"
    for attempt in range(retries):
        try:
            response = session.get(
                url,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=25,
            )
            if response.status_code == 200:
                data = response.json()
                if not data or not isinstance(data, list):
                    return "None"

                genes = []
                for entry in data:
                    consequences = entry.get("transcript_consequences", [])
                    for cons in consequences:
                        gene_id = cons.get("gene_symbol") or cons.get("gene_id")
                        if gene_id:
                            genes.append(gene_id)

                unique_genes = sorted(list(set(filter(None, genes))))
                return ", ".join(unique_genes) if unique_genes else "None"
            if response.status_code in (429, 500, 502, 503, 504):
                retry_after = float(response.headers.get("Retry-After", 1.5 * (attempt + 1)))
                time.sleep(retry_after)
                continue
        except Exception as e:
            if attempt == retries - 1:
                print(f"    [!] Ensembl VEP Network Error for {rsid}: {e}")
            time.sleep(1.0 * (attempt + 1))
    return "None"


_ENSEMBL_VARIANT_IDS = {}


def batch_get_ensembl_genes(rsids, session, chunk_size=10, max_workers=12):
    """Batch-fetch Ensembl VEP gene consequences for a list of rsIDs via POST /vep/human/id."""
    results = {}
    clean_list = [r for r in dict.fromkeys(rsids) if r]
    if not clean_list:
        return results

    url = "https://rest.ensembl.org/vep/human/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    chunks = [clean_list[i : i + chunk_size] for i in range(0, len(clean_list), chunk_size)]


    def _fetch_chunk(chunk):
        chunk_res = {}
        chunk_vids = {}
        for attempt in range(3):
            try:
                resp = session.post(url, headers=headers, json={"ids": chunk}, timeout=35)
                if resp.status_code == 200:
                    for r in chunk:
                        chunk_res[r] = "None"
                    data = resp.json()
                    if isinstance(data, list):
                        for entry in data:
                            qid = clean_rsid(entry.get("id") or entry.get("input"))
                            if not qid:
                                continue
                            seq_reg = str(entry.get("seq_region_name", "")).strip()
                            pos = entry.get("start")
                            alleles = str(entry.get("allele_string", "")).strip().split("/")
                            if (
                                seq_reg
                                and pos
                                and len(alleles) >= 2
                                and not seq_reg.startswith("CHR_")
                                and not seq_reg.startswith("GL")
                                and not seq_reg.startswith("KI")
                            ):
                                chrom_str = seq_reg if seq_reg.startswith("chr") else f"chr{seq_reg}"
                                ref_a = alleles[0]
                                vids = [
                                    f"{chrom_str}_{pos}_{ref_a}_{alt_a}_b38"
                                    for alt_a in alleles[1:]
                                    if ref_a.isalpha() and alt_a.isalpha()
                                ]
                                if vids:
                                    chunk_vids.setdefault(qid, []).extend(vids)
                            genes = []
                            for cons in entry.get("transcript_consequences", []):
                                g = cons.get("gene_symbol") or cons.get("gene_id")
                                if g:
                                    genes.append(g)
                            unique_genes = sorted(list(set(filter(None, genes))))
                            chunk_res[qid] = ", ".join(unique_genes) if unique_genes else "None"
                    return chunk_res, chunk_vids
                if resp.status_code in (429, 500, 502, 503, 504):
                    time.sleep(float(resp.headers.get("Retry-After", 1.5 * (attempt + 1))))
            except Exception:
                time.sleep(1.0 * (attempt + 1))
        return chunk_res, chunk_vids

    completed_chunks = 0
    with ThreadPoolExecutor(max_workers=min(max_workers, max(1, len(chunks)))) as executor:
        futures = [executor.submit(_fetch_chunk, ch) for ch in chunks]
        for future in as_completed(futures):
            c_res, c_vids = future.result()
            results.update(c_res)
            with cache_lock:
                _ENSEMBL_VARIANT_IDS.update(c_vids)
            completed_chunks += 1
            if completed_chunks % 20 == 0 or completed_chunks == len(chunks):
                print(
                    f"    [Ensembl VEP] {completed_chunks}/{len(chunks)} chunks complete ({len(results):,}/{len(clean_list):,} SNPs)",
                    flush=True,
                )

    missing = [r for r in clean_list if r not in results]
    if missing:
        with ThreadPoolExecutor(max_workers=8) as executor:
            futs = {executor.submit(get_ensembl_genes, r, session): r for r in missing}
            for fut in as_completed(futs):
                results[futs[fut]] = fut.result()
    return results





def get_gtex_brain_eqtls(rsid, session):
    """Fetch significant eQTL target genes across brain tissues using GTEx Portal API v2."""
    if not rsid:
        return "None"

    variant_url = f"https://gtexportal.org/api/v2/dataset/variant?snpId={rsid}&datasetId=gtex_v10"
    try:
        var_resp = session.get(variant_url, timeout=25)
        if var_resp.status_code == 200:
            var_data = var_resp.json().get("data", [])
            if not var_data:
                return "None"

            variant_id = var_data[0].get("variantId")
            if not variant_id:
                return "None"

            eqtl_url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
            params = {
                "variantId": variant_id,
                "datasetId": "gtex_v10",
                "tissueSiteDetailId": GTEX_BRAIN_TISSUES,
                "itemsPerPage": 500,
            }

            eqtl_resp = session.get(eqtl_url, params=params, timeout=25)
            if eqtl_resp.status_code == 200:
                eqtl_data = eqtl_resp.json().get("data", [])
                genes = {item.get("geneSymbol") for item in eqtl_data if item.get("geneSymbol")}
                unique_genes = sorted(list(genes))
                return ", ".join(unique_genes) if unique_genes else "None"
    except Exception as e:
        print(f"    [!] GTEx API Network Error for {rsid}: {e}")
    return "None"


def find_common_genes(ensembl_str, gtex_str):
    """Finds overlapping genes between Ensembl and GTEx annotations."""
    res = _find_common_genes_base(ensembl_str, gtex_str)
    return "No match" if res == "None" else res


def generate_summary_file(df, output_summary_path, comment_lines=None):
    """Creates a concise TSV containing SNP, trait recurrence count, and common genes."""
    snp_col = detect_column(df.columns, "SNP") or "SNP"
    cols_upper = {c.upper(): c for c in df.columns}
    trait_count_candidates = ["NUM_SIG_TRAITS", "NUM_TRAITS", "SIGNIFICANT_TRAITS", "TRAIT_COUNT"]
    trait_col = next((cols_upper[c] for c in trait_count_candidates if c in cols_upper), None)

    summary_records = {"SNP": df[snp_col]}
    if trait_col:
        summary_records["Num_Sig_Traits"] = df[trait_col]

    ens_list = df["Ensembl_Genes"].tolist() if "Ensembl_Genes" in df.columns else [""] * len(df)
    gtex_list = (
        df["GTEx_Brain_eQTL_Genes"].tolist()
        if "GTEx_Brain_eQTL_Genes" in df.columns
        else [""] * len(df)
    )
    summary_records["Common_Genes"] = [
        find_common_genes(ens_val, gtex_val) for ens_val, gtex_val in zip(ens_list, gtex_list)
    ]

    summary_df = pd.DataFrame(summary_records)

    os.makedirs(
        os.path.dirname(output_summary_path) if os.path.dirname(output_summary_path) else ".",
        exist_ok=True,
    )
    with open(output_summary_path, "w", encoding="utf-8") as f:
        if comment_lines:
            for cl in comment_lines:
                f.write(cl)
    summary_df.to_csv(output_summary_path, sep="\t", index=False, mode="a")
    print(f"[OK] Summary file written to: {output_summary_path}\n")


def fetch_worker(row_idx, raw_snp, total_snps, session, cache=None, use_cache=True):
    """Worker task that queries APIs for a single SNP, checking cache first."""
    clean_snp = clean_rsid(raw_snp)
    is_cached = False
    ensembl_genes = "None"
    gtex_genes = None

    if use_cache and cache is not None and clean_snp:
        with cache_lock:
            cached_data = cache.get(clean_snp)
        if cached_data and isinstance(cached_data, dict):
            ensembl_genes = cached_data.get("ensembl", "None")
            if "gtex" in cached_data and cached_data.get("gtex") is not None:
                gtex_genes = cached_data.get("gtex")
                is_cached = True

    if gtex_genes is None:
        if ensembl_genes == "None":
            ensembl_genes = get_ensembl_genes(clean_snp, session)
        gtex_genes = get_gtex_brain_eqtls(clean_snp, session)
        if use_cache and cache is not None and clean_snp:
            with cache_lock:
                cache[clean_snp] = {"ensembl": ensembl_genes, "gtex": gtex_genes}

    return row_idx, raw_snp, ensembl_genes, gtex_genes, is_cached


def batch_get_gtex_brain_eqtls(rsids, session, max_workers=16):
    """Batch-fetch GTEx v10 brain eQTL target genes across 13 brain tissues (40 variantIds/request)."""
    from lookup_gtex_brain_eqtls import (
        DEFAULT_CACHE_FILE as GTEX_EQTL_CACHE_FILE,
        fetch_eqtl_batch,
        fetch_variant_id,
    )

    results = {}
    clean_list = [r for r in dict.fromkeys(rsids) if r]
    if not clean_list:
        return results

    # 1. Reuse existing Step 5c GTEx cache if available
    gtex_cache = load_json_cache(GTEX_EQTL_CACHE_FILE)
    tissue_set = set(GTEX_BRAIN_TISSUES)
    uncached = []
    for r in clean_list:
        if r in gtex_cache and isinstance(gtex_cache[r], list):
            genes = sorted({
                str(h.get("geneSymbol", "")).strip()
                for h in gtex_cache[r]
                if str(h.get("geneSymbol", "")).strip()
                and str(h.get("tissueSiteDetailId", "")) in tissue_set
                and pd.notna(pd.to_numeric(h.get("pValue"), errors="coerce"))
                and float(h.get("pValue")) <= 0.05
            })
            results[r] = ", ".join(genes) if genes else "None"
        else:
            uncached.append(r)

    if not uncached:
        return results

    # 2. Resolve rsID -> GTEx variantId (using Ensembl VEP GRCh38 coordinates first, then GTEx API fallback)
    vid_to_rsids = {}
    need_vid_lookup = []
    for r in uncached:
        ens_vids = _ENSEMBL_VARIANT_IDS.get(r)
        if ens_vids:
            for v in ens_vids:
                vid_to_rsids.setdefault(v, []).append(r)
        else:
            need_vid_lookup.append(r)

    if need_vid_lookup:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(fetch_variant_id, r, session): r for r in need_vid_lookup}
            for future in as_completed(futures):
                r, vid = future.result()
                if vid:
                    vid_to_rsids.setdefault(vid, []).append(r)
                else:
                    results[r] = "None"


    # 3. Query singleTissueEqtl in batches of 40 variantIds
    rsid_genes = {r: set() for r in uncached}
    all_vids = list(vid_to_rsids.keys())
    batches = [all_vids[i : i + 40] for i in range(0, len(all_vids), 40)]

    with ThreadPoolExecutor(max_workers=min(8, max(1, len(batches)))) as executor:
        batch_futures = [
            executor.submit(fetch_eqtl_batch, batch, GTEX_BRAIN_TISSUES, session)
            for batch in batches
        ]
        for future in as_completed(batch_futures):
            items = future.result()
            for item in items:
                vid = item.get("variantId")
                snp_id = clean_rsid(item.get("snpId"))
                gene = str(item.get("geneSymbol", "")).strip()
                pval = pd.to_numeric(item.get("pValue"), errors="coerce")
                if not gene or pd.isna(pval) or float(pval) > 0.05:
                    continue
                target_rsids = vid_to_rsids.get(vid, [snp_id] if snp_id else [])
                for r in target_rsids:
                    if r in rsid_genes:
                        rsid_genes[r].add(gene)

    for r in uncached:
        genes_sorted = sorted(rsid_genes.get(r, set()))
        results[r] = ", ".join(genes_sorted) if genes_sorted else "None"

    return results


def annotate_file(
    input_path,
    output_path,
    create_summary=None,
    max_workers=5,
    cache_path=DEFAULT_CACHE_FILE,
    use_cache=True,
):
    """Reads input TSV, queries APIs concurrently across worker threads with caching, and saves outputs."""
    if not os.path.exists(input_path):
        print(f"Error: File '{input_path}' not found.")
        return None, None

    comment_lines = []
    with open(input_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                comment_lines.append(line)
            else:
                break

    df = pd.read_csv(input_path, sep="\t", comment="#")
    snp_col = detect_column(df.columns, "SNP")
    if not snp_col:
        print(f"Error: Could not locate a variant ID (SNP/RSID) column in {input_path}")
        return None, None

    raw_snps = df[snp_col].tolist()
    total_snps = len(raw_snps)
    print(f"\nAnnotating {total_snps:,} SNPs from: {os.path.basename(input_path)}")
    print(f"Querying with {max_workers} parallel workers: Ensembl VEP & GTEx v10 (13 Brain Tissues)")

    cache = load_json_cache(cache_path) if use_cache else {}
    if use_cache and cache:
        print(f"Loaded {len(cache):,} cached variant annotations from {cache_path}")

    session = create_http_session(pool_size=max(20, max_workers * 2))

    all_clean_snps = [clean_rsid(s) for s in raw_snps]
    unique_clean_snps = [s for s in dict.fromkeys(all_clean_snps) if s]

    need_ensembl = [
        s
        for s in unique_clean_snps
        if s not in cache
        or not isinstance(cache.get(s), dict)
        or cache[s].get("ensembl") in (None, "")
    ]
    need_gtex = [
        s
        for s in unique_clean_snps
        if s not in cache
        or not isinstance(cache.get(s), dict)
        or cache[s].get("gtex") is None
    ]

    cached_hits = total_snps - max(len(need_ensembl), len(need_gtex))

    try:
        if need_ensembl:
            print(f"  -> Batch-querying Ensembl VEP for {len(need_ensembl):,} uncached SNPs...")
            batch_ens = batch_get_ensembl_genes(need_ensembl, session)
            with cache_lock:
                for s, ens_val in batch_ens.items():
                    if s in cache and isinstance(cache[s], dict):
                        cache[s]["ensembl"] = ens_val
                    else:
                        cache[s] = {"ensembl": ens_val, "gtex": None}

        if need_gtex:
            print(f"  -> Batch-querying GTEx v10 (13 Brain Tissues) for {len(need_gtex):,} uncached SNPs...")
            batch_gtex = batch_get_gtex_brain_eqtls(
                need_gtex, session, max_workers=max(16, max_workers * 2)
            )
            with cache_lock:
                for s, gtex_val in batch_gtex.items():
                    if s in cache and isinstance(cache[s], dict):
                        cache[s]["gtex"] = gtex_val
                    else:
                        cache[s] = {"ensembl": "None", "gtex": gtex_val}
    finally:
        if use_cache and cache:
            save_json_cache(cache_path, cache)

    ensembl_results = []
    gtex_results = []
    for clean_snp in all_clean_snps:
        entry = cache.get(clean_snp, {}) if isinstance(cache, dict) else {}
        ens_val = entry.get("ensembl", "None") if isinstance(entry, dict) else "None"
        gtex_val = entry.get("gtex", "None") if isinstance(entry, dict) else "None"
        ensembl_results.append(ens_val or "None")
        gtex_results.append(gtex_val or "None")

    print(
        f"API Queries complete: {max(0, cached_hits)}/{total_snps} loaded from cache, "
        f"{total_snps - max(0, cached_hits)} newly fetched."
    )

    df["Ensembl_Genes"] = ensembl_results
    df["GTEx_Brain_eQTL_Genes"] = gtex_results

    os.makedirs(
        os.path.dirname(output_path) if os.path.dirname(output_path) else ".",
        exist_ok=True,
    )
    with open(output_path, "w", encoding="utf-8") as f:
        for cl in comment_lines:
            f.write(cl)
    df.to_csv(output_path, sep="\t", index=False, mode="a")

    print(f"\nAnnotation complete! Full results written to: {output_path}")

    summary_path = None
    if create_summary is None:
        choice = (
            input(
                "\nGenerate a summarized file (SNP, trait count, matching Ensembl/GTEx genes)? [y/N]: "
            )
            .strip()
            .lower()
        )
        should_summarize = choice in ["y", "yes"]
    else:
        should_summarize = create_summary

    if should_summarize:
        default_summary = output_path.replace("_annotated.tsv", "_summary.tsv")
        if default_summary == output_path:
            default_summary = output_path.replace(".tsv", "_summary.tsv")

        if create_summary is None:
            sum_input = input(f"Summary output path [default: {default_summary}]: ").strip()
            summary_path = sum_input if sum_input else default_summary
        else:
            summary_path = default_summary

        generate_summary_file(df, summary_path, comment_lines)
    else:
        print("Skipping summarized file generation.\n")

    return output_path, summary_path


def interactive_mode():
    """Interactive selector for comparison files."""
    files = sorted(glob.glob(os.path.join(COMPARISONS_DIR, "*.tsv")))
    files = [
        f
        for f in files
        if not f.endswith("_annotated.tsv") and not f.endswith("_summary.tsv")
    ]

    if not files:
        print(
            f"No comparison files found in '{COMPARISONS_DIR}/'. Run compare_overlap_snps.py first."
        )
        sys.exit(1)

    print("==================================================")
    print("      SNP Gene Annotation Tool (Ensembl & GTEx)   ")
    print("==================================================")
    print("\nAvailable comparison files:")
    for idx, f in enumerate(files, 1):
        print(f"  [{idx:2d}] {os.path.basename(f)}")

    choice = input(f"\nSelect file to annotate (1-{len(files)}): ").strip()
    try:
        selected_file = files[int(choice) - 1]
    except (ValueError, IndexError):
        print("Invalid selection.")
        sys.exit(1)

    default_out = selected_file.replace(".tsv", "_annotated.tsv")
    custom_out = input(f"Output file path [default: {default_out}]: ").strip()
    out_file = custom_out if custom_out else default_out

    annotate_file(selected_file, out_file)


def main():
    parser = argparse.ArgumentParser(
        description="Annotate SNPs with Ensembl VEP proximal genes and GTEx brain eQTLs."
    )
    parser.add_argument("--input", "-i", help="Path to input TSV file containing SNP column")
    parser.add_argument("--output", "-o", help="Path for output annotated TSV file")
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=5,
        help="Number of concurrent API worker threads (default: 5)",
    )
    parser.add_argument(
        "--summary",
        dest="create_summary",
        action="store_true",
        default=None,
        help="Generate the simplified summary TSV without prompting",
    )
    parser.add_argument(
        "--no-summary",
        dest="create_summary",
        action="store_false",
        help="Skip generating the simplified summary TSV without prompting",
    )
    parser.add_argument(
        "--cache-file",
        default=DEFAULT_CACHE_FILE,
        help=f"Path to local API JSON cache (default: {DEFAULT_CACHE_FILE})",
    )
    parser.add_argument(
        "--no-cache", action="store_true", help="Bypass local cache and force API queries"
    )

    args = parser.parse_args()

    if args.input:
        output_file = args.output or args.input.replace(".tsv", "_annotated.tsv")
        annotate_file(
            args.input,
            output_file,
            create_summary=args.create_summary,
            max_workers=args.threads,
            cache_path=args.cache_file,
            use_cache=not args.no_cache,
        )
    else:
        interactive_mode()


if __name__ == "__main__":
    main()