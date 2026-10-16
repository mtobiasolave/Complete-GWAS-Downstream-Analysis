#!/usr/bin/env python3
import os
import sys
import glob
import time
import argparse
import requests
from requests.adapters import HTTPAdapter
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import threading
import pandas as pd

COMPARISONS_DIR = "results/methylation/comparisons"
DEFAULT_CACHE_FILE = "results/methylation/.api_cache.json"
cache_lock = threading.Lock()

def load_api_cache(cache_path=DEFAULT_CACHE_FILE):
    """Loads previously queried variant annotations from a local JSON cache."""
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                if isinstance(data, dict):
                    return data
        except Exception as e:
            print(f"    [!] Warning: Could not read cache from {cache_path} ({e}). Starting with empty cache.")
    return {}

def save_api_cache(cache, cache_path=DEFAULT_CACHE_FILE):
    """Atomically persists the API cache dictionary to disk."""
    if not cache:
        return
    try:
        os.makedirs(os.path.dirname(cache_path) if os.path.dirname(cache_path) else ".", exist_ok=True)
        tmp_path = f"{cache_path}.tmp"
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(cache, f, indent=2)
        if os.path.exists(cache_path):
            os.replace(tmp_path, cache_path)
        else:
            os.rename(tmp_path, cache_path)
    except Exception as e:
        print(f"    [!] Warning: Could not save API cache to {cache_path} ({e})")

# 13 GTEx Brain Tissues
GTEX_BRAIN_TISSUES = [
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra"
]

def create_http_session(pool_size=10):
    """Creates a shared HTTP session with persistent connection pooling and retries."""
    session = requests.Session()
    # Passing an integer lets requests manage retries natively without namespace conflicts
    adapter = HTTPAdapter(pool_connections=pool_size, pool_maxsize=pool_size, max_retries=3)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session

def clean_rsid(val):
    """Sanitizes variant strings into a standard rsID format."""
    val_str = str(val).strip().lower()
    if not val_str or val_str in ['nan', 'none', '']:
        return None
    if '.' in val_str:
        val_str = val_str.split('.')[0]
    if not val_str.startswith('rs'):
        val_str = 'rs' + val_str
    return val_str

def get_ensembl_genes(rsid, session):
    """Fetch structurally associated/consequence genes via Ensembl VEP REST API."""
    if not rsid:
        return "None"
        
    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?"
    try:
        response = session.get(url, headers={"Content-Type": "application/json"}, timeout=25)
        if response.status_code == 200:
            data = response.json()
            if not data or not isinstance(data, list):
                return "None"
                
            genes = []
            for entry in data:
                consequences = entry.get('transcript_consequences', [])
                for cons in consequences:
                    gene_id = cons.get('gene_symbol') or cons.get('gene_id')
                    if gene_id:
                        genes.append(gene_id)
            
            unique_genes = sorted(list(set(filter(None, genes))))
            return ", ".join(unique_genes) if unique_genes else "None"
        elif response.status_code == 429:
            retry_after = int(response.headers.get("Retry-After", 1))
            time.sleep(retry_after)
            return get_ensembl_genes(rsid, session)
    except Exception as e:
        print(f"    [!] Ensembl VEP Network Error for {rsid}: {e}")
    return "None"

def get_gtex_brain_eqtls(rsid, session):
    """Fetch significant eQTL target genes across brain tissues using GTEx Portal API v2."""
    if not rsid:
        return "None"
        
    variant_url = f"https://gtexportal.org/api/v2/dataset/variant?snpId={rsid}&datasetId=gtex_v10"
    try:
        var_resp = session.get(variant_url, timeout=25)
        if var_resp.status_code == 200:
            var_data = var_resp.json().get('data', [])
            if not var_data:
                return "None"
                
            variant_id = var_data[0].get('variantId')
            if not variant_id:
                return "None"
                
            eqtl_url = "https://gtexportal.org/api/v2/association/singleTissueEqtl"
            params = {
                "variantId": variant_id,
                "datasetId": "gtex_v10",
                "tissueSiteDetailId": GTEX_BRAIN_TISSUES,
                "itemsPerPage": 500
            }
            
            eqtl_resp = session.get(eqtl_url, params=params, timeout=25)
            if eqtl_resp.status_code == 200:
                eqtl_data = eqtl_resp.json().get('data', [])
                genes = {item.get('geneSymbol') for item in eqtl_data if item.get('geneSymbol')}
                unique_genes = sorted(list(genes))
                return ", ".join(unique_genes) if unique_genes else "None"
    except Exception as e:
        print(f"    [!] GTEx API Network Error for {rsid}: {e}")
    return "None"

def find_common_genes(ensembl_str, gtex_str):
    """Finds overlapping genes between Ensembl and GTEx annotations."""
    if not ensembl_str or str(ensembl_str).strip() in ["None", "NA", "nan", ""]:
        ens_map = {}
    else:
        ens_map = {
            g.strip().upper(): g.strip() 
            for g in str(ensembl_str).split(",") 
            if g.strip() and g.strip().upper() != "NONE"
        }

    if not gtex_str or str(gtex_str).strip() in ["None", "NA", "nan", ""]:
        gtex_map = {}
    else:
        gtex_map = {
            g.strip().upper(): g.strip() 
            for g in str(gtex_str).split(",") 
            if g.strip() and g.strip().upper() != "NONE"
        }

    common_keys = sorted(list(set(ens_map.keys()).intersection(set(gtex_map.keys()))))
    if not common_keys:
        return "No match"
    
    return ", ".join([ens_map[k] for k in common_keys])

def generate_summary_file(df, output_summary_path, comment_lines=None):
    """Creates a concise TSV containing SNP, trait recurrence count, and common genes."""
    cols = {c.upper(): c for c in df.columns}
    snp_col = next((cols[c] for c in ["SNP", "RSID", "MARKERNAME", "ID"] if c in cols), "SNP")
    
    trait_count_candidates = ["NUM_SIG_TRAITS", "NUM_TRAITS", "SIGNIFICANT_TRAITS", "TRAIT_COUNT"]
    trait_col = next((cols[c] for c in trait_count_candidates if c in cols), None)

    summary_records = {
        "SNP": df[snp_col]
    }
    
    if trait_col:
        summary_records["Num_Sig_Traits"] = df[trait_col]

    summary_records["Common_Genes"] = df.apply(
        lambda row: find_common_genes(row.get('Ensembl_Genes', ''), row.get('GTEx_Brain_eQTL_Genes', '')),
        axis=1
    )

    summary_df = pd.DataFrame(summary_records)

    os.makedirs(os.path.dirname(output_summary_path) if os.path.dirname(output_summary_path) else ".", exist_ok=True)
    with open(output_summary_path, "w") as f:
        if comment_lines:
            for cl in comment_lines:
                f.write(cl)
    summary_df.to_csv(output_summary_path, sep="\t", index=False, mode="a")
    print(f"[OK] Summary file written to: {output_summary_path}\n")

def fetch_worker(row_idx, raw_snp, total_snps, session, cache=None, use_cache=True):
    """Worker task that queries both APIs for a single SNP, checking cache first."""
    clean_snp = clean_rsid(raw_snp)
    is_cached = False
    
    if use_cache and cache is not None and clean_snp:
        with cache_lock:
            cached_data = cache.get(clean_snp)
        if cached_data and isinstance(cached_data, dict):
            ensembl_genes = cached_data.get("ensembl", "None")
            gtex_genes = cached_data.get("gtex", "None")
            is_cached = True

    if not is_cached:
        ensembl_genes = get_ensembl_genes(clean_snp, session)
        gtex_genes = get_gtex_brain_eqtls(clean_snp, session)
        if use_cache and cache is not None and clean_snp:
            with cache_lock:
                cache[clean_snp] = {
                    "ensembl": ensembl_genes,
                    "gtex": gtex_genes
                }

    return row_idx, raw_snp, ensembl_genes, gtex_genes, is_cached

def annotate_file(input_path, output_path, create_summary=None, max_workers=5, cache_path=DEFAULT_CACHE_FILE, use_cache=True):
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

    df = pd.read_csv(input_path, sep='\t', comment='#')
    
    cols = {c.upper(): c for c in df.columns}
    snp_col = next((cols[c] for c in ["SNP", "RSID", "MARKERNAME", "ID"] if c in cols), None)
    if not snp_col:
        print(f"Error: Could not locate a variant ID (SNP/RSID) column in {input_path}")
        return None, None

    total_snps = len(df)
    print(f"\nAnnotating {total_snps:,} SNPs from: {os.path.basename(input_path)}")
    print(f"Querying with {max_workers} parallel workers: Ensembl VEP & GTEx v10 (13 Brain Tissues)")

    cache = load_api_cache(cache_path) if use_cache else {}
    if use_cache and cache:
        print(f"Loaded {len(cache):,} cached variant annotations from {cache_path}")

    session = create_http_session(pool_size=max_workers * 2)
    
    ensembl_results = [None] * total_snps
    gtex_results = [None] * total_snps
    
    completed_count = 0
    cached_hits = 0
    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(fetch_worker, idx, df.iloc[idx][snp_col], total_snps, session, cache, use_cache): idx
                for idx in range(total_snps)
            }
            
            for future in as_completed(futures):
                idx, raw_snp, ensembl_genes, gtex_genes, is_cached = future.result()
                ensembl_results[idx] = ensembl_genes
                gtex_results[idx] = gtex_genes
                completed_count += 1
                if is_cached:
                    cached_hits += 1
                cache_tag = " [cached]" if is_cached else ""
                print(f"[{completed_count:3d}/{total_snps}] {raw_snp} -> Ensembl: {ensembl_genes} | GTEx: {gtex_genes}{cache_tag}")
    finally:
        if use_cache and cache:
            save_api_cache(cache, cache_path)

    print(f"API Queries complete: {cached_hits}/{total_snps} loaded from cache, {total_snps - cached_hits} newly fetched.")

    df['Ensembl_Genes'] = ensembl_results
    df['GTEx_Brain_eQTL_Genes'] = gtex_results

    os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        for cl in comment_lines:
            f.write(cl)
    df.to_csv(output_path, sep="\t", index=False, mode="a")

    print(f"\nAnnotation complete! Full results written to: {output_path}")

    # Prompt user or evaluate flag for creating the summarized file
    summary_path = None
    if create_summary is None:
        choice = input("\nGenerate a summarized file (SNP, trait count, matching Ensembl/GTEx genes)? [y/N]: ").strip().lower()
        should_summarize = choice in ['y', 'yes']
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
    files = [f for f in files if not f.endswith("_annotated.tsv") and not f.endswith("_summary.tsv")]

    if not files:
        print(f"No comparison files found in '{COMPARISONS_DIR}/'. Run compare_overlap_snps.py first.")
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
    parser = argparse.ArgumentParser(description="Annotate SNPs with Ensembl VEP proximal genes and GTEx brain eQTLs.")
    parser.add_argument("--input", "-i", help="Path to input TSV file containing SNP column")
    parser.add_argument("--output", "-o", help="Path for output annotated TSV file")
    parser.add_argument("--threads", "-t", type=int, default=5, help="Number of concurrent API worker threads (default: 5)")
    parser.add_argument("--summary", dest="create_summary", action="store_true", default=None,
                        help="Generate the simplified summary TSV without prompting")
    parser.add_argument("--no-summary", dest="create_summary", action="store_false",
                        help="Skip generating the simplified summary TSV without prompting")
    parser.add_argument("--cache-file", default=DEFAULT_CACHE_FILE, help=f"Path to local API JSON cache (default: {DEFAULT_CACHE_FILE})")
    parser.add_argument("--no-cache", action="store_true", help="Bypass local cache and force API queries")
    
    args = parser.parse_args()

    if args.input:
        output_file = args.output or args.input.replace(".tsv", "_annotated.tsv")
        annotate_file(args.input, output_file, create_summary=args.create_summary, max_workers=args.threads,
                      cache_path=args.cache_file, use_cache=not args.no_cache)
    else:
        interactive_mode()

if __name__ == "__main__":
    main()