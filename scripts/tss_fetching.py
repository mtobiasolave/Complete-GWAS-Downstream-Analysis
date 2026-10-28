"""Extract transcript-level TSS positions from a GENCODE GTF file into a lightweight TSV.

Retained as a standalone reference preprocessing utility.
"""

import argparse
import csv
import gzip
import os
import re

from common_utils import PROJECT_ROOT, is_gzipped

GTF_INPUT = str(PROJECT_ROOT / "data" / "annotations" / "gencode.v19.annotation.gtf_withproteinids")
TSV_OUTPUT = str(PROJECT_ROOT / "data" / "annotations" / "gencode_v19_gene_tss.tsv")

RE_GENE_NAME = re.compile(r'gene_name "([^"]+)"')
RE_GENE_TYPE = re.compile(r'(?:gene_type|gene_biotype) "([^"]+)"')
RE_GENE_ID = re.compile(r'gene_id "([^"]+)"')


def extract_tss(gtf_path: str, output_path: str) -> None:
    """Parse transcript features from GTF and write unique (CHR, GENE_SYMBOL, TSS) rows."""
    gzipped = is_gzipped(gtf_path)
    open_fn = gzip.open if gzipped else open
    mode = "rt" if gzipped else "r"

    records = []
    seen = set()

    print(f"Reading {gtf_path}...")
    with open_fn(gtf_path, mode, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            tss = start if strand == "+" else end

            match_name = RE_GENE_NAME.search(attributes)
            match_type = RE_GENE_TYPE.search(attributes)
            match_id = RE_GENE_ID.search(attributes)

            gene_symbol = (
                match_name.group(1)
                if match_name
                else (match_id.group(1) if match_id else "UNKNOWN")
            )
            biotype = match_type.group(1) if match_type else "unknown"

            key = (chrom, gene_symbol, tss)
            if key not in seen:
                seen.add(key)
                records.append({
                    "CHR": chrom,
                    "GENE_SYMBOL": gene_symbol,
                    "STRAND": strand,
                    "TSS": tss,
                    "BIOTYPE": biotype,
                })

    records.sort(key=lambda x: (x["CHR"], x["TSS"], x["GENE_SYMBOL"]))

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(
            out_f,
            fieldnames=["CHR", "GENE_SYMBOL", "STRAND", "TSS", "BIOTYPE"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(records)

    print(f"Extraction complete! {len(records)} unique TSS positions saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract GENCODE TSS reference TSV from GTF.")
    parser.add_argument("--input", default=GTF_INPUT, help="Path to input GENCODE GTF.")
    parser.add_argument("--output", default=TSV_OUTPUT, help="Path to output TSS TSV.")
    args = parser.parse_args()
    extract_tss(args.input, args.output)