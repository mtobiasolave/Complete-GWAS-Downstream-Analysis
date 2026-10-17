# Script to extract TSS positions from GENCODE GTF file and save as a lightweight TSV for downstream analysis.

import csv
import gzip
import os
import re

# Update paths to match your folder structure
GTF_INPUT = "data/annotations/gencode.v19.annotation.gtf_withproteinids"
TSV_OUTPUT = "data/annotations/gencode_v19_gene_tss.tsv"

# Regex patterns for GTF attributes
re_gene_name = re.compile(r'gene_name "([^"]+)"')
re_gene_type = re.compile(r'(?:gene_type|gene_biotype) "([^"]+)"')
re_gene_id = re.compile(r'gene_id "([^"]+)"')


def is_gzipped(filepath):
  with open(filepath, "rb") as test_f:
    return test_f.read(2) == b"\x1f\x8b"


def extract_tss(gtf_path, output_path):
  open_fn = gzip.open if is_gzipped(gtf_path) else open
  mode = "rt" if is_gzipped(gtf_path) else "r"

  records = []
  seen = set()

  print(f"Reading {gtf_path}...")
  with open_fn(gtf_path, mode, encoding="utf-8") as f:
    for line in f:
      if line.startswith("#"):
        continue

      fields = line.strip().split("\t")
      if len(fields) < 9:
        continue

      feature = fields[2]
      # Using 'transcript' level ensures all alternative promoter TSSs are captured
      if feature != "transcript":
        continue

      chrom = fields[0]  # e.g., 'chr1'
      start = int(fields[3])  # 1-based start
      end = int(fields[4])  # 1-based end
      strand = fields[6]  # '+' or '-'
      attributes = fields[8]

      # Strand-aware TSS calculation
      tss = start if strand == "+" else end

      # Extract gene symbol and biotype
      match_name = re_gene_name.search(attributes)
      match_type = re_gene_type.search(attributes)
      match_id = re_gene_id.search(attributes)

      gene_symbol = (
          match_name.group(1)
          if match_name
          else (match_id.group(1) if match_id else "UNKNOWN")
      )
      biotype = match_type.group(1) if match_type else "unknown"

      # Store unique (CHR, GENE_SYMBOL, TSS) tuples
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

  # Sort by chromosome and position
  records.sort(key=lambda x: (x["CHR"], x["TSS"], x["GENE_SYMBOL"]))

  # Write clean, lightweight reference TSV
  os.makedirs(os.path.dirname(output_path), exist_ok=True)
  with open(output_path, "w", newline="", encoding="utf-8") as out_f:
    writer = csv.DictWriter(
        out_f,
        fieldnames=["CHR", "GENE_SYMBOL", "STRAND", "TSS", "BIOTYPE"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(records)

  print(f"Extraction complete! {len(records)} unique TSS positions saved to:")
  print(f" -> {output_path}")


if __name__ == "__main__":
  extract_tss(GTF_INPUT, TSV_OUTPUT)