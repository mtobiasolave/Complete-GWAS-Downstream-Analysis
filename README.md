# Complete-GWAS-Downstream-Analysis

An automated, scalable pipeline integrating psychiatric and behavioral Genome-Wide Association Study (GWAS) summary statistics with LDSC genetic correlations and dual-modality functional epigenetic variant panels (**DNA Allele-Specific Methylation [DNA ASM]** and **RNA $m^6\text{A}$ methylation**).

The workflow identifies cross-phenotype pleiotropic hub variants, calculates candidate-space False Discovery Rate (FDR) adjustments using Storey's $q$-value or Benjamini-Hochberg (`fdr_bh`) procedures, performs proximal TSS ($\le 5\text{ kb}$) and Ensembl VEP / GTEx v10 cis-eQTL concordance analysis across 13 human brain tissues, and provides both publication-grade figures and an interactive Streamlit explorer.

---

## Workflow Overview

```text
Raw GWAS Summary Statistics (data/raw_data/*_gwas.txt)
  │
  ├──► [Step 1–4: LDSC Genetic Correlation Branch]
  │      munge_sumstats ──► ldsc_h2 / ldsc_rg ──► plot_heatmap.py & plot_genetic_correlations.py
  │      Outputs: results/figures/rg_heatmap.pdf, results/figures/rg_network.pdf
  │
  ├──► [Step 5a–5e: DNA ASM CpG-to-TSS (<5 kb) & GTEx Brain eQTL Branch]
  │      overlap_dna_asm.py ──► map_cpg_tss.py ──► lookup_gtex_brain_eqtls.py ──► compare_tss_gtex.py
  │      ──► interactive_cross_trait_asm.py
  │      Outputs: results/methylation/dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv
  │               results/methylation/dna_asm_{trait}_concordant_genes_summary.tsv
  │               results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv
  │
  └──► [Step 5 / 5f: RNA m6A Host/VEP & GTEx Brain eQTL Branch]
         overlap_methylation.py ──► compare_overlap_snps.py ──► annotate_snps_genes.py
         ──► plot_pleiotropy_heatmap.py
         Outputs: results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps*.tsv
                  results/figures/cross_trait_rna_m6a_recurrent_snps_annotated_effect_heatmap.pdf
```

## Key Features

- **Memory-Inverted O(1) Variant Extraction:** Streams multi-gigabyte GWAS summary statistics against in-memory candidate sets (`scripts/common_utils.py`), preventing out-of-memory errors on HPC clusters.
- **Vectorized Binary-Search TSS Interval Mapping:** Uses `np.searchsorted` interval bounds ($[\text{CPG\_POS} - 5000, \text{CPG\_POS} + 5000]$) in `scripts/map_cpg_tss.py` to map CpGs to GENCODE v19 TSS positions without Cartesian memory blowups.
- **Dual FDR Framework:** Supports both Benjamini-Hochberg FDR (`fdr_bh`) and Storey’s $q$-value ($\pi_0$-adjusted empirical FDR).
- **Automated Brain Tissue Annotation & Caching:** Queries Ensembl VEP and GTEx v10 cis-eQTLs across 13 human brain regions with persistent, atomic JSON caching (`results/methylation/.gtex_eqtl_cache.json` and `results/methylation/.api_cache.json`).
- **Two-Tier Architecture:** Unattended Slurm batch execution via Snakemake 9 (`profiles/slurm`) paired with on-demand interactive CLI and Streamlit exploration (`dashboard/app.py`).

---

## Directory Structure

```text
Complete-GWAS-Downstream-Analysis/
├── .gitignore
├── AGENTS.md                           # Architecture & operating rules
├── DATA_SCHEMA.md                      # Exact file & column schemas
├── README.md
├── Snakefile                           # Master Snakemake 9 workflow
├── config.yaml                         # HPC & pipeline configuration
├── envs/
│   └── snakemake9.yaml                 # Conda environment specification
├── profiles/
│   └── slurm/config.yaml               # Snakemake 9 Slurm executor profile
├── scripts/
│   ├── common_utils.py                 # Consolidated utilities (chr norm, GWAS I/O, FDR, GTEx/VEP cache, LDSC parser)
│   ├── overlap_methylation.py          # Step 5: Stream GWAS overlap for DNA_ASM & RNA_m6A SNP lists
│   ├── overlap_dna_asm.py              # Step 5a: Intersect GWAS with DNA ASM Tag & LD proxy table
│   ├── map_cpg_tss.py                  # Step 5b: Binary-search interval mapping of CpGs to GENCODE v19 TSS (<=5 kb)
│   ├── lookup_gtex_brain_eqtls.py      # Step 5c: Batched GTEx v10 brain eQTL lookup with persistent cache
│   ├── compare_tss_gtex.py             # Step 5d: Classify CpG nearest-TSS vs. GTEx eQTL gene concordance
│   ├── interactive_cross_trait_asm.py  # Step 5e: Cross-trait DNA ASM FDR & recurrence matrix builder
│   ├── run_interactive_pipeline.py     # Step 5f / Interactive unified runner (RNA_m6A, DNA_ASM, or both)
│   ├── compare_overlap_snps.py         # Multi-trait RNA m6A / overlap FDR & direction concordance filter
│   ├── annotate_snps_genes.py          # Ensembl VEP & GTEx v10 eQTL annotator for overlap tables
│   ├── plot_pleiotropy_heatmap.py      # Publication-grade cross-trait pleiotropy heatmap
│   ├── plot_heatmap.py                 # LDSC pairwise genetic correlation matrix heatmap
│   ├── plot_genetic_correlations.py    # LDSC circular genetic correlation network (FDR < 0.05)
│   ├── ldsnps_filtering.py             # Standalone preprocessing utility: Excel -> DNA_ASM_mapping_table.tsv
│   └── tss_fetching.py                 # Standalone preprocessing utility: GENCODE GTF -> gencode_v19_gene_tss.tsv
├── dashboard/
│   ├── app.py                          # Two-tab Streamlit Explorer (RNA m6A & DNA ASM)
│   └── requirements.txt                # Dashboard dependencies
├── data/
│   ├── annotations/                    # Tracked: cpg_coordinates_hg19.tsv, gencode_v19_gene_tss.tsv
│   ├── methylation/                    # Tracked: DNA_ASM_mapping_table.tsv, DNA_ASM_SNPlist.txt, RNA_m6A_SNPlist.txt
│   ├── raw_data/                       # Ignored: Raw GWAS summary statistics (*_gwas.txt)
│   └── reference/                      # Ignored: HapMap3 snplist & LDscores
└── results/                            # Ignored: Generated munged sumstats, LDSC logs, methylation tables & figures
```

---

## Installation & Computational Environments

1. **Snakemake 9 Batch Pipeline (HPC or Local):**
```bash
conda env create -f envs/snakemake9.yaml
conda activate snakemake9_env
```

2. **Interactive Dashboard (Local PC or Interactive Node):**
```bash
pip install -r dashboard/requirements.txt
```

---

## Pipeline Execution

### Tier 1: Automated Batch Pipeline (Snakemake 9)
```bash
# Dry-run DAG validation
snakemake -n all

# Full batch pipeline on Slurm HPC (Munge + LDSC + Figures + DNA ASM + RNA m6A)
snakemake --profile profiles/slurm all

# Methylation branch only (DNA ASM + RNA m6A)
snakemake --profile profiles/slurm --jobs 10 methylation

# Local multi-core execution
snakemake --cores 2 all
```

### Tier 2: Interactive Post-Processing & Exploration
```bash
# Option A: Launch the two-tab Streamlit Dashboard (RNA m6A + DNA ASM)
streamlit run dashboard/app.py

# Option B: Unified interactive CLI runner (choose RNA_m6A, DNA_ASM, or both)
python scripts/run_interactive_pipeline.py -i

# Option C: Standalone DNA ASM cross-trait interactive explorer
python scripts/interactive_cross_trait_asm.py -i
```