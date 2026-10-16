# Complete-GWAS-Downstream-Analysis

An automated, scalable pipeline integrating psychiatric and behavioral Genome-Wide Association Study (GWAS) summary statistics with functional epigenetic variant panels (RNA $m^6A$ methylation and DNA allele-specific methylation [DNA ASM]).

The workflow identifies cross-phenotype pleiotropic hub variants, calculates candidate-space False Discovery Rate (FDR) adjustments using Storey's $q$-value or Benjamini-Hochberg procedures, performs functional and positional gene annotation across 13 human brain tissues, and generates publication-grade visualizations.

---

## Workflow Overview

```text
Raw GWAS Sumstats + Candidate Methylation Panels (DNA ASM / RNA m6A)
                              │
                              ▼
        [Snakemake Rule: methylation_overlap]
        Memory-inverted awk extraction with dynamic header identification
                              │
                              ▼
        Candidate Overlap Tables (results/methylation/*_overlap_fullheaders.txt)
                              │
                              ▼
        [Script: scripts/compare_overlap_snps.py]
        Trait-specific FDR (Storey q / BH <= 0.05), direction harmonization, recurrence filter (k >= 5)
                              │
                              ▼
        Cross-Trait Candidate TSV (results/methylation/comparisons/common_*_sig_snps.tsv)
                              │
                              ▼
        [Script: scripts/annotate_snps_genes.py]
        REST API queries: Ensembl VEP (proximal) & GTEx v10 (13 brain tissue cis-eQTLs)
                              │
                              ▼
        Annotated & Summarized TSVs (*_annotated.tsv & *_summary.tsv)
                              │
                              ▼
        [Script: scripts/plot_pleiotropy_heatmap.py]
        Log-odds harmonization (Beta / ln(OR)), 90th-percentile scaling, publication plots
                              │
                              ▼
        Publication Figures (results/figures/*_heatmap.png / .pdf)

Key Features

Memory-Inverted Variant Extraction: Streams multi-gigabyte summary statistics against in-memory candidate hash tables using Awk, preventing walltime limits and out-of-memory errors on HPC clusters.

Dual FDR Framework: Supports both conservative Benjamini-Hochberg FDR (fdr_bh) and Storey’s $q$-value (replicating Bioconductor qvalue in R via cubic smoothing splines for $\pi_0$ estimation).

Automated Brain Tissue Annotation: Directly interfaces with external REST APIs to extract Ensembl VEP transcript consequences and GTEx v10 cis-eQTL target genes across 13 distinct brain tissues.

Effector Prioritization: Automatically filters for convergent effector genes (genes shared between positional and eQTL annotations) while reporting divergent loci.

Harmonized Effect Visualization: Scales binary and quantitative GWAS traits onto a unified scale ($\beta_{\text{equivalent}} = \ln(\text{OR})$) with robust outlier handling.Directory Structure

Plaintext

Complete-GWAS-Downstream-Analysis/
├── .gitignore
├── README.md
├── Snakefile
├── config.yaml                      # Externalized HPC paths and reference files
├── envs/
│   └── snakemake9.yaml              # Snakemake 9 & Slurm plugin Conda environment
├── profiles/
│   └── slurm/config.yaml            # Snakemake 9 HPC Slurm profile configuration
├── scripts/
│   ├── run_interactive_pipeline.py  # End-to-end chained runner (compare -> annotate -> plot)
│   ├── compare_overlap_snps.py      # Trait comparison, Storey/BH FDR, pleiotropy filter
│   ├── annotate_snps_genes.py       # REST API integration: Ensembl VEP & GTEx v10 eQTLs + caching
│   ├── plot_pleiotropy_heatmap.py   # Publication-grade diverging heatmap visualization
│   ├── plot_genetic_correlations.py # Circular genetic correlation network
│   └── plot_heatmap.py              # LDSC genetic correlation matrix heatmap
├── dashboard/
│   ├── app.py                       # Streamlit Visual Explorer & in-app generator
│   └── requirements.txt             # Dashboard dependencies
├── data/                            # Ignored by git
│   ├── raw_data/                    # Raw GWAS summary statistics (.txt, .tsv)
│   ├── reference/                   # HapMap3 reference snplist
│   └── methylation/                 # Candidate variant lists (m6A, ASM)
└── results/                         # Ignored by git
    ├── munged/                      # Standardized .sumstats.gz outputs
    ├── ldsc_rg/                     # Bivariate LDSC genetic correlation logs
    ├── methylation/                 # Per-trait extracted overlap outputs
    │   └── comparisons/             # Multi-trait recurrent TSVs and annotations
    └── figures/                     # 300 DPI PNG and vector PDF outputs

Installation & Computational Environments

1. Snakemake 9 Batch Pipeline (HPC):
```bash
conda env create -f envs/snakemake9.yaml
conda activate snakemake9_env
```

2. Interactive Dashboard (Local PC or Interactive Node):
```bash
pip install -r dashboard/requirements.txt
```

Pipeline Execution

### Tier 1: Automated Slurm Batch Pipeline (Snakemake 9)
```bash
# Dry run verification
snakemake -n

# Full pipeline: Munge -> Genetic Correlations -> Heatmap & Network figures
snakemake --profile profiles/slurm all

# Candidate methylation overlap extraction
snakemake --profile profiles/slurm --jobs 10 methylation
```

### Tier 2: Interactive Post-Processing & Exploration

Option A: Streamlit Web Dashboard (In-App Generator + Visual Explorer + Figure Export)
```bash
streamlit run dashboard/app.py
```

Option B: End-to-End Guided Interactive CLI Runner
```bash
python scripts/run_interactive_pipeline.py -i
```

Option C: Step-by-Step CLI Execution
```bash
# 1. Multi-Trait Aggregation & FDR Filtering
python scripts/compare_overlap_snps.py -i

# 2. Functional & eQTL Brain Annotation (with persistent API caching)
python scripts/annotate_snps_genes.py --summary

# 3. Publication Figure Generation (300 DPI PNG & Vector PDF)
python scripts/plot_pleiotropy_heatmap.py

Data Availability & External Databases

Ensembl VEP REST API: https://rest.ensembl.org/
GTEx Portal API v2: https://gtexportal.org/api/v2/ (Dataset: gtex_v10)
Target GTEx Brain Tissues (13 Regions): Amygdala, Anterior cingulate cortex (BA24), Caudate, Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex (BA9), Hippocampus, Hypothalamus, Nucleus accumbens, Putamen, Cervical spinal cord, and Substantia nigra.