# Project Overview: GWAS, LDSC & Methylation Pipeline

## 1. Two-Tier Architecture (Batch HPC vs. Interactive Post-Processing)
This project uses a strict two-tier workflow. Do NOT convert interactive post-processing scripts into unattended Slurm batch jobs unless explicitly asked.

### Tier 1: Automated HPC Batch Pipeline (`Snakefile` + `sbatch`)
Managed by the master `Snakefile` and executed on a Slurm HPC cluster using target rules (`all` or `methylation`):
1. **GWAS Standardization (`munge`):** Reads raw GWAS summary statistics from `data/raw_data/`, harmonizes varied column names/formats using the custom dictionary inside the `Snakefile`, and outputs standardized `.sumstats.gz` files to `results/munged/`.
2. **Genetic Correlations (`ldsc_rg`):** Calculates pairwise LDSC genetic correlations (`results/ldsc_rg/`) and generates correlation heatmaps and network plots (`results/figures/`).
3. **Methylation Base Processing (`methylation`):** Runs the heavy methylation analysis steps and stores outputs in `results/methylation/`.

### Tier 2: On-Demand Interactive Exploration (Run Outside `sbatch`)
Executed manually by the user whenever needed; these scripts ONLY read existing output files from `results/`:
1. **Interactive Methylation Comparison Scripts:** Prompt the user interactively for filtering parameters (e.g., FDR engine, number of significant traits) and write comparison outputs to `results/methylation/comparisons/`.
2. **Streamlit Dashboard:** An interactive web app where the user selects and visualizes the generated methylation and comparison files.

## 2. Directory Structure (Relative to Project Root)
- `Snakefile`: Master workflow definition, target rules (`all`, `methylation`), and the custom GWAS column-mapping dictionary.
- `data/`: Contains `raw_data/`, `reference/`, and `methylation/` (excluded from Git/indexing).
- `results/`: Contains `munged/`, `ldsc_rg/`, `figures/`, `methylation/`, and `methylation/comparisons/` (excluded from Git/indexing).
- `logs/`: Slurm output and error logs (`logs/{rule}_%j.out` and `logs/{rule}_%j.err`).

## 3. HPC Execution Commands (Snakemake 9 Profile)
- **Environment Activation:**
  `conda activate snakemake9_env`
- **Full Batch Pipeline (`all`):**
  `snakemake --profile profiles/slurm all`
- **Methylation Batch Pipeline (`methylation`):**
  `snakemake --profile profiles/slurm --jobs 10 methylation`
- **Correlations Only (`correlations`):**
  `snakemake --profile profiles/slurm correlations`
- **Heritability Only (`heritability`):**
  `snakemake --profile profiles/slurm heritability`
- **GWAS Standardization Only (`munge_only`):**
  `snakemake --profile profiles/slurm munge_only`
- **Dry-Run Validation (`snakemake -n`) / Local PC Run (`--cores`):**
  `snakemake -n` / `snakemake --cores 2 all`

# Strict Agent Operating Rules for This Pipeline

1. **Never Read Heavy Data Files:**
   - Never run `cat`, `grep`, or file-reading tools on whole files inside `data/`, `results/`, or `.snakemake/`.
   - Always check `DATA_SCHEMA.md` first to see exact directory paths and file extensions.
   - If you must inspect a data file during testing, ONLY use `head -n 5` or `zcat <file> | head -n 5`.

2. **Relative Paths Only:**
   - All rules in `Snakefile` and scripts in `scripts/` must use relative paths (`data/...`, `results/...`, `logs/...`) or `config.yaml` keys so the repository runs identically on this local PC and on the Slurm HPC cluster.

3. **Dry-Run First:**
   - When modifying `Snakefile` rules, validate syntax and DAG resolution using `snakemake -n` before attempting to execute any scripts.

4. **GWAS Header Dictionary:** 
   - The master `Snakefile` already contains a custom dictionary mapping the specific column headers (`SNP`, `A1`, `A2`, `P`, `N`, etc.) and flags for each GWAS summary statistics dataset in `data/raw_data/`.
   - Always respect and use this existing dictionary as the source of truth when modifying rules, adding new traits, or generating local mock test files.