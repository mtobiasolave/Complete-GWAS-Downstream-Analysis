# Project Overview: GWAS x DNA & RNA Methylation Pipeline

## 1. Architecture
- **Batch pipeline (`Snakefile`):** Two tracks over user-selected GWAS summary statistics in `data/raw_data/`.
  - Track A (DNA methylation): A1 overlap DNA ASM tag SNPs with GWAS + BH FDR 5% (tags absent from the GWAS tested via their highest-LD proxy), then expand significant tags to all their LD proxies (kept even if absent from the GWAS) -> A2 CpG-to-gene TSS (<= 5 kb) -> A3 GTEx v10 brain eQTLs (13 tissues, BH-significant SNPs only) -> A4 gene concordance -> A5 pleiotropy table `results/dna/pleiotropy_dna.tsv`.
  - Track B (RNA methylation): B1 overlap RNA m6A SNP list with GWAS + BH FDR 5% -> B2 Ensembl VEP (GRCh37 REST, `rna.ensembl_server`) genes each BH-significant SNP falls inside, one job for all traits (`results/rna/rna_sig_genes.tsv`; up/downstream genes excluded) -> B3 pleiotropy table `results/rna/pleiotropy_rna.tsv` (`SNP, GENES, N_TRAITS, TRAITS_PRESENT, MIN_FDR_BH`).
  - Munging: `munge` target runs LDSC `../ldsc/munge_sumstats.py` -> `results/munged/{trait}.sumstats.gz` (sample size: `ldsc.n_override` > `ldsc.n_col_override` > detected N column).
  - Genetic correlation: `rg` target runs `../ldsc/ldsc.py --rg` once per trait pair (alphabetical, `results/ldsc_rg/{trait1}__{trait2}.log`; LD scores `ldsc.ld_ref_chr`/`ldsc.w_ld_chr`) for config `rg: {traits, target}` (empty target = all-vs-all) -> `results/ldsc_rg/rg_summary.tsv`. Pair selection and log parsing live in `scripts/ldsc_rg_summary.py`, shared with `app.py`.
  - Targets: `all`, `dna`, `rna`, `munge`, `rg`. Traits: `config.yaml`/`config_local.yaml` `traits`, or `--config traits=A,B`.
- **Front-end (`app.py`):** Streamlit app that launches runs in the background (local process via `scripts/run_pipeline.py`, or Slurm controller job via `scripts/run_pipeline.sbatch`), explores the pleiotropy tables with cached, in-memory filters, and runs/visualizes genetic correlations (Genetic Correlation tab; BH FDR across the pairs on screen).

## 2. Directory Structure (Relative to Project Root)
- `Snakefile`, `config.yaml`, `profiles/slurm/config.yaml`, `envs/snakemake9.yaml`, `app.py`
- `scripts/`: pipeline steps and shared helpers (`common_utils.py`).
- `data/`: `raw_data/` (GWAS), `methylation/` (reference SNP lists), `annotations/` (CpG coordinates, GENCODE TSS).
- `results/`: `dna/`, `rna/`, `munged/`, `ldsc_rg/`, `.cache/` (git-ignored).
- `patches/ldsc_python3.patch`: Python 3 fixes for the sibling `../ldsc` repository.
- `logs/`: `runs/` (app and controller logs), `slurm/` (per-job Slurm logs) (git-ignored).

## 3. Execution Commands (Snakemake 9)
- Environment: `conda activate snakemake9_env`
- Dry run: `snakemake -n all`
- Local: `snakemake --cores 2 all [--config traits=A,B]`
- Slurm: `snakemake --profile profiles/slurm all` or `sbatch scripts/run_pipeline.sbatch all A,B`
- App: `streamlit run app.py`

# Strict Agent Operating Rules for This Pipeline

1. **Never Read Heavy Data Files:**
   - Never run `cat`, `grep`, or file-reading tools on whole files inside `data/`, `results/`, or `.snakemake/`.
   - If you must inspect a data file during testing, ONLY use `head -n 5` or `zcat <file> | head -n 5`.

2. **Relative Paths Only:**
   - No absolute paths in any file. Snakefile paths are relative to the repository root (`workdir: workflow.basedir`); script defaults resolve from `common_utils.PROJECT_ROOT`; `app.py` resolves from its own location.

3. **Multiple Testing:**
   - Benjamini-Hochberg FDR (`common_utils.bh_fdr`) at `fdr_alpha` (default 0.05) is the only correction method, applied per trait at the overlap step (A1/B1). A1's testing family is one SNP per DNA ASM tag (not tags + proxies).

4. **Dry-Run First:**
   - When modifying `Snakefile` rules, validate syntax and DAG resolution using `snakemake -n` before executing.

5. **Slurm Profile:**
   - Do not put `--output`, `--error`, `--partition`, `--mem` or `--time` in `slurm_extra`; the executor plugin forbids them. Use resources and `slurm-logdir` instead.
