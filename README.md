# Complete-GWAS-Downstream-Analysis

An automated Snakemake pipeline that overlaps psychiatric and behavioral GWAS summary statistics with two functional methylation variant panels, **DNA allele-specific methylation (DNA ASM)** and **RNA m6A methylation**, and identifies cross-trait pleiotropic loci. A Streamlit app (`app.py`) configures and launches runs (locally or on Slurm) and explores the results.

Multiple testing is controlled with **Benjamini-Hochberg FDR at 5%**, applied per trait. DNA ASM: one test per tag SNP (a tag absent from the GWAS is tested via its highest-LD proxy, `dna.tag_fallback`), and every LD proxy of a significant tag is retained, even if absent from the GWAS. RNA m6A: across the reference SNPs found in that trait's GWAS.

---

## Workflow Overview

```text
GWAS summary statistics (data/raw_data/*_gwas.txt)
  │
  ├──► Track A: DNA Methylation (DNA ASM tag SNPs & LD proxies with CpG annotations)
  │      A1 overlap_dna_asm.py         Tag SNPs x GWAS + BH FDR 5%, then expand significant tags to LD proxies
  │      A2 map_cpg_tss.py             CpG -> genes with TSS within 5 kb (hg19, GENCODE v19)
  │      A3 lookup_gtex_brain_eqtls.py GTEx v10 eQTL target genes, 13 brain regions (significant SNPs only)
  │      A4 compare_tss_gtex.py        Concordance: CpG-proximal gene == GTEx brain eQTL gene
  │      A5 pleiotropy.py --mode dna   results/dna/pleiotropy_dna.tsv
  │
  ├──► Track B: RNA Methylation (RNA m6A SNP list)
  │      B1 overlap_rna_m6a.py         Overlap reference SNPs with GWAS + BH FDR 5%
  │      B2 annotate_rna_genes.py      Ensembl GRCh37 genes each significant SNP falls inside (VEP REST)
  │      B3 pleiotropy.py --mode rna   results/rna/pleiotropy_rna.tsv
  │
  ├──► Munging (target `munge`, independent of both tracks)
  │      ../ldsc/munge_sumstats.py     results/munged/{trait}.sumstats.gz (+ LDSC .log)
  │
  ├──► Genetic correlation (target `rg`, from munged sumstats)
  │      ../ldsc/ldsc.py --rg          results/ldsc_rg/{trait1}__{trait2}.log (one job per pair)
  │      ldsc_rg_summary.py            results/ldsc_rg/rg_summary.tsv (rg, SE, Z, P, BH FDR across the run's pairs)
  │
  └──► SNP heritability (target `h2`, from munged sumstats)
         ../ldsc/ldsc.py --h2          results/ldsc_h2/{trait}.log (one job per trait; observed or liability scale)
         ldsc_h2_summary.py            results/ldsc_h2/h2_summary.tsv (observed & liability h2, SE, Z, P, intercept, ratio)
```

A SNP is *present* in a trait when it is BH-significant in that trait. The pleiotropy tables keep every SNP significant in at least one trait; the trait-sharing threshold is applied interactively in the app.

| Output | Columns |
|---|---|
| `results/dna/pleiotropy_dna.tsv` | `SNP, N_TRAITS, TRAITS_PRESENT, CONCORDANT, CONCORDANT_GENES, CPG_GENES_5KB, GTEX_BRAIN_GENES, CPG_IDS, GTEX_TISSUES, MIN_FDR_BH` |
| `results/rna/pleiotropy_rna.tsv` | `SNP, GENES, N_TRAITS, TRAITS_PRESENT, MIN_FDR_BH` (`GENES`: overlapping Ensembl genes, or `intergenic` / `not in Ensembl` / `no primary-assembly locus`) |
| `results/rna/rna_sig_genes.tsv` | One row per BH-significant RNA m6A SNP: `SNP, GENES, ENSEMBL_GENE_IDS, STATUS, CHR, POS` (GRCh37), `CONSEQUENCES, TRANSCRIPT_BIOTYPES, ENSEMBL_RELEASE` |
| `results/dna/{trait}_overlap_all.tsv` | Every mapping row (tag + LD proxies) of each tested tag with the tag-level `GWAS_P, GWAS_BETA, GWAS_FDR_BH, BH_SIGNIFICANT, TESTED_SNP` and row-level `LD_R2, QUERY_IN_GWAS, QUERY_GWAS_P, QUERY_GWAS_BETA` |
| `results/rna/{trait}_overlap_bh.tsv` | Every overlapping RNA m6A SNP with the same BH columns |
| `results/ldsc_rg/rg_summary.tsv` | `Trait 1, Trait 2, rg, SE, Z, P-value, FDR`, h2 per trait, `gcov_int`, `Note` (why rg is NA) |
| `results/ldsc_h2/h2_summary.tsv` | `Trait, Scale, h2 (observed), SE (observed), h2 (liability), SE (liability), Z, P-value` (one-sided, h2 > 0), `Pop prev (K), Samp prev (P), Intercept, Intercept SE, Lambda GC, Mean chi2, Ratio, Ratio SE, SNPs, Note` |

---

## Directory Structure

All paths are relative to the repository root, so the folder can be moved between a workstation and the HPC unchanged.

```text
├── app.py                         # Streamlit front-end: run configuration, live progress, result explorer
├── Snakefile                      # Workflow: targets all / dna / rna
├── config.yaml                    # Pipeline configuration (HPC defaults)
├── config_local.yaml              # Optional, git-ignored local overrides (e.g. reduced trait subset)
├── envs/snakemake9.yaml           # Conda environment (pipeline + app)
├── profiles/slurm/config.yaml     # Snakemake 9 Slurm executor profile (partition irbio01)
├── patches/ldsc_python3.patch     # Python 3 fixes for the sibling ../ldsc repository (CBIIT/ldsc 6c67395)
├── scripts/
│   ├── common_utils.py            # GWAS discovery & streaming, BH FDR, HTTP/cache helpers
│   ├── overlap_dna_asm.py         # A1
│   ├── map_cpg_tss.py             # A2
│   ├── lookup_gtex_brain_eqtls.py # A3
│   ├── compare_tss_gtex.py        # A4
│   ├── overlap_rna_m6a.py         # B1
│   ├── annotate_rna_genes.py      # B2
│   ├── pleiotropy.py              # A5 / B3
│   ├── ldsc_rg_summary.py         # LDSC --rg log parser, pair selection, BH FDR across pairs
│   ├── ldsc_h2_summary.py         # LDSC --h2 flags & validation, log parser (observed/liability), summary
│   ├── run_pipeline.py            # Launcher used by app.py and run_pipeline.sbatch
│   ├── run_pipeline.sbatch        # Slurm controller job
│   ├── ldsnps_filtering.py        # Reference builder: Excel -> DNA_ASM_mapping_table.tsv
│   └── tss_fetching.py            # Reference builder: GENCODE GTF -> gencode_v19_gene_tss.tsv
├── data/
│   ├── annotations/               # cpg_coordinates_hg19.tsv, gencode_v19_gene_tss.tsv
│   ├── methylation/               # DNA_ASM_mapping_table.tsv, RNA_m6A_SNPlist.txt
│   └── raw_data/                  # GWAS summary statistics (git-ignored)
├── results/                       # dna/, rna/, munged/, ldsc_rg/ & ldsc_h2/ (+ stdout/ live logs), .cache/ (git-ignored)
└── logs/                          # runs/ (app & controller logs), munge/ (live munge output), slurm/ (git-ignored)
../ldsc/                           # Sibling LDSC repository providing munge_sumstats.py and ldsc.py
```

---

## Installation

```bash
git clone git@github.com:mtobiasolave/Complete-GWAS-Downstream-Analysis.git
cd Complete-GWAS-Downstream-Analysis

# Strict channel priority keeps the Anaconda defaults channel out of the solve
CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/snakemake9.yaml
conda activate snakemake9_env
conda list | grep -c pkgs/main   # expect 0
```
`envs/snakemake9.yaml` pins the versions the pipeline was tested with (Snakemake 9.24, pandas 3.0, Python 3.12). If `pkgs/main` packages still appear, create the environment with `conda create -n snakemake9_env --override-channels -c conda-forge -c bioconda --strict-channel-priority` and the package list from the YAML.

The `munge`, `rg` and `h2` targets also need the sibling LDSC repository (see [LDSC requirements](#ldsc-requirements)) and `data/reference/w_hm3.snplist` plus `data/reference/LDscores/` (git-ignored; copy them separately).

---

## Usage

### Streamlit app (local or HPC login node)
```bash
streamlit run app.py
```
**Methylation tabs:** select traits, target mode (Both / DNA only / RNA only) and execution engine (Local or HPC; HPC is preselected when `sbatch` is available), then press **Run pipeline**. The run continues in the background; the app shows live progress and the log, and refreshes the result tables when it finishes. The RNA tab lists each SNP's overlapping Ensembl gene(s) and can be searched by rsID or gene. On the HPC, open the app through SSH port forwarding, e.g. `ssh -L 8502:localhost:8502 user@cluster`.

**🧹 Munge Sumstats tab:** select unmunged traits (or tick *Include already munged traits* to re-munge), check the reference alleles file and the per-trait sample size (manual N > N column override > detected N column; the resolved LDSC flags are shown), then press **Run Munging**. Each trait shows Queued / Running / Done with a live log viewer, and the inventory lists `results/munged/` with size, date, SNPs with Z, mean χ² and λGC, plus file previews.

**📊 Genetic Correlation tab:** check the munged-sumstats QC table (SNPs retained, mean χ², λGC), choose **All-vs-All Matrix** (N ≥ 2 traits, N(N−1)/2 pairs) or **Target vs. Panel** (one target, K panel traits), then press **Run Genetic Correlation**. Already computed pairs are reused (tick *Re-run* to recompute). A status panel shows per-pair progress and a live LDSC log. Results show a symmetric rg heatmap (matrix) or a forest plot with 95% CI (target vs. panel), plus a table with BH FDR computed across the pairs on screen (`*` = FDR < 5%) and a TSV export.

**🌱 Heritability (h²) tab:** select one or more munged traits (*Select all* for a batch). For case-control traits, switch on **Liability scale** and enter each trait's population and sample prevalence in the per-trait table (see [Heritability: liability scale](#heritability-h2-liability-scale)); traits left empty stay on the observed scale. **⚙️ Advanced options** shows the LD reference in use and sets the intercept: *Free* (LDSC default), *Constrained to 1*, or a *Fixed value*. The tab counts the traits already computed with exactly these settings and only re-runs the rest (tick *Re-run* to force). Results offer two views of the selected traits, plus an **Observed / Liability** scale switch when liability results exist:
- **Multi-Trait Comparison (Panel):** a forest plot of h² with 95% CI or ± 1 SE error bars, sortable by h², trait name or Z; traits with BH FDR < 5% across the traits shown are highlighted, h² outside 0–1 is drawn as an open diamond, and non-converged / invalid-SE traits are listed instead of plotted. Below it, the comparative table (with FDR) and TSV export.
- **Single Trait Deep-Dive (Target):** metric cards (h², SE, Z, one-sided P, intercept), the trait's 95% CI highlighted among the other selected traits, an **LDSC diagnostic checklist** (mean χ² ≥ 1.02, intercept ≈ 1 or lower if GC-corrected, ratio ≈ 0 with 10–20% not uncommon, ≥ 200,000 regression SNPs, 0 ≤ h² ≤ 1; benchmarks quoted from the LDSC wiki, which defines no minimum h² Z-score), and the raw LDSC log.

Results computed with different settings than selected are flagged.

### Command line
```bash
# Dry-run DAG validation
snakemake -n all

# Local run (traits from config / config_local.yaml, or explicit)
snakemake --cores 2 all
snakemake --cores 2 dna --config traits=SUD2023,BipolarDisorder2024

# LDSC munging (all traits, or a subset; manual N per trait lives in config.yaml ldsc.n_override)
snakemake --cores 2 munge
snakemake --cores 2 munge --config "traits=SUD2023" --forcerun munge_sumstats

# Genetic correlation: all pairs of rg.traits, or rg.target vs each of rg.traits
snakemake --cores 2 rg --config traits=A,B,C 'rg={"traits": ["A", "B", "C"], "target": ""}'
snakemake --cores 2 rg --config traits=A,B,C 'rg={"traits": ["B", "C"], "target": "A"}'

# SNP heritability: observed scale for all h2.traits ...
snakemake --cores 2 h2 --config traits=A,B 'h2={"traits": ["A", "B"]}'
# ... A on the liability scale (both prevalences), B observed, intercept fixed at 1
snakemake --cores 2 h2 --config traits=A,B 'h2={"traits": ["A", "B"], "intercept_h2": 1, "prevalence": {"A": {"pop_prev": 0.01, "samp_prev": 0.5}}}'

# Slurm, directly from the login node
snakemake --profile profiles/slurm all

# Slurm, as a controller job that survives logout
mkdir -p logs/runs
sbatch scripts/run_pipeline.sbatch all SUD2023,BipolarDisorder2024
```

### LDSC requirements
Munging runs `../ldsc/munge_sumstats.py` with the `snakemake9_env` interpreter (validated with pandas 3: Z, alleles and N match the raw sumstats). The sibling LDSC copy must contain the Python 3 fixes; on a fresh machine, from the repository root:
```bash
git clone https://github.com/CBIIT/ldsc.git ../ldsc
git -C ../ldsc checkout 6c67395      # CBIIT/ldsc base commit the patch applies to
git -C ../ldsc apply --check "$PWD/patches/ldsc_python3.patch" && git -C ../ldsc apply "$PWD/patches/ldsc_python3.patch"
```
(`git -C` resolves relative paths from `../ldsc`, hence `$PWD`.)

### Heritability (h²): liability scale
LDSC estimates SNP heritability on the **observed scale**. For a case-control trait it can be transformed to the **liability scale** with two prevalences, which must be given together (both or neither) and each strictly between 0 and 1 (enforced by the app, the Snakefile and LDSC itself):

- `--pop-prev` (**K**): population (lifetime) prevalence of the disorder.
- `--samp-prev` (**P**): proportion of cases in the GWAS sample.
  - Use **P = 0.5** when the trait was munged with an **effective sample size** (sum of per-cohort effective N, e.g. a `Neff` / `N_EFF` column), as recommended by Grotzinger et al. (2023, *Biological Psychiatry*) for meta-analysed GWAS.
  - Otherwise (total N or a manual N of cases + controls) use the actual case fraction, cases / (cases + controls).
  - The tab shows which N each trait was munged with (from its munge log), pre-fills P = 0.5 for effective-N traits, and warns when P = 0.5 is used with a total N, or when a half effective N (`Neff_half`) was used (which would inflate h² about 2×).

LDSC prints either the observed or the liability estimate; for liability runs the observed-scale h² is recovered exactly with LDSC's own conversion factor. Z = h² / SE is identical on both scales and the P-value is one-sided (h² > 0). LDSC reports h² and SE with 4 decimals, so Z and P carry that rounding (relevant only for weak signals, Z ≈ 2).

### Configuration fields
| File | Field | Purpose |
|---|---|---|
| `config.yaml` | `h2.traits` | Traits for the `h2` target (set per run by the app or `--config`) |
| `config.yaml` | `h2.intercept_h2` | Empty = free intercept; a positive value is passed as `--intercept-h2` (1 = constrained) |
| `config.yaml` | `h2.prevalence` | Optional per-trait defaults `{trait: {pop_prev, samp_prev}}`, pre-filled in the app's liability table |
| `profiles/slurm/config.yaml` | `set-resources.ldsc_h2` | HPC runtime (30 min) and memory (8 GB) per h2 job |
| `config_local.yaml` | — | No new keys needed: the LD reference (`ldsc.ld_ref_chr` / `ldsc.w_ld_chr`) is the same relative path on both machines. Any `h2.*` key set here overrides `config.yaml` locally. |

### Notes
- Step A3 queries the GTEx Portal API and needs outbound internet access on the node that runs it. Results are cached in `results/.cache/gtex_eqtl_cache.json`; failed requests are not cached and make the job fail so it can be retried.
- Step B2 queries the Ensembl VEP REST API (`rna.ensembl_server`, GRCh37 by default to match Track A's hg19 coordinates) and also needs internet access. A gene is reported when the SNP falls inside one of its transcripts (any biotype); genes VEP reports only as up/downstream (≤ 5 kb) are excluded. Results are cached in `results/.cache/ensembl_grch37_genes_cache.json`; failed requests are not cached and make the job fail. Use a different cache file if you switch server/assembly.
- `config_local.yaml` overrides `config.yaml` whenever it is present; do not copy it to the HPC unless you want its trait subset there.
- Always run the pipeline and the app from the **activated** `snakemake9_env` (`conda activate snakemake9_env` or `conda run -n snakemake9_env ...`). Calling the environment's `python.exe` directly on Windows leaves its DLL folder off `PATH`, and LDSC regressions (`rg`, `h2`) then crash in NumPy's linear algebra with `Windows fatal exception: code 0xc06d007f`.
