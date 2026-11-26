# ==============================================================================
# DNA & RNA Methylation x GWAS Pipeline (VIBRA HPC)
# Compatible with Snakemake 9 & Slurm Executor Plugin
#
# Targets:
#   all : Track A (DNA methylation) + Track B (RNA methylation)
#   dna : Track A only -> results/dna/pleiotropy_dna.tsv
#   rna : Track B only -> results/rna/pleiotropy_rna.tsv (with overlapping hg19/GRCh37 genes; config rna_gene_source)
#   munge : LDSC munge_sumstats.py -> results/munged/{trait}.sumstats.gz
#   rg    : LDSC genetic correlation for the pairs in config `rg`
#           -> results/ldsc_rg/{trait1}__{trait2}.log + results/ldsc_rg/rg_summary.tsv
#   h2    : LDSC SNP heritability for the traits in config `h2` (observed or liability scale)
#           -> results/ldsc_h2/{h2_trait}.log + results/ldsc_h2/h2_summary.tsv
#
# Local:  snakemake --cores 2 all [--config traits=TraitA,TraitB]
# Slurm:  snakemake --profile profiles/slurm all
# rg:     snakemake --cores 2 rg --config traits=A,B,C 'rg={"traits": ["A", "B", "C"], "target": ""}'
# h2:     snakemake --cores 2 h2 --config traits=A,B 'h2={"traits": ["A", "B"], "prevalence": {"A": {"pop_prev": 0.01, "samp_prev": 0.5}}}'
# ==============================================================================

import re
import sys
from pathlib import Path

# Resolve every relative path below from the repository root, wherever snakemake is launched
workdir: workflow.basedir

sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from common_utils import discover_gwas, munge_column_flags, read_gwas_header
from ldsc_rg_summary import pair_stem, rg_pairs
from ldsc_h2_summary import h2_flags

# Load configuration file (merges local overrides if config_local.yaml exists)
configfile: "config.yaml"
if Path("config_local.yaml").exists():
    configfile: "config_local.yaml"

# Run scripts with the interpreter running Snakemake (the activated snakemake9_env);
# Slurm jobs are submitted with --export=ALL and resolve the same environment
PYTHON = f'"{Path(sys.executable).as_posix()}"'
FDR_ALPHA = float(config.get("fdr_alpha", 0.05))
LDSC = config.get("ldsc") or {}

# ----- TRAIT DETECTION ----- #
# 1. Discover all GWAS summary statistics in raw_data_dir, plus any extra configured paths
RAW_DATA_DIR = config.get("raw_data_dir", "data/raw_data")
GWAS_FILES = discover_gwas(RAW_DATA_DIR)
GWAS_FILES.update(config.get("gwas_sumstats") or {})

# 2. Restrict to the requested traits (list in config, or "A,B,C" from --config traits=A,B,C)
def requested_traits():
    selected = config.get("traits") or []
    if isinstance(selected, str):
        selected = selected.split(",")
    selected = sorted({str(t).strip() for t in selected if str(t).strip()})
    return selected or sorted(GWAS_FILES)

TRAITS = requested_traits()
MISSING_TRAITS = [t for t in TRAITS if t not in GWAS_FILES]
if MISSING_TRAITS:
    raise ValueError(
        f"No GWAS file found for trait(s) {MISSING_TRAITS} in '{RAW_DATA_DIR}/'. "
        f"Available traits: {sorted(GWAS_FILES)}"
    )

# 3. Genetic correlation pairs: all-vs-all of rg.traits, or rg.target vs each of rg.traits
RG = config.get("rg") or {}
RG_PAIRS = rg_pairs(RG.get("traits") or [], RG.get("target") or None)

# 4. SNP heritability traits: config h2.traits, with per-trait prevalences (liability scale)
#    and an optional fixed intercept; flags are validated when the DAG is built
H2 = config.get("h2") or {}
H2_TRAITS = sorted({str(t).strip() for t in (H2.get("traits") or []) if str(t).strip()})

wildcard_constraints:
    trait="|".join(re.escape(t) for t in TRAITS) or "[^/]+",
    # Pair outputs are {trait1}__{trait2}; trait names must not contain "__"
    trait1="(?:(?!__)[^/])+",
    trait2="(?:(?!__)[^/])+",
    h2_trait="[^/]+"
# -----                ----- #

localrules: all, dna, rna, munge, rg, rg_summary, h2, h2_summary, pleiotropy_dna, pleiotropy_rna

# ==========================================
# PIPELINE TARGETS
# ==========================================

rule all:
    input:
        "results/dna/pleiotropy_dna.tsv",
        "results/rna/pleiotropy_rna.tsv"

rule dna:
    input:
        "results/dna/pleiotropy_dna.tsv"

rule rna:
    input:
        "results/rna/pleiotropy_rna.tsv"

rule munge:
    input:
        expand("results/munged/{trait}.sumstats.gz", trait=TRAITS)

rule rg:
    input:
        "results/ldsc_rg/rg_summary.tsv"

rule h2:
    input:
        "results/ldsc_h2/h2_summary.tsv"

# ------------------------------------------------------------------------------
# Munge GWAS Summary Statistics with LDSC munge_sumstats.py (Per-Trait)
# ------------------------------------------------------------------------------
def munge_flags(wildcards, input):
    """Column flags from the GWAS header plus the trait's manual N / N column overrides."""
    return " ".join(munge_column_flags(
        read_gwas_header(input.sumstats),
        n_override=(LDSC.get("n_override") or {}).get(wildcards.trait),
        n_col_override=(LDSC.get("n_col_override") or {}).get(wildcards.trait),
    ))

rule munge_sumstats:
    input:
        sumstats=lambda wildcards: GWAS_FILES[wildcards.trait],
        # ancient(): updating the LDSC repo or reference must not force re-munging
        merge_alleles=ancient(LDSC.get("merge_alleles", "data/reference/w_hm3.snplist")),
        script=ancient(LDSC.get("munge_script", "../ldsc/munge_sumstats.py"))
    output:
        "results/munged/{trait}.sumstats.gz"
    log:
        # LDSC writes {out}.log itself (buffered); stdout carries the same messages live
        ldsc="results/munged/{trait}.log",
        stdout="logs/munge/{trait}.log"
    params:
        out_prefix="results/munged/{trait}",
        col_flags=munge_flags,
        chunksize=LDSC.get("chunksize", 500000)
    shell:
        """
        {PYTHON} -u {input.script}             --sumstats {input.sumstats}             --merge-alleles {input.merge_alleles}             --out {params.out_prefix}             {params.col_flags}             --chunksize {params.chunksize}             > {log.stdout} 2>&1
        """

# ------------------------------------------------------------------------------
# Genetic Correlation (rg) with LDSC ldsc.py --rg (Per Trait Pair)
# ------------------------------------------------------------------------------
rule ldsc_rg:
    input:
        sumstats1="results/munged/{trait1}.sumstats.gz",
        sumstats2="results/munged/{trait2}.sumstats.gz",
        # ancient(): updating the LDSC repo must not force re-running every pair
        script=ancient(LDSC.get("ldsc_script", "../ldsc/ldsc.py"))
    output:
        # LDSC writes {out}.log itself (buffered); stdout carries the same messages live
        "results/ldsc_rg/{trait1}__{trait2}.log"
    log:
        "results/ldsc_rg/stdout/{trait1}__{trait2}.log"
    params:
        out_prefix="results/ldsc_rg/{trait1}__{trait2}",
        ld_ref=LDSC.get("ld_ref_chr", "data/reference/LDscores/"),
        w_ld=LDSC.get("w_ld_chr", "data/reference/LDscores/")
    shell:
        """
        {PYTHON} -u {input.script}             --rg {input.sumstats1},{input.sumstats2}             --ref-ld-chr {params.ld_ref}             --w-ld-chr {params.w_ld}             --out {params.out_prefix}             > {log} 2>&1
        """

def rg_summary_inputs(wildcards):
    if not RG_PAIRS:
        raise ValueError(
            "No genetic correlation pairs: set config rg.traits (>= 2 traits, or rg.target plus >= 1 trait)."
        )
    return [f"results/ldsc_rg/{pair_stem(*pair)}.log" for pair in RG_PAIRS]

rule rg_summary:
    input:
        rg_summary_inputs
    output:
        "results/ldsc_rg/rg_summary.tsv"
    shell:
        """
        {PYTHON} scripts/ldsc_rg_summary.py --logs {input} --alpha {FDR_ALPHA} --output {output}
        """

# ------------------------------------------------------------------------------
# SNP Heritability (h2) with LDSC ldsc.py --h2 (Per Trait)
# ------------------------------------------------------------------------------
rule ldsc_h2:
    input:
        sumstats="results/munged/{h2_trait}.sumstats.gz",
        # ancient(): updating the LDSC repo must not force re-running every trait
        script=ancient(LDSC.get("ldsc_script", "../ldsc/ldsc.py"))
    output:
        # LDSC writes {out}.log itself (buffered); stdout carries the same messages live
        "results/ldsc_h2/{h2_trait}.log"
    log:
        "results/ldsc_h2/stdout/{h2_trait}.log"
    params:
        out_prefix="results/ldsc_h2/{h2_trait}",
        ld_ref=LDSC.get("ld_ref_chr", "data/reference/LDscores/"),
        w_ld=LDSC.get("w_ld_chr", "data/reference/LDscores/"),
        # --samp-prev/--pop-prev only for the liability scale; --intercept-h2 only if fixed
        extra=lambda wildcards: " ".join(h2_flags(wildcards.h2_trait, H2))
    shell:
        """
        {PYTHON} -u {input.script} \
            --h2 {input.sumstats} \
            --ref-ld-chr {params.ld_ref} \
            --w-ld-chr {params.w_ld} \
            {params.extra} \
            --out {params.out_prefix} \
            > {log} 2>&1
        """

def h2_summary_inputs(wildcards):
    if not H2_TRAITS:
        raise ValueError("No heritability traits: set config h2.traits (>= 1 trait).")
    return [f"results/ldsc_h2/{trait}.log" for trait in H2_TRAITS]

rule h2_summary:
    input:
        h2_summary_inputs
    output:
        "results/ldsc_h2/h2_summary.tsv"
    shell:
        """
        {PYTHON} scripts/ldsc_h2_summary.py --logs {input} --output {output}
        """

# ------------------------------------------------------------------------------
# Step A1: DNA ASM Tag SNP Overlap with GWAS + BH FDR, then LD-Proxy Expansion (Per-Trait)
# ------------------------------------------------------------------------------
rule overlap_dna_asm:
    input:
        gwas=lambda wildcards: GWAS_FILES[wildcards.trait],
        mapping_table=config["dna"]["mapping_table"]
    output:
        all="results/dna/{trait}_overlap_all.tsv",
        sig="results/dna/{trait}_overlap_sig.tsv"
    params:
        alpha=FDR_ALPHA,
        tag_fallback=config["dna"].get("tag_fallback", "none")
    shell:
        """
        {PYTHON} scripts/overlap_dna_asm.py \
            --mapping-table {input.mapping_table} \
            --gwas {input.gwas} \
            --trait {wildcards.trait} \
            --alpha {params.alpha} \
            --tag-fallback {params.tag_fallback} \
            --output-all {output.all} \
            --output-sig {output.sig}
        """

# ------------------------------------------------------------------------------
# Step A2: Map BH-Significant CpGs to GENCODE v19 Gene TSS (< 5 kb, Per-Trait)
# ------------------------------------------------------------------------------
rule map_cpg_to_tss_genes:
    input:
        overlap="results/dna/{trait}_overlap_sig.tsv",
        cpg_coordinates=config["dna"]["cpg_coordinates"],
        gene_tss=config["dna"]["gene_tss"]
    output:
        "results/dna/{trait}_cpg_tss.tsv"
    params:
        cutoff=config["dna"].get("tss_distance_cutoff", 5000)
    shell:
        """
        {PYTHON} scripts/map_cpg_tss.py \
            --overlap {input.overlap} \
            --cpg-coords {input.cpg_coordinates} \
            --gene-tss {input.gene_tss} \
            --cutoff {params.cutoff} \
            --output {output}
        """

# ------------------------------------------------------------------------------
# Step A3: Query GTEx v10 Brain eQTLs (13 Tissues) for BH-Significant SNPs (Per-Trait)
# ------------------------------------------------------------------------------
rule lookup_gtex_brain_eqtls:
    input:
        overlap="results/dna/{trait}_overlap_sig.tsv"
    output:
        "results/dna/{trait}_gtex_brain_eqtls.tsv"
    params:
        tissues=" ".join(config["gtex"]["brain_tissues"]),
        p_threshold=config["gtex"].get("p_threshold", 0.05),
        cache_file=config["gtex"].get("cache_file", "results/.cache/gtex_eqtl_cache.json")
    shell:
        """
        {PYTHON} scripts/lookup_gtex_brain_eqtls.py \
            --overlap {input.overlap} \
            --output {output} \
            --tissues {params.tissues} \
            --p-threshold {params.p_threshold} \
            --cache-file {params.cache_file}
        """

# ------------------------------------------------------------------------------
# Step A4: CpG-Proximal Gene vs. GTEx Brain eQTL Target Gene Concordance (Per-Trait)
# ------------------------------------------------------------------------------
rule compare_tss_and_gtex_genes:
    input:
        tss="results/dna/{trait}_cpg_tss.tsv",
        gtex="results/dna/{trait}_gtex_brain_eqtls.tsv"
    output:
        "results/dna/{trait}_concordance.tsv"
    shell:
        """
        {PYTHON} scripts/compare_tss_gtex.py \
            --tss {input.tss} \
            --gtex {input.gtex} \
            --trait {wildcards.trait} \
            --output {output}
        """

# ------------------------------------------------------------------------------
# Step A5: Cross-Trait DNA Methylation Pleiotropy Summary
# ------------------------------------------------------------------------------
rule pleiotropy_dna:
    input:
        expand("results/dna/{trait}_concordance.tsv", trait=TRAITS)
    output:
        "results/dna/pleiotropy_dna.tsv"
    shell:
        """
        {PYTHON} scripts/pleiotropy.py --mode dna --inputs {input} --output {output}
        """

# ------------------------------------------------------------------------------
# Step B1: Overlap RNA m6A SNP List with GWAS + BH FDR (Per-Trait)
# ------------------------------------------------------------------------------
rule overlap_rna_m6a:
    input:
        gwas=lambda wildcards: GWAS_FILES[wildcards.trait],
        snp_list=config["rna"]["snp_list"]
    output:
        "results/rna/{trait}_overlap_bh.tsv"
    params:
        alpha=FDR_ALPHA
    shell:
        """
        {PYTHON} scripts/overlap_rna_m6a.py \
            --snp-list {input.snp_list} \
            --gwas {input.gwas} \
            --trait {wildcards.trait} \
            --alpha {params.alpha} \
            --output {output}
        """

# ------------------------------------------------------------------------------
# Step B2: Genes Overlapping the BH-Significant RNA m6A SNPs (hg19/GRCh37, All Traits)
# Source (config rna_gene_source): "ucsc" = UCSC REST API, GENCODE V19 (default);
# "ensembl" = Ensembl VEP REST, GRCh37. Each source has its own server and cache
# (config rna.{source}_server / rna.{source}_cache_file).
# ------------------------------------------------------------------------------
RNA_GENE_SOURCE = str(config.get("rna_gene_source") or "ucsc").strip().lower()
if RNA_GENE_SOURCE not in ("ucsc", "ensembl"):
    raise ValueError(f"config rna_gene_source must be 'ucsc' or 'ensembl', not {config.get('rna_gene_source')!r}")
RNA_GENE_DEFAULTS = {
    "ucsc": ("https://api.genome.ucsc.edu", "results/.cache/ucsc_hg19_genes_cache.json"),
    "ensembl": ("https://grch37.rest.ensembl.org", "results/.cache/ensembl_grch37_genes_cache.json"),
}

rule annotate_rna_genes:
    input:
        expand("results/rna/{trait}_overlap_bh.tsv", trait=TRAITS)
    output:
        "results/rna/rna_sig_genes.tsv"
    params:
        source=RNA_GENE_SOURCE,
        server=config["rna"].get(f"{RNA_GENE_SOURCE}_server", RNA_GENE_DEFAULTS[RNA_GENE_SOURCE][0]),
        cache_file=config["rna"].get(f"{RNA_GENE_SOURCE}_cache_file", RNA_GENE_DEFAULTS[RNA_GENE_SOURCE][1])
    shell:
        """
        {PYTHON} scripts/annotate_rna_genes.py \
            --source {params.source} \
            --inputs {input} \
            --output {output} \
            --server {params.server} \
            --cache-file {params.cache_file}
        """

# ------------------------------------------------------------------------------
# Step B3: Cross-Trait RNA Methylation Pleiotropy Summary
# ------------------------------------------------------------------------------
rule pleiotropy_rna:
    input:
        overlaps=expand("results/rna/{trait}_overlap_bh.tsv", trait=TRAITS),
        genes="results/rna/rna_sig_genes.tsv"
    output:
        "results/rna/pleiotropy_rna.tsv"
    shell:
        """
        {PYTHON} scripts/pleiotropy.py --mode rna --inputs {input.overlaps} --genes {input.genes} --output {output}
        """
