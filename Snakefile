# ==============================================================================
# LDSC & GC Automated Pipeline (VIBRA HPC)
# Compatible with Snakemake 9 & Slurm Executor Plugin
# ==============================================================================

import os
import sys
import itertools

# Ensure logs directory exists for cluster execution
os.makedirs("logs", exist_ok=True)

# Load configuration file (merges local overrides if config_local.yaml exists)
configfile: "config.yaml"
if os.path.exists("config_local.yaml"):
    configfile: "config_local.yaml"

# ----- TRAIT & PAIR DETECTION ----- #
# 1. Detect traits from config['traits'] (with automatic fallback to raw_data folder)
DISCOVERED_TRAITS, = glob_wildcards("data/raw_data/{trait}_gwas.txt")
TRAITS = config.get("traits", DISCOVERED_TRAITS) or DISCOVERED_TRAITS

def get_gwas_input(wildcards):
    """Resolve raw GWAS summary statistics path for {trait} from config or default pattern."""
    gwas_map = config.get("gwas_sumstats", {})
    return gwas_map.get(wildcards.trait, f"data/raw_data/{wildcards.trait}_gwas.txt")

# 2. Automatically generate every unique pair (e.g., A_B, A_C, B_C)
PAIRS = [f"{t1}_{t2}" for t1, t2 in itertools.combinations(TRAITS, 2)]
# -----                        ----- #

# ----- METHYLATION CONFIGURATION ----- #
METH_TYPES = config.get("meth_types", ["DNA_ASM", "RNA_m6A"])

# ==========================================
# PIPELINE SHORTCUTS (DUMMY RULES)
# ==========================================

# 1. Run everything (Munge -> Correlations -> Plots + Per-Trait & Cross-Trait DNA ASM & RNA m6A Methylation)
# Snakemake 9 Profile Command: snakemake --profile profiles/slurm all
rule all:
    input:
        "results/figures/rg_heatmap.pdf",
        "results/figures/rg_network.pdf",
        expand(
            "results/methylation/{meth}_{trait}_overlap.txt",
            meth=METH_TYPES,
            trait=TRAITS
        ),
        expand(
            "results/methylation/{meth}_{trait}_overlap_fullheaders.txt",
            meth=METH_TYPES,
            trait=TRAITS
        ),
        expand("results/methylation/dna_asm_{trait}_cpg_gene_tss_matches.tsv", trait=TRAITS),
        expand("results/methylation/dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv", trait=TRAITS),
        expand("results/methylation/dna_asm_{trait}_concordant_genes_summary.tsv", trait=TRAITS),
        "results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv",
        "results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps.tsv",
        "results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps_annotated.tsv",
        "results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps_summary.tsv",
        "results/figures/cross_trait_rna_m6a_recurrent_snps_annotated_effect_heatmap.pdf"

# 2. Run correlations and stop (No plotting)
# Snakemake 9 Profile Command: snakemake --profile profiles/slurm correlations
rule correlations:
    input:
        expand("results/ldsc_rg/{pair}.log", pair=PAIRS)

# 3. Run single-trait SNP heritability (h2) only
# Snakemake 9 Profile Command: snakemake --profile profiles/slurm heritability
rule heritability:
    input:
        expand("results/ldsc_h2/{trait}_h2.log", trait=TRAITS)

# 4. Only prep the raw GWAS data and stop
# Snakemake 9 Profile Command: snakemake --profile profiles/slurm munge_only
rule munge_only:
    input:
        expand("results/munged/{trait}.sumstats.gz", trait=TRAITS)

# 5. Run methylation overlap & functional convergence analysis step only (DNA ASM + RNA m6A)
# Snakemake 9 Profile Command: snakemake --profile profiles/slurm --jobs 10 methylation
rule methylation:
    input:
        expand(
            "results/methylation/{meth}_{trait}_overlap.txt",
            meth=METH_TYPES,
            trait=TRAITS
        ),
        expand(
            "results/methylation/{meth}_{trait}_overlap_fullheaders.txt",
            meth=METH_TYPES,
            trait=TRAITS
        ),
        expand(
            "results/methylation/dna_asm_{trait}_cpg_gene_tss_matches.tsv",
            trait=TRAITS
        ),
        expand(
            "results/methylation/dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv",
            trait=TRAITS
        ),
        expand(
            "results/methylation/dna_asm_{trait}_concordant_genes_summary.tsv",
            trait=TRAITS
        ),
        "results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv",
        "results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps.tsv",
        "results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps_annotated.tsv",
        "results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps_summary.tsv",
        "results/figures/cross_trait_rna_m6a_recurrent_snps_annotated_effect_heatmap.pdf"

def get_col_flags(wildcards, input):
    """Reads the first line of the GWAS file and matches headers to LDSC flags."""
    # Open the file and read just the first line, splitting it by tabs or spaces
    with open(input.sumstats, 'r') as f:
        header = f.readline().strip().split()
        
    # Your master dictionary of all possible column names
    aliases = {
        "--snp": ["SNP", "snp", "rsID", "rsid"],
        "--a1": ["A1", "a1", "Allele1", "effect_allele", "alt", "Aallele"],
        "--a2": ["A2", "a2", "Allele2", "other_allele", "ref", "Ballele"],
        "--p": ["P", "p", "pvalue", "p_value", "P_value", "P.value"],
        "--N-col": ["N", "n", "N_eff", "Neff", "N_EFF", "N_TOTAL", "Weight"]
    }

    # Dictionary for sumstats that lack an N column, add more if needed
    constant_n = {
        "SUD2023": "1025550"
    }
    
    flags = []
    for ldsc_arg, possible_names in aliases.items():
        # Look for the first name in our list that actually exists in the file's header
        match = next((col for col in possible_names if col in header), None)
        
        if match:
            flags.extend([ldsc_arg, match]) # e.g., adds "--snp" and "rsid"
        else:
            # Fallback logic for missing N
            if ldsc_arg == "--N-col" and wildcards.trait in constant_n:
                flags.extend(["--N", constant_n[wildcards.trait]]) # Use --N instead of --N-col
            else:
                # If no match is found, crash with a helpful error message
                raise ValueError(
                    f"\nERROR: Could not find a match for '{ldsc_arg}' in {input.sumstats}.\n"
                    f"File header contains: {header}\n"
                    f"Please add the correct column name to the aliases dictionary in the Snakefile.\n"
                    f"Or if 'N' is missing, add the trait to the constant_n dictionary with its sample size.\n"
                )
            
    # Joins the list into a single string: "--snp rsid --a1 Allele1 --a2 Allele2 ..."
    return " ".join(flags)

# Helper to format conda environment activation if enabled
def get_env_setup(env_key):
    if config.get("use_conda_activation", True):
        conda_sh = config.get("conda_sh", "/aplic/anaconda/2024.02/etc/profile.d/conda.sh")
        target_env = config.get(env_key, "ldsc")
        return f"source {conda_sh}\n        conda activate {target_env}\n\n        "
    return ""

# Helper to resolve the Python interpreter (uses active Snakemake Python when running locally)
def get_python_bin():
    py_cfg = config.get("python_bin", "python")
    if not config.get("use_conda_activation", True) and py_cfg == "python":
        return f'"{sys.executable.replace(chr(92), "/")}"'
    return py_cfg

# ------------------------------------------------------------------------------
# Step 1: Munge Summary Statistics
# ------------------------------------------------------------------------------
rule munge_sumstats:
    input:
        sumstats=get_gwas_input,
        snplist=config.get("snplist", "data/reference/w_hm3.snplist")
    output:
        "results/munged/{trait}.sumstats.gz"
    priority: 50
    params:
        # LDSC expects the output prefix without the '.sumstats.gz' extension
        out_prefix="results/munged/{trait}",
        col_flags=get_col_flags,
        env_setup=lambda wildcards: get_env_setup("ldsc_env"),
        python_bin=get_python_bin(),
        munge_script=config.get("munge_script", "../ldsc/munge_sumstats.py"),
        # Legacy parameter compatibility
        time="06:00:00",
        partition=config.get("slurm_partition", "irbio01")
    resources:
        mem_mb=16000,
        cpus=1,
        runtime=360,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/munge_sumstats_%j.out --error=logs/munge_sumstats_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} {params.munge_script} \
            --sumstats {input.sumstats} \
            --merge-alleles {input.snplist} \
            --out {params.out_prefix} \
            {params.col_flags} \
            --chunksize 500000
        """

# ------------------------------------------------------------------------------
# Step 1b: Single-Trait SNP Heritability (h2)
# ------------------------------------------------------------------------------
rule ldsc_h2:
    input:
        sumstats="results/munged/{trait}.sumstats.gz",
        all_munged=expand("results/munged/{t}.sumstats.gz", t=TRAITS)
    output:
        "results/ldsc_h2/{trait}_h2.log"
    priority: 20
    params:
        out_prefix="results/ldsc_h2/{trait}_h2",
        ld_ref=config.get("ld_ref", "data/reference/LDscores"),
        env_setup=lambda wildcards: get_env_setup("ldsc_env"),
        python_bin=get_python_bin(),
        ldsc_script=config.get("ldsc_script", "../ldsc/ldsc.py"),
        time="01:00:00",
        partition=config.get("slurm_partition", "irbio01")
    resources:
        mem_mb=16000,
        cpus=1,
        runtime=60,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/ldsc_h2_%j.out --error=logs/ldsc_h2_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} {params.ldsc_script} \
            --h2 {input.sumstats} \
            --ref-ld-chr {params.ld_ref}/ \
            --w-ld-chr {params.ld_ref}/ \
            --out {params.out_prefix}
        """

# ------------------------------------------------------------------------------
# Step 2: Perform LD Score Regression to Estimate Genetic Correlations
# ------------------------------------------------------------------------------
rule ldsc_rg:
    input:
        # Using wildcards here allows you to pair any two traits dynamically
        # all_munged ensures all traits finish Step 1 (munge) before Step 2 starts
        trait1="results/munged/{trait1}.sumstats.gz",
        trait2="results/munged/{trait2}.sumstats.gz",
        all_munged=expand("results/munged/{t}.sumstats.gz", t=TRAITS)
    output:
        "results/ldsc_rg/{trait1}_{trait2}.log"
    priority: 20
    params:
        out_prefix="results/ldsc_rg/{trait1}_{trait2}",
        ld_ref=config.get("ld_ref", "data/reference/LDscores"),
        env_setup=lambda wildcards: get_env_setup("ldsc_env"),
        python_bin=get_python_bin(),
        ldsc_script=config.get("ldsc_script", "../ldsc/ldsc.py"),
        # Legacy parameter compatibility
        time="02:00:00",
        partition=config.get("slurm_partition", "irbio01")
    resources:
        mem_mb=16000,
        cpus=1,
        runtime=120,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/ldsc_rg_%j.out --error=logs/ldsc_rg_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} {params.ldsc_script} \
            --rg {input.trait1},{input.trait2} \
            --ref-ld-chr {params.ld_ref}/ \
            --w-ld-chr {params.ld_ref}/ \
            --out {params.out_prefix}
        """

# ------------------------------------------------------------------------------
# Step 3: Plot Heatmap for Genetic Correlations
# ------------------------------------------------------------------------------
rule plot_heatmap:
    input:
        # Require every single correlation log to be finished before plotting
        logs=expand("results/ldsc_rg/{pair}.log", pair=PAIRS)
    output:
        "results/figures/rg_heatmap.pdf"
    params:
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        # Legacy parameter compatibility
        time="00:30:00",
        partition=config.get("slurm_partition", "irbio01")
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=30,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/plot_heatmap_%j.out --error=logs/plot_heatmap_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/plot_heatmap.py {input.logs} {output}
        """

# ------------------------------------------------------------------------------
# Step 4: Plot Network for Genetic Correlations
# ------------------------------------------------------------------------------
rule plot_network:
    input:
        logs=expand("results/ldsc_rg/{pair}.log", pair=PAIRS)
    output:
        "results/figures/rg_network.pdf"
    params:
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        # Legacy parameter compatibility
        time="00:30:00",
        partition=config.get("slurm_partition", "irbio01")
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=30,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/plot_network_%j.out --error=logs/plot_network_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/plot_genetic_correlations.py {input.logs} {output}
        """

# ------------------------------------------------------------------------------
# Step 5: Methylation analysis, overlap DNA & RNA SNP lists with GWAS
# ------------------------------------------------------------------------------
rule methylation_overlap:
    input:
        gwas=get_gwas_input,
        snp_list="data/methylation/{meth}_SNPlist.txt"
    output:
        overlap_ids="results/methylation/{meth}_{trait}_overlap.txt",
        overlap_full="results/methylation/{meth}_{trait}_overlap_fullheaders.txt"
    params:
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        # Legacy parameter compatibility
        partition=config.get("slurm_partition", "irbio01"),
        time="00:15:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=15,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/methylation_overlap_%j.out --error=logs/methylation_overlap_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/overlap_methylation.py {input.snp_list} {input.gwas} {output.overlap_full} {output.overlap_ids}
        """

# ------------------------------------------------------------------------------
# Step 5a: Intersect GWAS Summary Statistics with DNA ASM Tag SNPs & LD Proxies (Per-Trait)
# ------------------------------------------------------------------------------
rule overlap_dna_asm:
    input:
        gwas=get_gwas_input,
        mapping_table=config["methylation"]["dna_mapping_table"]
    output:
        "results/methylation/dna_asm_{trait}_gwas_overlap.tsv"
    params:
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        partition=config.get("slurm_partition", "irbio01"),
        time="00:15:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=15,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/overlap_dna_asm_%j.out --error=logs/overlap_dna_asm_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/overlap_dna_asm.py \
            --mapping-table {input.mapping_table} \
            --gwas {input.gwas} \
            --output {output}
        """

# ------------------------------------------------------------------------------
# Step 5b: Map Overlapping CpGs to hg19 Coordinates & GENCODE v19 Gene TSS (< 5 kb, Per-Trait)
# ------------------------------------------------------------------------------
rule map_cpg_to_tss_genes:
    input:
        overlap="results/methylation/dna_asm_{trait}_gwas_overlap.tsv",
        cpg_coordinates=config["methylation"]["cpg_coordinates"],
        gene_tss=config["methylation"]["gene_tss"]
    output:
        "results/methylation/dna_asm_{trait}_cpg_gene_tss_matches.tsv"
    params:
        tss_distance_cutoff=config["methylation"].get("tss_distance_cutoff", 5000),
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        partition=config.get("slurm_partition", "irbio01"),
        time="00:15:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=15,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/map_cpg_to_tss_genes_%j.out --error=logs/map_cpg_to_tss_genes_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/map_cpg_tss.py \
            --overlap {input.overlap} \
            --cpg-coords {input.cpg_coordinates} \
            --gene-tss {input.gene_tss} \
            --cutoff {params.tss_distance_cutoff} \
            --output {output}
        """

# ------------------------------------------------------------------------------
# Step 5c: Query GTEx v10 Across 13 Brain Tissues for DNA ASM Overlapping Variants
# ------------------------------------------------------------------------------
rule lookup_gtex_brain_eqtls:
    input:
        overlap="results/methylation/dna_asm_{trait}_gwas_overlap.tsv"
    output:
        "results/methylation/dna_asm_{trait}_gtex_brain_eqtls.tsv"
    params:
        tissues=" ".join(config.get("gtex", {}).get("brain_tissues", [])),
        p_threshold=config.get("gtex", {}).get("p_threshold", 0.05),
        gwas_p_threshold=config.get("gtex", {}).get("gwas_p_threshold", 0.001),
        cache_file=config.get("gtex", {}).get("cache_file", "results/methylation/.gtex_eqtl_cache.json"),
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        partition=config.get("slurm_partition", "irbio01"),
        time="00:30:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=30,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/lookup_gtex_brain_eqtls_%j.out --error=logs/lookup_gtex_brain_eqtls_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/lookup_gtex_brain_eqtls.py \
            --overlap {input.overlap} \
            --output {output} \
            --tissues {params.tissues} \
            --p-threshold {params.p_threshold} \
            --gwas-p-threshold {params.gwas_p_threshold} \
            --cache-file {params.cache_file}
        """

# ------------------------------------------------------------------------------
# Step 5d: Compare Nearest CpG TSS Genes (< 5 kb) with GTEx Brain eQTL Target Genes
# ------------------------------------------------------------------------------
rule compare_tss_and_gtex_genes:
    input:
        tss="results/methylation/dna_asm_{trait}_cpg_gene_tss_matches.tsv",
        gtex="results/methylation/dna_asm_{trait}_gtex_brain_eqtls.tsv"
    output:
        full="results/methylation/dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv",
        concordant="results/methylation/dna_asm_{trait}_concordant_genes_summary.tsv"
    params:
        trait="{trait}",
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        partition=config.get("slurm_partition", "irbio01"),
        time="00:15:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=15,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/compare_tss_and_gtex_genes_%j.out --error=logs/compare_tss_and_gtex_genes_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/compare_tss_gtex.py \
            --tss {input.tss} \
            --gtex {input.gtex} \
            --trait {params.trait} \
            --output-full {output.full} \
            --output-concordant {output.concordant}
        """

# ------------------------------------------------------------------------------
# Step 5e: Cross-Trait DNA ASM FDR Recurrence & Functional Concordance Summary
# ------------------------------------------------------------------------------
rule cross_trait_asm_summary:
    input:
        expand("results/methylation/dna_asm_{trait}_tss_vs_gtex_full_comparison.tsv", trait=TRAITS)
    output:
        tsv="results/methylation/comparisons/cross_trait_asm_recurrent_snps.tsv"
    params:
        method=config.get("cross_trait", {}).get("fdr_method", "storey"),
        k=config.get("cross_trait", {}).get("min_recurrence", 2),
        threshold=config.get("cross_trait", {}).get("fdr_threshold", 0.05),
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        partition=config.get("slurm_partition", "irbio01"),
        time="00:15:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=15,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/cross_trait_asm_summary_%j.out --error=logs/cross_trait_asm_summary_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/interactive_cross_trait_asm.py \
            --input-dir results/methylation \
            --method {params.method} \
            --threshold {params.threshold} \
            --min-traits {params.k} \
            --output {output.tsv}
        """

# ------------------------------------------------------------------------------
# Step 5f: Cross-Trait RNA m6A FDR Recurrence, Ensembl/GTEx Annotation & Heatmap
# ------------------------------------------------------------------------------
rule cross_trait_rna_m6a_summary:
    input:
        expand("results/methylation/RNA_m6A_{trait}_overlap_fullheaders.txt", trait=TRAITS)
    output:
        sig_tsv="results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps.tsv",
        annotated="results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps_annotated.tsv",
        summary="results/methylation/comparisons/cross_trait_rna_m6a_recurrent_snps_summary.tsv",
        heatmap_pdf="results/figures/cross_trait_rna_m6a_recurrent_snps_annotated_effect_heatmap.pdf",
        heatmap_png="results/figures/cross_trait_rna_m6a_recurrent_snps_annotated_effect_heatmap.png"
    params:
        method=config.get("cross_trait", {}).get("fdr_method", "storey"),
        k=config.get("cross_trait", {}).get("min_recurrence", 2),
        env_setup=lambda wildcards: get_env_setup("snakemake_env"),
        python_bin=get_python_bin(),
        partition=config.get("slurm_partition", "irbio01"),
        time="00:30:00"
    resources:
        mem_mb=4000,
        cpus=1,
        runtime=30,
        slurm_partition=config.get("slurm_partition", "irbio01"),
        slurm_extra="'--output=logs/cross_trait_rna_m6a_summary_%j.out --error=logs/cross_trait_rna_m6a_summary_%j.err'"
    shell:
        """
        {params.env_setup}{params.python_bin} scripts/run_interactive_pipeline.py \
            --meth RNA_m6A \
            --traits all \
            --min-k {params.k} \
            --fdr-method {params.method} \
            --metric effect \
            --threads 8 \
            --output {output.sig_tsv}
        """