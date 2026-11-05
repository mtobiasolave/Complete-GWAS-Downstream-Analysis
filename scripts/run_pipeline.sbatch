#!/bin/bash
# ==============================================================================
# Snakemake Controller Job for VIBRA HPC
#
# Runs the Snakemake main process as a lightweight Slurm job, which in turn submits
# every pipeline step as its own Slurm job via profiles/slurm. The controller survives
# logouts and closed Streamlit sessions.
#
# Usage (from the repository root, with snakemake9_env activated):
#   mkdir -p logs/runs
#   sbatch scripts/run_pipeline.sbatch [target] [traits] [status_file] [run_pipeline.py options]
#   sbatch scripts/run_pipeline.sbatch all SUD2023,BipolarDisorder2024
#   sbatch scripts/run_pipeline.sbatch munge SUD2023 "" --forcerun munge_sumstats
#   sbatch scripts/run_pipeline.sbatch h2 SUD2023 "" --extra-config 'h2={"traits": ["SUD2023"]}'
#
# target: all | dna | rna | munge | rg | h2 (default: all); traits: comma-separated (default: config.yaml)
# ==============================================================================
#SBATCH --job-name=methyl_pipeline
#SBATCH --partition=irbio01
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/runs/controller_%j.log

TARGET="${1:-all}"
TRAITS="${2:-}"
STATUS_FILE="${3:-}"
shift $(( $# < 3 ? $# : 3 ))

python scripts/run_pipeline.py --engine slurm --target "$TARGET" --traits "$TRAITS" \
    ${STATUS_FILE:+--status-file "$STATUS_FILE"} "$@"
