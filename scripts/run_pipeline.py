#!/usr/bin/env python3
"""
Pipeline Launcher for app.py (Local Background Runs) and scripts/run_pipeline.sbatch (HPC).

Runs Snakemake from the repository root with the interpreter of the active environment,
streams its output to stdout, and writes the exit code to --status-file when finished.

Usage:
    python scripts/run_pipeline.py --target all --traits SUD2023,Anxiety --engine local --cores 4
    python scripts/run_pipeline.py --target dna --engine slurm
    python scripts/run_pipeline.py --target munge --traits SUD2023 --forcerun munge_sumstats
    python scripts/run_pipeline.py --target munge --extra-config 'ldsc={"n_override": {"SUD2023": 1025550}}'
    python scripts/run_pipeline.py --target rg --traits A,B --extra-config 'rg={"traits": ["A", "B"], "target": ""}'
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from common_utils import PROJECT_ROOT


def build_command(
    target: str, traits: str, engine: str, cores: int, extra_config: list = None, forcerun: str = None
) -> list:
    cmd = [sys.executable, "-m", "snakemake", target, "--rerun-incomplete"]
    if engine == "slurm":
        cmd += ["--profile", "profiles/slurm"]
    else:
        cmd += ["--cores", str(cores)]
    config_items = ([f"traits={traits}"] if traits else []) + list(extra_config or [])
    if config_items:
        cmd += ["--config", *config_items]
    if forcerun:
        cmd += ["--forcerun", forcerun]
    return cmd


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the methylation Snakemake pipeline.")
    parser.add_argument("--target", choices=["all", "dna", "rna", "munge", "rg"], default="all", help="Snakemake target rule")
    parser.add_argument("--traits", default="", help="Comma-separated trait names (default: config)")
    parser.add_argument("--engine", choices=["local", "slurm"], default="local", help="Execution engine")
    parser.add_argument("--cores", type=int, default=2, help="Local cores (ignored for slurm)")
    parser.add_argument("--status-file", default=None, help="File receiving the exit code when done")
    parser.add_argument(
        "--extra-config", action="append", default=[], help="Additional Snakemake KEY=VALUE config (repeatable)"
    )
    parser.add_argument("--forcerun", default=None, help="Rule whose jobs are re-run even if up to date")
    args = parser.parse_args()

    cmd = build_command(
        args.target, args.traits, args.engine, args.cores, args.extra_config, args.forcerun
    )
    print("[run_pipeline] " + " ".join(cmd), flush=True)

    env = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")
    exit_code = subprocess.call(cmd, cwd=PROJECT_ROOT, env=env)
    print(f"[run_pipeline] Snakemake exited with code {exit_code}", flush=True)

    if args.status_file:
        Path(args.status_file).write_text(str(exit_code), encoding="utf-8")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
