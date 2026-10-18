#!/usr/bin/env python3
"""Generate pairwise LDSC genetic correlation lower-triangle heatmap PDF."""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from common_utils import parse_ldsc_rg_logs

log_files = sys.argv[1:-1]
output_file = sys.argv[-1]

parsed = parse_ldsc_rg_logs(log_files)
if parsed.empty:
    sys.stderr.write("ERROR: No genetic correlation results parsed from log files.\n")
    sys.exit(1)

p_vals = parsed["p"].fillna(1.0)
forward_df = pd.DataFrame({
    "Trait1": parsed["trait1"],
    "Trait2": parsed["trait2"],
    "rg": parsed["rg"],
    "p_value": p_vals,
})
reverse_df = pd.DataFrame({
    "Trait1": parsed["trait2"],
    "Trait2": parsed["trait1"],
    "rg": parsed["rg"],
    "p_value": p_vals,
})

traits = pd.concat([forward_df["Trait1"], forward_df["Trait2"]]).unique()
diag_df = pd.DataFrame({
    "Trait1": traits,
    "Trait2": traits,
    "rg": 1.0,
    "p_value": 0.0,
})

df = pd.concat([forward_df, reverse_df, diag_df], ignore_index=True)

# Vectorized significance asterisk assignment
df["asterisks"] = np.select(
    [df["p_value"] < 0.001, df["p_value"] < 0.01, df["p_value"] < 0.05],
    ["***", "**", "*"],
    default="",
)

rg_matrix = df.pivot(index="Trait2", columns="Trait1", values="rg").reindex(
    index=traits, columns=traits
)
asterisk_matrix = df.pivot(index="Trait2", columns="Trait1", values="asterisks").reindex(
    index=traits, columns=traits
)

# Hide upper triangle and diagonal, and trim empty first row / last column
mask = np.triu(np.ones_like(rg_matrix, dtype=bool))
rg_matrix = rg_matrix.iloc[1:, :-1]
asterisk_matrix = asterisk_matrix.iloc[1:, :-1]
mask = mask[1:, :-1]

plt.figure(figsize=(18, 14))
ax = sns.heatmap(
    rg_matrix,
    mask=mask,
    annot=asterisk_matrix,
    fmt="",
    cmap="coolwarm",
    center=0,
    vmin=-1.0,
    vmax=1.0,
    square=True,
    annot_kws={"size": 20},
    cbar_kws={"label": "Genetic Correlation ($r_g$)"},
)

plt.xlabel(None)
plt.ylabel(None)

ax.set_xticklabels(
    ax.get_xticklabels(),
    fontsize=20,
    fontweight="bold",
    rotation=45,
    ha="right",
)

ax.set_yticklabels(
    ax.get_yticklabels(),
    fontsize=20,
    fontweight="bold",
    rotation=0,
)

cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)
cbar.set_label(
    "Genetic Correlation ($r_g$)", fontsize=22, fontweight="bold", labelpad=15
)

plt.tight_layout()
os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.close()