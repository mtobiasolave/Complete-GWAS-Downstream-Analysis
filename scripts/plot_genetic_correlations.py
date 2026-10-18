#!/usr/bin/env python3
"""Generate circular LDSC genetic correlation network PDF (FDR < 0.05)."""

import os
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from common_utils import compute_fdr, parse_ldsc_rg_logs

log_files = sys.argv[1:-1]
output_file = sys.argv[-1]

parsed = parse_ldsc_rg_logs(log_files)

if not parsed.empty:
    df = pd.DataFrame({
        "Trait1": parsed["trait1"],
        "Trait2": parsed["trait2"],
        "rg": parsed["rg"],
        "p_value": parsed["p"].fillna(1.0),
    })

    df["p_adjusted"] = compute_fdr(df["p_value"], method="bh")
    df["FDR_significant"] = df["p_adjusted"] <= 0.05
    sig_df = df[df["FDR_significant"]].copy()

    G = nx.Graph()
    traits = sorted(pd.concat([df["Trait1"], df["Trait2"]]).unique())
    G.add_nodes_from(traits)

    for row in sig_df.itertuples(index=False):
        G.add_edge(row.Trait1, row.Trait2, weight=row.rg, p=row.p_value)

    fig, ax = plt.subplots(figsize=(16, 16))
    pos = nx.circular_layout(G)

    edges = G.edges(data=True)
    weights = [data["weight"] for _, _, data in edges]

    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=600,
        node_color="lightblue",
        edgecolors="black",
        linewidths=1.5,
        ax=ax,
    )

    if edges:
        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=[(u, v) for u, v, _ in edges],
            width=[max(abs(w) * 4, 1.0) for w in weights],
            edge_color=weights,
            edge_cmap=plt.cm.coolwarm,
            edge_vmin=-1.0,
            edge_vmax=1.0,
            alpha=0.75,
            ax=ax,
        )

    label_ratio = 1.10
    for node, (x, y) in pos.items():
        angle = np.arctan2(y, x)
        lx = label_ratio * np.cos(angle)
        ly = label_ratio * np.sin(angle)

        if abs(x) < 0.08:
            ha = "center"
        elif x > 0:
            ha = "left"
        else:
            ha = "right"

        if y > 0.85:
            va = "bottom"
        elif y < -0.85:
            va = "top"
        else:
            va = "center"

        ax.text(
            lx,
            ly,
            s=node,
            fontsize=12,
            fontweight="bold",
            ha=ha,
            va=va,
            rotation=0,
        )

    ax.set_title(
        "Genetic Correlation Network (FDR < 0.05)",
        fontsize=18,
        fontweight="bold",
        pad=40,
    )
    ax.axis("off")
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)

    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=-1.0, vmax=1.0)
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.05)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(
        "Genetic Correlation ($r_g$)", fontsize=14, fontweight="bold", labelpad=12
    )

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close()