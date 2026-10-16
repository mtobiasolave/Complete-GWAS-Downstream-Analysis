import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from statsmodels.stats.multitest import multipletests

# 1. Load files from Snakemake arguments
log_files = sys.argv[1:-1]
output_file = sys.argv[-1]

data = []

# 2. Parse LDSC .log files
for log in log_files:
    with open(log, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        if line.startswith("Summary of Genetic Correlation Results"):
            values = lines[i+2].split()
            
            t1 = os.path.basename(values[0]).replace(".sumstats.gz", "")
            t2 = os.path.basename(values[1]).replace(".sumstats.gz", "")
            rg = float(values[2]) if values[2] not in ("NA", "nan", "NaN") else np.nan
            pval = float(values[5]) if values[5] not in ("NA", "nan", "NaN") else 1.0
            
            data.append({"Trait1": t1, "Trait2": t2, "rg": rg, "p_value": pval})
            break

df = pd.DataFrame(data)

if not df.empty:
    # 3. Apply 5% FDR correction (Benjamini-Hochberg)
    reject, pvals_corrected, _, _ = multipletests(df['p_value'], alpha=0.05, method='fdr_bh')
    df['FDR_significant'] = reject
    df['p_adjusted'] = pvals_corrected

    # Filter for significant correlations at 5% FDR
    sig_df = df[df['FDR_significant']].copy()

    # 4. Build NetworkX Graph
    G = nx.Graph()
    
    # Sort traits alphabetically so the circle has a clean, deterministic order
    traits = sorted(pd.concat([df['Trait1'], df['Trait2']]).unique())
    G.add_nodes_from(traits)

    for _, row in sig_df.iterrows():
        G.add_edge(row['Trait1'], row['Trait2'], weight=row['rg'], p=row['p_value'])

    # 5. Perfect Circular Layout
    fig, ax = plt.subplots(figsize=(16, 16))
    pos = nx.circular_layout(G)

    edges = G.edges(data=True)
    weights = [data['weight'] for u, v, data in edges]
    
    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos, 
        node_size=600, 
        node_color='lightblue', 
        edgecolors='black', 
        linewidths=1.5,
        ax=ax
    )

    # Draw edges with rg-based color gradient and width
    if edges:
        nx.draw_networkx_edges(
            G, pos, 
            edgelist=[(u, v) for u, v, d in edges],
            width=[max(abs(w) * 4, 1.0) for w in weights],
            edge_color=weights,
            edge_cmap=plt.cm.coolwarm,
            edge_vmin=-1.0,
            edge_vmax=1.0,
            alpha=0.75,
            ax=ax
        )

# 6. Completely horizontal radial label placement
    label_ratio = 1.10  # pushes labels 10% outside the unit circle radius
    for node, (x, y) in pos.items():
        angle = np.arctan2(y, x)
        lx = label_ratio * np.cos(angle)
        ly = label_ratio * np.sin(angle)
        
        # Horizontal alignment: left for right side, right for left side, center for exact top/bottom
        if abs(x) < 0.08:
            ha = 'center'
        elif x > 0:
            ha = 'left'
        else:
            ha = 'right'

        # Vertical alignment: bottom for top node, top for bottom node, center for middle nodes
        if y > 0.85:
            va = 'bottom'
        elif y < -0.85:
            va = 'top'
        else:
            va = 'center'
        
        ax.text(
            lx, ly, s=node, 
            fontsize=12, 
            fontweight='bold', 
            ha=ha, va=va,
            rotation=0  # All labels strictly horizontal
        )

    ax.set_title("Genetic Correlation Network (FDR < 0.05)", fontsize=18, fontweight="bold", pad=40)
    ax.axis('off')

    # Expand plot limits so outer text labels aren't clipped
    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.45, 1.45)
    
    # Colorbar for edge correlations
    sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=plt.Normalize(vmin=-1.0, vmax=1.0))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.05)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("Genetic Correlation ($r_g$)", fontsize=14, fontweight="bold", labelpad=12)

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, format="pdf", bbox_inches="tight")
    plt.close()