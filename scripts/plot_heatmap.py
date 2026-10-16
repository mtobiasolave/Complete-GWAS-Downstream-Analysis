import sys
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np  # 1. ADDED NUMPY

log_files = sys.argv[1:-1]
output_file = sys.argv[-1]

data = []

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
            data.append({"Trait1": t2, "Trait2": t1, "rg": rg, "p_value": pval})
            break

df = pd.DataFrame(data)

traits = pd.concat([df["Trait1"], df["Trait2"]]).unique()
for t in traits:
    df = pd.concat([df, pd.DataFrame([{"Trait1": t, "Trait2": t, "rg": 1.0, "p_value": 0.0}])], ignore_index=True)

def get_asterisks(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return ''

df['asterisks'] = df['p_value'].apply(get_asterisks)

rg_matrix = df.pivot(index='Trait2', columns='Trait1', values='rg').reindex(index=traits, columns=traits)
asterisk_matrix = df.pivot(index='Trait2', columns='Trait1', values='asterisks').reindex(index=traits, columns=traits)

# 2. CREATE THE MASK (hides the upper triangle and the diagonal)
mask = np.triu(np.ones_like(rg_matrix, dtype=bool))

# 3. TRIM THE GHOST LABELS (Slice off the empty first row and last column)
rg_matrix = rg_matrix.iloc[1:, :-1]
asterisk_matrix = asterisk_matrix.iloc[1:, :-1]
mask = mask[1:, :-1]

plt.figure(figsize=(18, 14))
ax = sns.heatmap(
    rg_matrix, 
    mask=mask,              # 3. APPLY THE MASK HERE
    annot=asterisk_matrix,  
    fmt="",                 
    cmap="coolwarm",       
    center=0,              
    vmin=-1.0,             
    vmax=1.0,               
    square=True,            
    annot_kws={'size': 20}, 
    cbar_kws={'label': 'Genetic Correlation ($r_g$)'}
)

plt.xlabel(None)
plt.ylabel(None)

ax.set_xticklabels(
    ax.get_xticklabels(),
    fontsize=20,  
    fontweight="bold",
    rotation=45,  
    ha="right"
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

os.makedirs(os.path.dirname(output_file), exist_ok=True)
plt.savefig(output_file, format="pdf", bbox_inches="tight")