# ── CELL 1: Imports ──────────────────────────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ── CELL 2: Load Data ─────────────────────────────────────────
filepath = r"C:\data d\CEBS\Vs code python\NDA2\GSE5281_series_matrix.txt.gz"

df_raw = pd.read_csv(filepath,
                     sep='\t',
                     comment='!',
                     compression='gzip',
                     index_col=0)

print("Shape:", df_raw.shape)

# ── CELL 3: Split Groups & Statistics ─────────────────────────
sample_names = df_raw.columns.tolist()
control_samples = sample_names[:74]
ad_samples = sample_names[74:]

df_raw['Control_Mean'] = df_raw[control_samples].mean(axis=1)
df_raw['AD_Mean'] = df_raw[ad_samples].mean(axis=1)
df_raw['Fold_Change'] = df_raw['AD_Mean'] / df_raw['Control_Mean']
df_raw['Log2_FC'] = np.log2(df_raw['Fold_Change'])

pvalues = []
for idx in df_raw.index:
    ctrl_vals = df_raw.loc[idx, control_samples].values.astype(float)
    ad_vals = df_raw.loc[idx, ad_samples].values.astype(float)
    _, p = stats.ttest_ind(ctrl_vals, ad_vals)
    pvalues.append(p)

df_raw['p_value'] = pvalues
df_raw['-log10_p'] = -np.log10(df_raw['p_value'])

# ── CELL 4: Volcano Plot ───────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))

colors = []
for idx in df_raw.index:
    fc = df_raw.loc[idx, 'Log2_FC']
    p = df_raw.loc[idx, 'p_value']
    if p < 0.05 and fc < -1:
        colors.append('steelblue')
    elif p < 0.05 and fc > 1:
        colors.append('coral')
    else:
        colors.append('lightgray')

ax.scatter(df_raw['Log2_FC'], df_raw['-log10_p'],
           c=colors, alpha=0.4, s=5)

ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
ax.axvline(x=1, color='black', linestyle='--', linewidth=0.8)
ax.axvline(x=-1, color='black', linestyle='--', linewidth=0.8)

ax.set_xlabel('Log2 Fold Change (AD vs Control)', fontsize=12)
ax.set_ylabel('-log10(p-value)', fontsize=12)
ax.set_title("Volcano Plot — Alzheimer's vs Control Brain (GSE5281)", fontsize=13)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='coral',
           markersize=8, label='Upregulated in AD'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='steelblue',
           markersize=8, label='Downregulated in AD'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgray',
           markersize=8, label='Not significant'),
]
ax.legend(handles=legend_elements, fontsize=10)
plt.tight_layout()
plt.show()

# ── CELL 5: Top Downregulated Genes ───────────────────────────
top_down = df_raw[df_raw['Log2_FC'] < -1].nsmallest(20, 'p_value')[['Log2_FC', 'p_value']]

probe_to_gene = {
    '227404_s_at': 'HSPA2',
    '230656_s_at': 'NPTX1',
    '203122_at': 'SNAP25',
    '201400_at': 'ATP2B2',
    '219549_s_at': 'STXBP6',
    '202517_at': 'SYT1',
    '213726_x_at': 'GABRB2',
    '202712_s_at': 'VAMP2',
    '207232_s_at': 'KCNC2',
    '221476_s_at': 'CBLN1',
    '221772_s_at': 'NRXN3',
    '208977_x_at': 'GABRA1',
    '211750_x_at': 'GABRB2',
    '209251_x_at': 'GABRA1',
    '213646_x_at': 'NRXN1',
    '209755_at': 'KCNA1',
    '208678_at': 'SYP',
    '201199_s_at': 'CALM2',
    '227944_at': 'NRXN3',
    '201387_s_at': 'GRIN2B'
}

top_down['Gene'] = top_down.index.map(probe_to_gene)

# ── CELL 6: Final Bar Plot ─────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))
top_down_named = top_down.dropna(subset=['Gene'])

colors = ['steelblue' if fc < -1.5 else 'cornflowerblue'
          for fc in top_down_named['Log2_FC']]

ax.barh(top_down_named['Gene'],
        top_down_named['Log2_FC'],
        color=colors, edgecolor='white')

ax.axvline(x=0, color='black', linewidth=0.8)
ax.axvline(x=-1, color='red', linestyle='--',
           linewidth=0.8, label='FC threshold (-1)')

ax.set_xlabel("Log2 Fold Change (AD vs Control)", fontsize=12)
ax.set_title("Top Downregulated Synaptic Genes in Alzheimer's Disease\n(GSE5281 — Human Brain Tissue)",
             fontsize=12)

for i, (idx, row) in enumerate(top_down_named.iterrows()):
    p = row['p_value']
    if p < 0.001: stars = '***'
    elif p < 0.01: stars = '**'
    else: stars = '*'
    ax.text(row['Log2_FC'] - 0.05, i, stars,
            va='center', ha='right', fontsize=9, color='white', fontweight='bold')

ax.legend()
ax.invert_yaxis()
plt.tight_layout()
plt.show()

print("\nAnalysis complete.")