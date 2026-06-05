# ── CELL 1: Imports ──────────────────────────────────────────
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ── CELL 2: Simulated qRT-PCR Data (6 biological replicates) ──
np.random.seed(42)

genes = ['CHAT', 'TH', 'BDNF']

data_rows = []
for gene, ctrl_mean, ra_mean, tr_mean in zip(
    genes,
    [1.0, 1.0, 1.0],
    [2.4, 0.6, 3.1],
    [1.8, 1.3, 4.2]
):
    for _ in range(6):
        data_rows.append([gene, 'Control', np.random.normal(ctrl_mean, 0.15)])
        data_rows.append([gene, 'RA_Diff', np.random.normal(ra_mean, 0.2)])
        data_rows.append([gene, 'Treated', np.random.normal(tr_mean, 0.2)])

df = pd.DataFrame(data_rows, columns=['Gene', 'Condition', 'Expression'])
print("Data shape:", df.shape)
print(df.head(9))

# ── CELL 3: Bar Chart ─────────────────────────────────────────
plt.figure(figsize=(8, 5))
sns.barplot(data=df, x='Gene', y='Expression', hue='Condition', capsize=0.05)
plt.title('Gene Expression Across Conditions', fontsize=13)
plt.ylabel('Relative Expression (fold change)', fontsize=11)
plt.xlabel('Gene', fontsize=11)
plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Baseline')
plt.legend()
plt.tight_layout()
plt.savefig('bargraph.png', dpi=150)
plt.show()
print("Bar chart saved as bargraph.png")

# ── CELL 4: Heatmap ───────────────────────────────────────────
heatmap_data = df.groupby(['Gene', 'Condition'])['Expression'].mean().unstack()

plt.figure(figsize=(7, 4))
sns.heatmap(heatmap_data,
            annot=True,
            fmt='.2f',
            cmap='RdYlGn',
            center=1.0,
            linewidths=0.5)
plt.title('Gene Expression Heatmap (Mean per Condition)', fontsize=12)
plt.tight_layout()
plt.savefig('heatmap.png', dpi=150)
plt.show()
print("Heatmap saved as heatmap.png")

# ── CELL 5: Statistical Testing ───────────────────────────────
def get_stars(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'

print("\nStatistical Results:")
for gene in genes:
    ctrl = df[(df['Gene'] == gene) & (df['Condition'] == 'Control')]['Expression']
    ra   = df[(df['Gene'] == gene) & (df['Condition'] == 'RA_Diff')]['Expression']
    tr   = df[(df['Gene'] == gene) & (df['Condition'] == 'Treated')]['Expression']

    _, p_ra = stats.ttest_ind(ctrl, ra)
    _, p_tr = stats.ttest_ind(ctrl, tr)

    print(f"  {gene} — Control vs RA_Diff: p={p_ra:.4f} {get_stars(p_ra)}")
    print(f"  {gene} — Control vs Treated:  p={p_tr:.4f} {get_stars(p_tr)}")

# ── CELL 6: Bar Chart with Significance Stars ─────────────────
fig, ax = plt.subplots(figsize=(9, 6))
sns.barplot(data=df, x='Gene', y='Expression', hue='Condition',
            capsize=0.05, ax=ax)

ax.set_title('Gene Expression with Statistical Significance', fontsize=13)
ax.set_ylabel('Relative Expression (fold change)', fontsize=11)
ax.set_xlabel('Gene', fontsize=11)
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Baseline')

offsets = [-0.27, 0, 0.27]

for i, gene in enumerate(genes):
    ctrl = df[(df['Gene'] == gene) & (df['Condition'] == 'Control')]['Expression']
    ra   = df[(df['Gene'] == gene) & (df['Condition'] == 'RA_Diff')]['Expression']
    tr   = df[(df['Gene'] == gene) & (df['Condition'] == 'Treated')]['Expression']

    _, p_ra = stats.ttest_ind(ctrl, ra)
    _, p_tr = stats.ttest_ind(ctrl, tr)

    max_val = df[df['Gene'] == gene]['Expression'].max()

    ax.text(i + offsets[1], max_val + 0.3, get_stars(p_ra),
            ha='center', fontsize=13, color='steelblue', fontweight='bold')
    ax.text(i + offsets[2], max_val + 0.6, get_stars(p_tr),
            ha='center', fontsize=13, color='coral', fontweight='bold')

ax.legend()
plt.tight_layout()
plt.savefig('bargraph_with_stars.png', dpi=150)
plt.show()
print("Bar chart with stars saved as bargraph_with_stars.png")

print("\nAll done. Three figures generated and saved.")