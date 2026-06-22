import json
import os
import sys

nb_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts', 'oxygen_tension_analysis.ipynb')
with open(nb_path, 'r') as f:
    nb = json.load(f)

new_code = """plt.figure(figsize=(8, 6))

# Clean NaN if any
df_plot_hif1 = df_res.dropna(subset=['O2_Tension_Tumour_Pct', 'Mean_HIF1_LFC']).copy()

# Scatter plot
sns.scatterplot(
    data=df_plot_hif1,
    x='O2_Tension_Tumour_Pct',
    y='Mean_HIF1_LFC',
    hue='Cancer',
    s=200,
    palette='Set1'
)

# Fit line
sns.regplot(
    data=df_plot_hif1,
    x='O2_Tension_Tumour_Pct',
    y='Mean_HIF1_LFC',
    scatter=False,
    color='grey',
    line_kws={"linestyle":"--"}
)

plt.title("Oxygen Tension vs Mean HIF-1 Log2FC", fontsize=14, weight='bold')
plt.xlabel("Metastatic Site Oxygen Tension (%)", fontsize=12)
plt.ylabel("Mean HIF-1 LFC (Metastasis vs Primary)", fontsize=12)
plt.axhline(0.0, color='black', linewidth=1, linestyle=':')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Calculate stats
r_hif1, p_hif1 = stats.pearsonr(df_plot_hif1['O2_Tension_Tumour_Pct'], df_plot_hif1['Mean_HIF1_LFC'])
plt.annotate(f"Pearson r: {r_hif1:.2f}\\np-value: {p_hif1:.3f}", 
             xy=(0.05, 0.85), xycoords='axes fraction', 
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

plt.tight_layout()

# output path
out_dir = os.path.join(os.path.abspath(os.path.join(os.getcwd(), '..')), "output", "oxygen_tension")
os.makedirs(out_dir, exist_ok=True)
filename_hif1 = "oxygen_tension_hif1_correlation"+ ANALYSIS_SUFFIX + ".png"
output_plot_hif1 = os.path.join(out_dir, filename_hif1)
plt.savefig(output_plot_hif1, dpi=300, bbox_inches='tight')
print(f"Saved plot to {output_plot_hif1}")
plt.show()
"""

# Find the cell that has plt.show() and sns.scatterplot
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code' and 'sns.scatterplot(' in "".join(cell.get('source', [])):
        # We append the new code to this cell
        source = cell['source']
        if not source[-1].endswith('\n'):
            source[-1] += '\n'
        source.append("\n")
        # Split new code into lines with newlines
        lines = [line + '\n' for line in new_code.split('\n')]
        # remove last newline to match format
        lines[-1] = lines[-1].strip('\n')
        source.extend(lines)
        break

with open(nb_path, 'w') as f:
    json.dump(nb, f, indent=1)
    
print("Added HIF1 plot to notebook!")
