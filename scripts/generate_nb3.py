import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
import nbformat as nbf
import os

def create_notebook():
    nb = nbf.v4.new_notebook()
    
    nb.cells.append(nbf.v4.new_markdown_cell("""# Serotonin Axis Spatial Mapping in Ovarian Metastasis
## Investigating Juxtacrine vs Paracrine T-Cell Silencing
This notebook addresses **Priority 2** (Step 6) from the AI Summary (Version 5).
**Purpose**: Compute "spatial proximity scores" between TPH1-expressing tumor clusters and HTR2A-expressing T-cell clusters in the ovarian metastatic niche.
**Interpretation**: Determines whether serotonin-mediated immune evasion requires direct cell-cell contact (juxtacrine) or acts via diffusion (paracrine).
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import os
import sys
import pandas as pd

sys.path.append(os.path.abspath("."))
from compute_serotonin_spatial import compute_serotonin_spatial

compute_serotonin_spatial()

s_out = os.path.join("..", "output", "serotonin_axis_spatial_mapping", f"serotonin_proximity_results{ANALYSIS_SUFFIX}.csv")
if os.path.exists(s_out):
    df = pd.read_csv(s_out)
    display(df)
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import matplotlib.pyplot as plt
import seaborn as sns

if os.path.exists(s_out):
    counts = [df['TPH1_Tumor_Count'].iloc[0], df['HTR2A_Tcell_Count'].iloc[0]]
    labels = ['TPH1+ Tumor Cells', 'HTR2A+ T Cells']
    plt.figure(figsize=(6,4))
    sns.barplot(x=labels, y=counts, palette='pastel')
    plt.title('Serotonin Axis: Key Cell Type Counts in Ovarian Microenvironment')
    plt.ylabel('Number of Cells')
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import subprocess
import sys
import os

notebook_filename = 'serotonin_axis_spatial_mapping.ipynb'
output_base = 'serotonin_axis_spatial_mapping'
output_dir = os.path.join('..', 'output', 'jupyter_reports')
os.makedirs(output_dir, exist_ok=True)

cmd_html = [sys.executable, "-m", "jupyter", "nbconvert", "--to", "html", "--execute", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{os.path.join(output_dir, output_base)}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""))

    out_path = os.path.join(os.path.dirname(__file__), "serotonin_axis_spatial_mapping.ipynb")
    if os.path.exists(out_path):
        print(f"⚠️ {out_path} already exists. Skipping overwrite to protect user edits.")
    else:
        with open(out_path, 'w') as f:
            nbf.write(nb, f)
        print(f"Created {out_path}")

if __name__ == "__main__":
    create_notebook()
