import nbformat as nbf
import os
import sys
import base64
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output', 'ovarian_results')
H5AD_FILE = os.path.join(DATA_DIR, 'ovary_abdomen_omentum_uterus_100k_whole_transcriptome_2025-11-08.h5ad')

def run_analysis_and_plot():
    print("Loading single-cell data...")
    adata = sc.read_h5ad(H5AD_FILE)
    genes_of_interest = ['HTR2A', 'HTR2C', 'TPH1', 'IDO1', 'TDO2']
    genes_present = [g for g in genes_of_interest if g in adata.var_names]
    
    adata_filtered = adata[adata.obs['cell_type'].isin(['malignant cell', 'T cell'])].copy()
    adata_filtered.obs['group'] = adata_filtered.obs['tissue_general'].astype(str) + " - " + adata_filtered.obs['cell_type'].astype(str)
    
    sc.pl.dotplot(adata_filtered, var_names=genes_present, groupby='group', 
                  standard_scale='var', cmap='Reds', title="Serotonin Axis Expression in Ovarian Cancer", 
                  show=False)
    plt.tight_layout()
    plot_path = os.path.join(DATA_DIR, 'serotonin_axis_dotplot.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    return plot_path, genes_present

def create_notebook():
    plot_path, genes_present = run_analysis_and_plot()
    
    nb = nbf.v4.new_notebook()
    
    nb.cells.append(nbf.v4.new_markdown_cell("""# Ovarian Serotonin Immune Evasion

### Goal
Determine whether ovarian peritoneal metastases exploit serotonin signaling to suppress local immune responses.

### Purpose
To computationally validate the hypothesis that HTR2A/HTR2C (serotonin receptors) and serotonin synthesis pathways (e.g., TPH1) are specifically upregulated in the ovarian metastatic niche (omentum/peritoneum), and to identify which cell types (Malignant cells vs. T cells) express these receptors.

### Interpretation
- **Upregulation in Metastasis:** If serotonin receptors and synthesis enzymes are highly expressed in the omentum/peritoneum compared to the primary ovary, it supports the serotonin-immune evasion hypothesis.
- **Cell-Type Specificity:** If T cells express HTR2A/HTR2C, they may be the direct targets of tumor-derived serotonin, leading to T cell suppression. If malignant cells express them, it may be an autocrine growth loop.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os

BASE_DIR = os.path.dirname(os.path.abspath('.'))
DATA_DIR = os.path.join(BASE_DIR, 'output', 'ovarian_results')
H5AD_FILE = os.path.join(DATA_DIR, 'ovary_abdomen_omentum_uterus_100k_whole_transcriptome_2025-11-08.h5ad')
"""))

    nb.cells.append(nbf.v4.new_markdown_cell("""### 1. Load Data and Define Genes of Interest
We are analyzing **HTR2A** and **HTR2C** (Serotonin Receptors) along with **TPH1** (Tryptophan Hydroxylase 1, the rate-limiting enzyme in serotonin synthesis). We also check **IDO1** / **TDO2** for general tryptophan depletion.
"""))

    code_load = """print("Loading single-cell data...")
adata = sc.read_h5ad(H5AD_FILE)
genes_of_interest = ['HTR2A', 'HTR2C', 'TPH1', 'IDO1', 'TDO2']

genes_present = [g for g in genes_of_interest if g in adata.var_names]
print(f"Genes found in dataset: {genes_present}")
"""
    cell_load = nbf.v4.new_code_cell(code_load)
    cell_load.outputs.append(nbf.v4.new_output("stream", text=f"Loading single-cell data...\\nGenes found in dataset: {genes_present}\\n", name="stdout"))
    nb.cells.append(cell_load)

    nb.cells.append(nbf.v4.new_markdown_cell("""### 2. Expression Across Cell Types and Metastatic Sites
We visualize the mean expression of these genes in Malignant cells vs. T cells, separated by the primary site (Ovary) and metastatic sites (Omentum / Abdomen).
"""))

    code_plot = """# Filter to cell types of interest
adata_filtered = adata[adata.obs['cell_type'].isin(['malignant cell', 'T cell'])].copy()

# Group by tissue and cell type
adata_filtered.obs['group'] = adata_filtered.obs['tissue_general'].astype(str) + " - " + adata_filtered.obs['cell_type'].astype(str)

sc.pl.dotplot(adata_filtered, var_names=genes_present, groupby='group', 
              standard_scale='var', cmap='Reds', title="Serotonin Axis Expression in Ovarian Cancer", 
              show=False)
plt.tight_layout()
plt.show()
"""
    cell_plot = nbf.v4.new_code_cell(code_plot)
    with open(plot_path, 'rb') as f:
        png_data = base64.b64encode(f.read()).decode('utf-8')
    cell_plot.outputs.append(nbf.v4.new_output("display_data", data={"image/png": png_data}))
    nb.cells.append(cell_plot)
    
    nb.cells.append(nbf.v4.new_markdown_cell("""### 3. Conclusion
The dotplot visualizes the expression patterns of serotonin receptors and synthesis enzymes across T cells and malignant cells in different anatomical sites. Note any significant upregulation of HTR2A/HTR2C in the omental/peritoneal T cells, which would validate the hypothesis that ovarian metastases exploit serotonin to suppress local T cell activation.
"""))

    export_code = """import subprocess
import sys

notebook_filename = 'ovarian_serotonin_immune_evasion.ipynb'
output_base = 'ovarian_serotonin_immune_evasion_5MetCan_100k'
output_dir = os.path.join('..', 'output', 'ovarian_results')
os.makedirs(output_dir, exist_ok=True)

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", "--execute", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{os.path.join(output_dir, output_base)}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""
    nb.cells.append(nbf.v4.new_code_cell(export_code))

    out_path = os.path.join(BASE_DIR, 'ovarian_serotonin_immune_evasion.ipynb')
    if os.path.exists(out_path):
        print(f"⚠️ {out_path} already exists. Skipping overwrite to protect user edits.")
    else:
        with open(out_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print(f"Created {out_path}")

if __name__ == '__main__':
    create_notebook()
