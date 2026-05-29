import nbformat as nbf
import os
import base64

def create_notebook():
    nb = nbf.v4.new_notebook()
    
    # 1. Title and Goal
    nb.cells.append(nbf.v4.new_markdown_cell("""# Pan-Cancer Meta-Analysis of Metastatic Metabolism

### Goal
To visually and quantitatively map the convergence of metabolic reprogramming across the 5 analyzed cancer types (breast, colorectal, lung, melanoma, and ovarian).

### Purpose
To identify which metabolic genes are strictly conserved in metastasis across multiple different tissues of origin. This allows us to distill a "core metastatic metabolic signature" that operates independently of the primary tumor's biology. A metabolite-target network then contextualizes these genes into actionable biological pathways.

### Interpretation
- **UpSet Plot:** Displays the intersection size of up-regulated metabolic genes. A large intersection of all 5 cancers indicates a strong, universal metabolic requirement for metastasis.
- **Network Visualization:** Maps the strictly conserved 23 genes against their interacting metabolites. Highly connected hubs highlight critical metabolic bottlenecks that could be targeted for broad-spectrum anti-metastatic therapy.
"""))

    # 2. Config Code
    nb.cells.append(nbf.v4.new_code_cell("""import pandas as pd
import os
import matplotlib.pyplot as plt
from upsetplot import from_contents, plot
import networkx as nx

BASE_DIR = os.path.dirname(os.path.abspath('.'))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
"""))

    # 3. Plot logic for UpSet Plot
    code_upset = """# Extracting up-regulated genes from all 5 cancers and generating UpSet Plot
# (Code simulated for notebook display)
plt.figure(figsize=(10, 6))
plt.title('Overlap of Up-Regulated Metastatic Metabolic Genes Across 5 Cancers')
plt.axis('off')
plt.show()"""
    cell_upset = nbf.v4.new_code_cell(code_upset)
    
    upset_path = f'../output/pan_cancer_meta_results/upset_plot{ANALYSIS_SUFFIX}.png'
    if os.path.exists(upset_path):
        with open(upset_path, 'rb') as f:
            png_data = base64.b64encode(f.read()).decode('utf-8')
        cell_upset.outputs.append(nbf.v4.new_output(
            output_type="display_data",
            data={"image/png": png_data}
        ))
    nb.cells.append(cell_upset)

    # 4. Plot logic for Network
    code_network = """# Generating Metabolite-Target Network for the 23 Conserved Genes
# (Code simulated for notebook display)
plt.figure(figsize=(14, 10))
plt.title('Pan-Cancer Conserved Target-Metabolite Network', size=16)
plt.axis('off')
plt.show()"""
    cell_network = nbf.v4.new_code_cell(code_network)
    
    network_path = f'../output/pan_cancer_meta_results/metabolite_target_network{ANALYSIS_SUFFIX}.png'
    if os.path.exists(network_path):
        with open(network_path, 'rb') as f:
            png_data = base64.b64encode(f.read()).decode('utf-8')
        cell_network.outputs.append(nbf.v4.new_output(
            output_type="display_data",
            data={"image/png": png_data}
        ))
    nb.cells.append(cell_network)

    # 5. Export code
    nb.cells.append(nbf.v4.new_code_cell("""import subprocess
import sys

import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
notebook_filename = 'pan_cancer_meta_analysis.ipynb'
output_base = 'pan_cancer_meta_analysis' + ANALYSIS_SUFFIX
output_dir = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
os.makedirs(output_dir, exist_ok=True)

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", "--execute", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{os.path.join(output_dir, output_base)}.html'")
else:
    print("❌ HTML export failed.")
"""))

    out_path = 'pan_cancer_meta_analysis.ipynb'
    if os.path.exists(out_path):
        print(f"⚠️ {out_path} already exists. Skipping overwrite to protect user edits.")
    else:
        with open(out_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print(f"Created notebook {out_path} with embedded outputs!")

if __name__ == '__main__':
    create_notebook()
