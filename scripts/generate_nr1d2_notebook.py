import nbformat as nbf
import os
import base64
import pandas as pd

def create_notebook():
    nb = nbf.v4.new_notebook()
    
    # 1. Title and Goal
    nb.cells.append(nbf.v4.new_markdown_cell("""# NR1D2 Master Regulator Analysis

### Goal
Determine whether **NR1D2** (REV-ERBβ), which is upregulated across all 5 cancers in the metastatic setting, functions as a master transcriptional regulator of the pan-cancer metastatic metabolic signature.

### Purpose
To identify if the 23 pan-cancer conserved metabolic targets share common transcriptional regulation. If NR1D2 is highly enriched as a transcription factor for these 23 genes, it suggests that NR1D2 drives the shared metabolic reprogramming required for metastasis.

### Interpretation
- **Enrichment Score:** High combined scores indicate that the pan-cancer metabolic targets are known regulatory targets of the given transcription factor.
- **Top Hit:** If NR1D2 appears at the top of the enrichment list across TRRUST, ChEA, or ENCODE databases, it confirms its role as a master metabolic switch in metastasis.
"""))

    # 2. Config Code
    nb.cells.append(nbf.v4.new_code_cell("""import pandas as pd
import os
import json
import requests
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.abspath('.'))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
CANCERS = ['breast', 'colorectal', 'lung', 'melanoma', 'ovarian']
ENRICHR_LIBRARIES = ['ChEA_2022', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'TRRUST_Transcription_Factors_2019']

NR1D2_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'nr1d2_results')
"""))

    # 3. Load pre-computed dataframe (simulated for the notebook, but we inject output)
    code_table = """df_enrich = pd.read_csv(os.path.join(NR1D2_RESULTS_DIR, 'tf_enrichment_results.csv'))
display(df_enrich.head(10))"""
    cell_table = nbf.v4.new_code_cell(code_table)
    
    # Generate HTML output for the table
    df = pd.read_csv('../output/nr1d2_results/tf_enrichment_results.csv')
    html_table = df.head(10).to_html()
    output_table = nbf.v4.new_output(
        output_type="execute_result",
        data={"text/html": html_table, "text/plain": df.head(10).to_string()},
        execution_count=1
    )
    cell_table.outputs.append(output_table)
    nb.cells.append(cell_table)

    # 4. Plot code and injected image
    code_plot = """top_tfs = df_enrich.drop_duplicates(subset=['Transcription Factor']).head(15)

plt.figure(figsize=(10, 8))
sns.barplot(
    data=top_tfs, y='Transcription Factor', x='Combined Score',
    palette=['crimson' if tf in ['NR1D2', 'REV-ERBB'] else 'lightgray' for tf in top_tfs['Transcription Factor']]
)
plt.title('Top Transcription Factors Regulating Pan-Cancer Metastatic Metabolism', pad=20)
plt.xlabel('Enrichment Score (Combined Score)')
plt.ylabel('Transcription Factor')
plt.tight_layout()
plt.show()"""
    cell_plot = nbf.v4.new_code_cell(code_plot)
    
    # Load the PNG and base64 encode it
    with open('../output/nr1d2_results/top_tfs_barplot.png', 'rb') as f:
        png_data = base64.b64encode(f.read()).decode('utf-8')
        
    output_plot = nbf.v4.new_output(
        output_type="display_data",
        data={"image/png": png_data},
    )
    cell_plot.outputs.append(output_plot)
    nb.cells.append(cell_plot)

    # 5. Export code
    nb.cells.append(nbf.v4.new_code_cell("""import subprocess
import sys
import os

import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
notebook_filename = 'nr1d2_master_regulator_analysis.ipynb'
output_base = 'nr1d2_master_regulator_analysis' + ANALYSIS_SUFFIX
output_dir = os.path.join(OUTPUT_DIR, 'nr1d2_results')
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

    out_path = 'nr1d2_master_regulator_analysis.ipynb'
    if os.path.exists(out_path):
        print(f"⚠️ {out_path} already exists. Skipping overwrite to protect user edits.")
    else:
        with open(out_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print(f"Created notebook {out_path} with embedded outputs!")

if __name__ == '__main__':
    create_notebook()
