import nbformat as nbf
import os

def create_notebook():
    nb = nbf.v4.new_notebook()
    
    nb.cells.append(nbf.v4.new_markdown_cell("""# MITF Regulon Expansion Across Cancer Types
## Expanding the Metabolic Regulon of MITF
This notebook addresses **Priority 1** (Step 3) from the AI Summary (Version 5).
**Purpose**: Assess MITF binding across the 1,669 metabConnectomeDB target universe using ChEA 2022 ChIP-Seq data.
**Interpretation**: Reveals whether MITF controls a broader metabolic network active in non-melanoma cancers.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import os
import sys
import pandas as pd

sys.path.append(os.path.abspath("."))
from compute_mitf_regulon import compute_mitf_regulon

compute_mitf_regulon()
# Load output
mitf_out = os.path.join("..", "output", "mitf_regulon", "mitf_metabolic_regulon_pairs.csv")
if os.path.exists(mitf_out):
    df = pd.read_csv(mitf_out)
    display(df.head())
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import matplotlib.pyplot as plt
import seaborn as sns

if os.path.exists(mitf_out) and 'Metabolite' in df.columns:
    plt.figure(figsize=(10,6))
    top_metabolites = df['Metabolite'].value_counts().head(15)
    sns.barplot(y=top_metabolites.index, x=top_metabolites.values, palette='viridis')
    plt.title('Top 15 Metabolites Interacting with MITF Transcriptional Targets')
    plt.xlabel('Number of Target Genes')
    plt.ylabel('Metabolite')
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import subprocess
import sys
import os

notebook_filename = 'mitf_regulon_expansion.ipynb'
output_base = 'mitf_regulon_expansion'
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

    out_path = os.path.join(os.path.dirname(__file__), "mitf_regulon_expansion.ipynb")
    if os.path.exists(out_path):
        print(f"⚠️ {out_path} already exists. Skipping overwrite to protect user edits.")
    else:
        with open(out_path, 'w') as f:
            nbf.write(nb, f)
        print(f"Created {out_path}")

if __name__ == "__main__":
    create_notebook()
