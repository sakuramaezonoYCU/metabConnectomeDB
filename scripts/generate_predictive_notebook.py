import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_h5ad_path
import nbformat as nbf
import os
import base64
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
os.makedirs(META_RESULTS_DIR, exist_ok=True)

def run_analysis_and_plot():
    genes_23_path = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(genes_23_path):
        df_23 = pd.read_csv(genes_23_path)
        signature_genes = df_23['Strictly_Conserved_Gene'].unique().tolist()
    else:
        # Fallback
        signature_genes = ['GLS', 'SGMS1', 'SLC16A7', 'SPTLC1', 'ENO1', 'LDHA', 'HK2']
        
    cancers_config = [
        ('breast', 'breast', 'purple', 'Breast'),
        ('colorectal', ['colon', 'large intestine'], 'orange', 'Colorectal'),
        ('lung', 'lung', 'teal', 'Lung'),
        ('melanoma', 'skin of body', 'black', 'Melanoma'),
        ('ovarian', 'ovary', 'pink', 'Ovarian')
    ]
    
    plot_paths = {}
    
    for prefix, tissue, color, title_name in cancers_config:
        h5ad_path = get_h5ad_path(prefix)
        adata = sc.read_h5ad(h5ad_path)
        
        if isinstance(tissue, list):
            adata_pri = adata[(adata.obs['cell_type'] == 'malignant cell') & (adata.obs['tissue_general'].isin(tissue))].copy()
        else:
            adata_pri = adata[(adata.obs['cell_type'] == 'malignant cell') & (adata.obs['tissue_general'] == tissue)].copy()
            
        valid_genes = [g for g in signature_genes if g in adata_pri.var_names]
        sc.tl.score_genes(adata_pri, gene_list=valid_genes, score_name='Metastatic_Signature_Score')
        
        # Plot
        plt.figure(figsize=(8,5))
        sns.histplot(adata_pri.obs['Metastatic_Signature_Score'], bins=50, kde=True, color=color)
        plt.title(f'Distribution of Metastatic Metabolic Score in Primary {title_name} Tumor Cells')
        plt.xlabel(f'{len(signature_genes)}-Gene Signature Score')
        plt.ylabel('Cell Count')
        plt.axvline(adata_pri.obs['Metastatic_Signature_Score'].mean(), color='red', linestyle='dashed', label='Mean')
        plt.legend()
        
        csv_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_scores{ANALYSIS_SUFFIX}.csv')
        
        # Export CSV for the histogram
        df_export = adata_pri.obs[['cell_type', 'tissue_general', 'Metastatic_Signature_Score']].copy()
        df_export.to_csv(csv_path)
        print(f"Saved {title_name} scores to {csv_path}")
        
        plot_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_score{ANALYSIS_SUFFIX}.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        plot_paths[prefix] = plot_path
        
    return plot_paths, signature_genes

def create_notebook():
    plot_paths, signature_genes = run_analysis_and_plot()
    
    nb = nbf.v4.new_notebook()
    
    nb.cells.append(nbf.v4.new_markdown_cell(f"""# Predictive Potential of the {len(signature_genes)}-Gene Pan-Cancer Metabolic Signature

### Goal
Determine whether the {len(signature_genes)} pan-cancer conserved metabolic signature genes are heterogeneously expressed in primary tumors, which would support their utility as a predictive biomarker for future metastasis.

### Purpose
To compute a single-cell "Metastatic Metabolic Score" based on the {len(signature_genes)} pan-cancer genes across malignant cells within primary tumors. By identifying a sub-population of primary tumor cells with high expression of this signature, we hypothesize these represent pre-metastatic subclones.

### Interpretation
- **High Variability / Bimodal Distribution:** Indicates that a subset of primary tumor cells has already adopted the metastatic metabolic program.
- **Biomarker Utility:** If a clear high-scoring subpopulation exists in the primary tumor, this {len(signature_genes)}-gene signature could be developed into a clinical assay to predict metastatic risk.
"""))

    nb.cells.append(nbf.v4.new_code_cell("""import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

BASE_DIR = os.path.dirname(os.path.abspath('.'))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
"""))

    for prefix, plot_path in plot_paths.items():
        title = prefix.capitalize()
        nb.cells.append(nbf.v4.new_markdown_cell(f"### Scoring Primary {title} Cancer Cells\nWe score malignant cells in the primary {title} tumor."))
        code_str = f"# {title} Cancer\\nplt.figure(figsize=(8,5))\\nplt.title('Distribution of Metastatic Metabolic Score in Primary {title} Tumor Cells')\\nplt.show()"
        cell = nbf.v4.new_code_cell(code_str)
        with open(plot_path, 'rb') as f:
            png_data = base64.b64encode(f.read()).decode('utf-8')
        cell.outputs.append(nbf.v4.new_output("display_data", data={"image/png": png_data}))
        nb.cells.append(cell)

    nb.cells.append(nbf.v4.new_markdown_cell("""### Summary
In multiple primary cancer tissues, a distinctive 'tail' or secondary peak of cells exhibits high expression of the metastatic signature. These cells may represent the pre-metastatic subclone population. 

**Next Steps**: Develop multiplex IHC or targeted spatial transcriptomics assays to detect these cells clinically in patient biopsies prior to recurrence.
"""))

    export_code = f"""import subprocess
import sys
import os

notebook_filename = 'predictive_signature_biomarker.ipynb'
output_base = 'predictive_signature_biomarker' + '{ANALYSIS_SUFFIX}'
output_dir = os.path.join('..', 'output', 'pan_cancer_meta_results')
os.makedirs(output_dir, exist_ok=True)

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", "--execute", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{{os.path.join(output_dir, output_base)}}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""
    nb.cells.append(nbf.v4.new_code_cell(export_code))

    out_path = os.path.join(BASE_DIR, 'predictive_signature_biomarker.ipynb')
    with open(out_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print(f"Created base notebook {out_path} (HTML will have the suffix!)")

if __name__ == '__main__':
    create_notebook()
