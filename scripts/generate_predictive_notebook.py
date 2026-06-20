import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
import nbformat as nbf
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def create_notebook():
    nb = nbf.v4.new_notebook()
    
    nb.cells.append(nbf.v4.new_markdown_cell("""# Predictive Potential of the Conserved Metabolic Signatures

### Inputs / Parameters
*This section documents explicit inputs for reproducibility.*
- **Configs:** `pan_cancer_config.py`
- **Data files:** Dynamically discovered pan-cancer combinations, specific cancer signatures, and primary tumor `.h5ad` files.

### Goal
Determine whether the dynamically discovered conserved metabolic signature genes (and cancer-specific metabolic signatures) are heterogeneously expressed in primary tumors. 

### Purpose
To compute single-cell "Metastatic Metabolic Scores" based on these signatures across malignant cells within primary tumors. By identifying sub-populations of primary tumor cells with high expression of these signatures, we hypothesize these represent pre-metastatic subclones. This notebook systematically tests the strict intersection, the relaxed subset combinations, and the cancer's own specific signature.
"""))

    code_setup = """import os
import sys
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_h5ad_path, get_de_csv_path

BASE_DIR = os.path.dirname(os.path.abspath('.'))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
os.makedirs(META_RESULTS_DIR, exist_ok=True)

conserved_genes_path = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
signature_genes_strict = []
if os.path.exists(conserved_genes_path):
    df_conserved = pd.read_csv(conserved_genes_path)
    if not df_conserved.empty and 'Strictly_Conserved_Gene' in df_conserved.columns:
        signature_genes_strict = df_conserved['Strictly_Conserved_Gene'].dropna().unique().tolist()

import glob
combo_files = glob.glob(os.path.join(META_RESULTS_DIR, f'pan_cancer_signature_*{ANALYSIS_SUFFIX}.csv'))
combo_signatures = {}
for cf in combo_files:
    cname = os.path.basename(cf).replace('pan_cancer_signature_', '').replace(f'{ANALYSIS_SUFFIX}.csv', '')
    df_cf = pd.read_csv(cf)
    if not df_cf.empty and 'Gene' in df_cf.columns:
        genes = df_cf['Gene'].dropna().unique().tolist()
        if len(genes) > 0:
            combo_signatures[cname] = genes

if len(signature_genes_strict) > 0:
    print(f"Loaded strict 5-cancer signature with {len(signature_genes_strict)} genes.")
else:
    print(f"Strict 5-cancer signature is empty. Falling back to {len(combo_signatures)} subset combinations.")
    if not combo_signatures:
        raise ValueError("No subset combinations found either. Cannot proceed.")
"""
    print_config_code = """
print('--- INJECTED PIPELINE CONFIGURATION ---')
from pan_cancer_config import CANCER_CAP, SKEW_THRESHOLD, SUBCLONE_SD_MULTIPLIER, ANALYSIS_SUFFIX
print(f'CANCER_CAP: {CANCER_CAP}')
print(f'SKEW_THRESHOLD: {SKEW_THRESHOLD}')
print(f'SUBCLONE_SD_MULTIPLIER: {SUBCLONE_SD_MULTIPLIER}')
print(f'ANALYSIS_SUFFIX: {ANALYSIS_SUFFIX}')
"""
    nb.cells.append(nbf.v4.new_code_cell(print_config_code))
    nb.cells.append(nbf.v4.new_code_cell(code_setup))

    code_analysis = """# Configuration for cancer datasets.
from pan_cancer_config import CANCER_CAP, CANCER_PRIMARY_TISSUE, CANCER_COLORS

cancers_config = []
for cancer in CANCER_CAP.keys():
    tissue = CANCER_PRIMARY_TISSUE.get(cancer, cancer)
    color = CANCER_COLORS.get(cancer, 'blue')
    title_name = cancer.capitalize()
    cancers_config.append((cancer, tissue, color, title_name))

summary_data = []

for prefix, tissue, color, title_name in cancers_config:
    print(f"\\n--- Analyzing {title_name} Cancer ---")
    h5ad_path = get_h5ad_path(prefix)
    if not os.path.exists(h5ad_path):
        raise FileNotFoundError(f"Missing {h5ad_path}")
        
    adata = sc.read_h5ad(h5ad_path)
    
    if isinstance(tissue, list):
        adata_pri = adata[(adata.obs['cell_type'] == 'malignant cell') & (adata.obs['tissue_general'].isin(tissue))].copy()
    else:
        adata_pri = adata[(adata.obs['cell_type'] == 'malignant cell') & (adata.obs['tissue_general'] == tissue)].copy()
        
    df_export = adata_pri.obs[['cell_type', 'tissue_general']].copy()

    # Determine which signatures apply to this cancer
    sigs_to_test = {}
    
    # 1. Load its OWN specific signature
    de_csv = get_de_csv_path(prefix)
    if os.path.exists(de_csv):
        df_de = pd.read_csv(de_csv)
        up_genes = df_de[df_de['Significance'] == 'Up in Metastasis']['names'].tolist()
        if len(up_genes) > 0:
            sigs_to_test['Own_Specific_Signature'] = up_genes
    else:
        print(f"Warning: No specific DE results found for {title_name} at {de_csv}.")

    # 2. Add pan-cancer signature subsets
    if len(signature_genes_strict) > 0:
        sigs_to_test['Strict_5Cancer'] = signature_genes_strict
    else:
        for cname, genes in combo_signatures.items():
            if title_name.lower() in cname.lower():
                sigs_to_test[cname] = genes
                
    if not sigs_to_test:
        print(f"No applicable signature combinations found for {title_name}. Skipping.")
        continue
        
    for sig_name, genes in sigs_to_test.items():
        valid_genes = [g for g in genes if g in adata_pri.var_names]
        if not valid_genes:
            print(f"No valid signature genes in dataset for {title_name} against {sig_name}.")
            continue
            
        score_col = f'Metastatic_Signature_Score_{sig_name}' if sig_name != 'Strict_5Cancer' else 'Metastatic_Signature_Score'
        
        sc.tl.score_genes(adata_pri, gene_list=valid_genes, score_name=score_col)
        df_export[score_col] = adata_pri.obs[score_col]
        
        plt.figure(figsize=(8,5))
        sns.histplot(adata_pri.obs[score_col], bins=50, kde=True, color=color)
        plt.title(f'Score Dist in Primary {title_name} vs {sig_name} ({len(valid_genes)} genes)')
        plt.xlabel('Signature Score')
        plt.ylabel('Cell Count')
        plt.axvline(adata_pri.obs[score_col].mean(), color='red', linestyle='dashed', label='Mean')
        plt.legend()
        
        plot_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_score_{sig_name}.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        scores = df_export[score_col]
        skew_val = scores.skew()
        if skew_val > SKEW_THRESHOLD: dist_type = "Right-skewed"
        elif skew_val < -SKEW_THRESHOLD: dist_type = "Left-skewed"
        else: dist_type = "Symmetric"
            
        pct_gt = (scores > (scores.mean() + (SUBCLONE_SD_MULTIPLIER * scores.std()))).mean() * 100
        
        summary_data.append({
            'Cancer': f"**{title_name}**",
            'Signature Combo': sig_name,
            'Primary Cells Scored': f"{len(adata_pri):,}",
            'Score Distribution': dist_type,
            'Pre-Metastatic Subclone (%)': f"~{pct_gt:.1f}%"
        })
        
    csv_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_scores{ANALYSIS_SUFFIX}.csv')
    df_export.to_csv(csv_path)
    print(f"Saved scores to {csv_path}")
    
df_subclone_summary = pd.DataFrame(summary_data)
summary_csv_path = os.path.join(META_RESULTS_DIR, f'pre_metastatic_subclone_summary{ANALYSIS_SUFFIX}.csv')
df_subclone_summary.to_csv(summary_csv_path, index=False)
print(f"Saved subclone summary to {summary_csv_path}")
"""
    nb.cells.append(nbf.v4.new_markdown_cell("### Scoring All Primary Cancer Cells"))
    nb.cells.append(nbf.v4.new_code_cell(code_analysis))

    code_summary = """summary_df = pd.DataFrame(summary_data)
summary_csv_path = os.path.join(META_RESULTS_DIR, f'pre_metastatic_subclone_summary{ANALYSIS_SUFFIX}.csv')
summary_df.to_csv(summary_csv_path, index=False)
print(f"\\nSaved summary table to {summary_csv_path}")
display(summary_df)
"""
    nb.cells.append(nbf.v4.new_markdown_cell("### Generating Final Summary Table"))
    nb.cells.append(nbf.v4.new_code_cell(code_summary))

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
print("To automatically export to HTML, run the above command in your terminal.")
"""
    nb.cells.append(nbf.v4.new_code_cell(export_code))

    out_path = os.path.join(BASE_DIR, 'predictive_signature_biomarker.ipynb')
    with open(out_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print(f"Created base notebook {out_path} with REAL analytical code!")

if __name__ == '__main__':
    create_notebook()
