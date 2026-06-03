import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_h5ad_path
import os
import sys
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
        # Fallback if not run yet
        signature_genes = ['GLS', 'SGMS1', 'SLC16A7', 'SPTLC1', 'ENO1', 'LDHA', 'HK2']
        
    cancers_config = [
        ('breast', 'breast', 'purple', 'Breast'),
        ('lung', 'lung', 'teal', 'Lung'),
        ('colorectal', ['colon', 'large intestine'], 'orange', 'Colorectal'),
        ('melanoma', 'skin of body', 'black', 'Melanoma'),
        ('ovarian', 'ovary', 'pink', 'Ovarian')
    ]
    
    for prefix, tissue, color, title_name in cancers_config:
        h5ad_path = get_h5ad_path(prefix)
        adata = sc.read_h5ad(h5ad_path)
        
        tumor_mask = adata.obs['cell_type'].str.contains('malignant|tumor|epithelial|cancer', case=False, na=False)
        
        if isinstance(tissue, list):
            adata_pri = adata[tumor_mask & (adata.obs['tissue_general'].isin(tissue))].copy()
        else:
            adata_pri = adata[tumor_mask & (adata.obs['tissue_general'] == tissue)].copy()
            
        valid_genes = [g for g in signature_genes if g in adata_pri.var_names]
        sc.tl.score_genes(adata_pri, gene_list=valid_genes, score_name='Metastatic_Signature_Score')
        
        # Export CSV for the histogram
        df_export = adata_pri.obs[['cell_type', 'tissue_general', 'Metastatic_Signature_Score']].copy()
        csv_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_scores.csv')
        df_export.to_csv(csv_path)
        print(f"Saved {title_name} scores to {csv_path}")
        
        # Plot
        plt.figure(figsize=(8,5))
        sns.histplot(adata_pri.obs['Metastatic_Signature_Score'], bins=50, kde=True, color=color)
        plt.title(f'Distribution of Metastatic Metabolic Score in Primary {title_name} Tumor Cells')
        plt.xlabel('21-Gene Signature Score')
        plt.ylabel('Cell Count')
        plt.axvline(adata_pri.obs['Metastatic_Signature_Score'].mean(), color='red', linestyle='dashed', label='Mean')
        plt.legend()

        plot_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_score.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    run_analysis_and_plot()
