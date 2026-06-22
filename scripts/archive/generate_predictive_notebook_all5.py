"""
Purpose: Scores the primary tumor cells of each cancer against their OWN specific metastatic signature to identify highly metastatic ("Right-Skewed") subclones present before dissemination.
"""
import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_h5ad_path
import os
import glob
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
os.makedirs(META_RESULTS_DIR, exist_ok=True)

def run_analysis_and_plot():
    from pan_cancer_config import get_de_csv_path
    
    # We will score each cancer against its OWN upregulated metastatic genes, 
    # not the pan-cancer signature, as requested by the user.

    # Configuration for cancer datasets.
    # Format: (h5ad_prefix, tissue_general_filter, plot_color, display_name)
    # The display_name (4th value) is the proper capitalized name used in plots and reports.
    from pan_cancer_config import CANCER_CAP, CANCER_PRIMARY_TISSUE, CANCER_COLORS
    
    cancers_config = []
    for cancer in CANCER_CAP.keys():
        tissue = CANCER_PRIMARY_TISSUE.get(cancer, cancer)
        color = CANCER_COLORS.get(cancer, 'blue')
        title_name = cancer.capitalize()
        cancers_config.append((cancer, tissue, color, title_name))
    
    for prefix, tissue, color, title_name in cancers_config:
        h5ad_path = get_h5ad_path(prefix)
        adata = sc.read_h5ad(h5ad_path)
        
        tumor_mask = adata.obs['cell_type'].str.contains('malignant|tumor|epithelial|cancer', case=False, na=False)
        
        if isinstance(tissue, list):
            adata_pri = adata[tumor_mask & (adata.obs['tissue_general'].isin(tissue))].copy()
        else:
            adata_pri = adata[tumor_mask & (adata.obs['tissue_general'] == tissue)].copy()
            
        df_export = adata_pri.obs[['cell_type', 'tissue_general']].copy()
        
        # Load THIS cancer's specific metastatic DE genes
        de_csv = get_de_csv_path(prefix)
        if not os.path.exists(de_csv):
            print(f"No DE results found for {title_name} at {de_csv}. Skipping.")
            continue
            
        df_de = pd.read_csv(de_csv)
        up_genes = df_de[df_de['Significance'] == 'Up in Metastasis']['names'].tolist()
        
        valid_genes = [g for g in up_genes if g in adata_pri.var_names]
        if not valid_genes:
            print(f"No valid genes found in {title_name} dataset from its own signature. Skipping.")
            continue
            
        sc.tl.score_genes(adata_pri, gene_list=valid_genes, score_name='Metastatic_Signature_Score')
        
        df_export['Metastatic_Signature_Score'] = adata_pri.obs['Metastatic_Signature_Score']
        
        # Plot
        plt.figure(figsize=(8,5))
        sns.histplot(adata_pri.obs['Metastatic_Signature_Score'], bins=50, kde=True, color=color)
        plt.title(f'{title_name} Primary Tumor against its {len(valid_genes)}-Gene Specific Signature')
        plt.xlabel(f'Signature Score ({len(valid_genes)} genes)')
        plt.ylabel('Cell Count')
        plt.axvline(adata_pri.obs['Metastatic_Signature_Score'].mean(), color='red', linestyle='dashed', label='Mean')
        plt.legend()

        plot_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_score.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        csv_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_scores{ANALYSIS_SUFFIX}.csv')
        df_export.to_csv(csv_path)
        print(f"Saved {title_name} scores to {csv_path}")

if __name__ == '__main__':
    run_analysis_and_plot()
