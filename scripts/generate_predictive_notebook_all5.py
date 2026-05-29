import nbformat as nbf
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
    genes_23_path = os.path.join(META_RESULTS_DIR, 'pan_cancer_23_genes.csv')
    if os.path.exists(genes_23_path):
        df_23 = pd.read_csv(genes_23_path)
        signature_genes = df_23['Strictly_Conserved_Gene'].unique().tolist()
    else:
        # Fallback if not run yet
        signature_genes = ['GLS', 'SGMS1', 'SLC16A7', 'SPTLC1', 'ENO1', 'LDHA', 'HK2']
        
    cancers_config = [
        ('breast', 'breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad', 'breast', 'purple', 'Breast'),
        ('lung', 'lung_lymph-node_brain_pleural-fluid_100k_whole_transcriptome_2025-11-08.h5ad', 'lung', 'teal', 'Lung'),
        ('colorectal', 'colon_large-intestine_liver_intestine_lung_100k_whole_transcriptome_2025-11-08.h5ad', ['colon', 'large intestine'], 'orange', 'Colorectal'),
        ('melanoma', 'skin-of-body_brain_abdomen_paracolic-gutter_100k_whole_transcriptome_2025-11-08.h5ad', 'skin of body', 'black', 'Melanoma'),
        ('ovarian', 'ovary_abdomen_omentum_uterus_100k_whole_transcriptome_2025-11-08.h5ad', 'ovary', 'pink', 'Ovarian')
    ]
    
    for prefix, h5ad_name, tissue, color, title_name in cancers_config:
        h5ad_path = os.path.join(OUTPUT_DIR, f"{prefix}_results", h5ad_name)
        adata = sc.read_h5ad(h5ad_path)
        
        if isinstance(tissue, list):
            adata_pri = adata[(adata.obs['cell_type'] == 'malignant cell') & (adata.obs['tissue_general'].isin(tissue))].copy()
        else:
            adata_pri = adata[(adata.obs['cell_type'] == 'malignant cell') & (adata.obs['tissue_general'] == tissue)].copy()
            
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
        plt.xlabel('23-Gene Signature Score')
        plt.ylabel('Cell Count')
        plt.axvline(adata_pri.obs['Metastatic_Signature_Score'].mean(), color='red', linestyle='dashed', label='Mean')
        plt.legend()
        plot_path = os.path.join(META_RESULTS_DIR, f'{prefix}_primary_signature_score.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    run_analysis_and_plot()
