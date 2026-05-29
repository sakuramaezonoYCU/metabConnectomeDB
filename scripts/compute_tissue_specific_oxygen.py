import pandas as pd
import numpy as np
import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from oxygen_tension_config import GLYCOLYSIS_GENES, OXPHOS_GENES, BASE_DIR

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "oxygen_tension")
os.makedirs(OUTPUT_DIR, exist_ok=True)

H5AD_FILE = os.path.join(BASE_DIR, "output", "breast_results", "breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad")

def compute_tissue_specific_ratios():
    print(f"Loading {H5AD_FILE}...")
    adata = sc.read_h5ad(H5AD_FILE)
    
    print("Filtering for malignant cells...")
    adata_malig = adata[adata.obs['cell_type'] == 'malignant cell'].copy()
    
    # Check if we have multiple tissues
    tissues = adata_malig.obs['tissue_general'].unique().tolist()
    print(f"Tissues found: {tissues}")
    
    if 'breast' not in tissues:
        print("Error: 'breast' not found in tissues. Cannot compute DE against primary.")
        return
        
    print("Running Differential Expression (Wilcoxon) vs primary (breast)...")
    sc.tl.rank_genes_groups(adata_malig, groupby='tissue_general', reference='breast', method='wilcoxon')
    
    results = []
    
    met_tissues = [t for t in tissues if t != 'breast']
    
    for tissue in met_tissues:
        print(f"Processing {tissue}...")
        de_df = sc.get.rank_genes_groups_df(adata_malig, group=tissue)
        de_df['names'] = de_df['names'].str.upper()
        
        glyco_df = de_df[de_df['names'].isin(GLYCOLYSIS_GENES)]
        oxphos_df = de_df[de_df['names'].isin(OXPHOS_GENES)]
        
        mean_glyco_lfc = glyco_df['logfoldchanges'].mean() if not glyco_df.empty else 0
        mean_oxphos_lfc = oxphos_df['logfoldchanges'].mean() if not oxphos_df.empty else 0
        
        lfc_diff = mean_oxphos_lfc - mean_glyco_lfc
        oxphos_glyco_ratio = 2 ** lfc_diff
        
        results.append({
            "Tissue": tissue.capitalize(),
            "Mean_Glycolysis_LFC": mean_glyco_lfc,
            "Mean_OXPHOS_LFC": mean_oxphos_lfc,
            "OXPHOS_Glycolysis_Ratio": oxphos_glyco_ratio,
            "Glycolysis_Genes_Found": len(glyco_df),
            "OXPHOS_Genes_Found": len(oxphos_df)
        })
        
    df_res = pd.DataFrame(results)
    csv_path = os.path.join(OUTPUT_DIR, "tissue_specific_oxygen_ratios.csv")
    df_res.to_csv(csv_path, index=False)
    print(f"Saved {csv_path}")
    
    # Plotting
    plt.figure(figsize=(8, 6))
    sns.barplot(data=df_res, x="Tissue", y="OXPHOS_Glycolysis_Ratio", palette="viridis")
    plt.axhline(1, color='red', linestyle='--', alpha=0.7, label='Equal Ratio (1.0)')
    plt.yscale('log')
    plt.ylabel('OXPHOS / Glycolysis Ratio (Log Scale)')
    plt.title('Tissue-Specific Metabolic Adaptation in Breast Cancer Metastasis')
    plt.legend()
    
    plot_path = os.path.join(OUTPUT_DIR, "tissue_specific_correlation_5MetCan_100k.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {plot_path}")

if __name__ == "__main__":
    compute_tissue_specific_ratios()
