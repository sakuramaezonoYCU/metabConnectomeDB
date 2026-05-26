import scanpy as sc
import numpy as np
import pandas as pd

h5ad_path = "output/breast_results/breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad"
adata = sc.read_h5ad(h5ad_path)

PRIMARY_TISSUES = ['breast', 'mammary gland']
adata.obs['site'] = np.where(adata.obs['tissue_general'].isin(PRIMARY_TISSUES), 'Primary', 'Metastasis')

print("Site counts:")
print(adata.obs['site'].value_counts())

try:
    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
    print("DE Success")
except Exception as e:
    print(f"DE Failed: {e}")
