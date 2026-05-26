import scanpy as sc
import pandas as pd
import numpy as np

h5ad_path = 'output/breast_cancer/breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad'
adata = sc.read_h5ad(h5ad_path)
adata = adata[np.random.choice(adata.n_obs, 500, replace=False)].copy()

# Deliberately mismatched tissue
PRIMARY_TISSUES = ['lung']
adata.obs['site'] = np.where(adata.obs['tissue_general'].isin(PRIMARY_TISSUES), 'Primary', 'Metastasis')

print("Cells by site:")
print(adata.obs['site'].value_counts())

# Test Main DE
try:
    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
    print("Main DE passed!")
except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"Main DE failed: {type(e).__name__}: {e}")

