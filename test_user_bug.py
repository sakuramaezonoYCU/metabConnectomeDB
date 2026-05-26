import scanpy as sc
import numpy as np
import pandas as pd

h5ad_path = 'output/breast_cancer/breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad'
PRIMARY_TISSUES = ['breast']

adata = sc.read_h5ad(h5ad_path)

print(f"Original site dtype (if exists): {adata.obs['site'].dtype if 'site' in adata.obs else 'does not exist'}")

adata.obs['site'] = np.where(adata.obs['tissue_general'].isin(PRIMARY_TISSUES), 'Primary', 'Metastasis')

print(f"New site dtype: {adata.obs['site'].dtype}")
print(f"New site value counts:\n{adata.obs['site'].value_counts()}")

try:
    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
    print("Main DE PASSED!")
except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"Main DE FAILED: {e}")
