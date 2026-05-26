import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.AnnData(np.random.randn(10, 10))
adata.obs['site'] = ['Primary']*5 + ['Metastasis']*5
try:
    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
    print("DE passed")
except Exception as e:
    print(f"DE failed: {e}")

adata.obs['site'] = adata.obs['site'].astype('category')
try:
    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
    print("DE with category passed")
except Exception as e:
    print(f"DE with category failed: {e}")
