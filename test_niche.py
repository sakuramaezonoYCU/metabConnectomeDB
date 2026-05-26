import scanpy as sc
import pandas as pd
import numpy as np
import os

# Load the file
h5ad_path = 'output/cancer_breast-cancer_breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad'
adata = sc.read_h5ad(h5ad_path)

PRIMARY_TISSUES = ['breast', 'mammary gland']
adata.obs['site'] = np.where(adata.obs['tissue_general'].isin(PRIMARY_TISSUES), 'Primary', 'Metastasis')

meta_sites = [t for t in adata[adata.obs['site'] == 'Metastasis'].obs['tissue_general'].unique() if pd.notnull(t)]
print(f"Meta sites: {meta_sites}")

site_significant_genes = {}
target_genes = ["GAPDH", "HK1"] # Dummy targets
for site in meta_sites:
    adata_site = adata[(adata.obs['tissue_general'] == site) | (adata.obs['site'] == 'Primary')].copy()
    adata_site.obs['comparison_group'] = adata_site.obs.apply(lambda x: site if x['site'] == 'Metastasis' else 'Primary', axis=1)
    print(f"Site: {site}, Cells: {adata_site.shape[0]}")
    sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')
    result = adata_site.uns['rank_genes_groups']
    print(f"Computed DE for {site}")
