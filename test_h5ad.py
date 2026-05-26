import scanpy as sc
import sys

adata = sc.read_h5ad("output/breast_results/breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad")
print("Columns:", adata.obs.columns.tolist())
if 'tissue' in adata.obs.columns:
    print("Unique tissues:", adata.obs['tissue'].unique().tolist())
if 'tissue_general' in adata.obs.columns:
    print("Unique tissue_general:", adata.obs['tissue_general'].unique().tolist())
