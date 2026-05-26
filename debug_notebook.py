import scanpy as sc
import pandas as pd
import numpy as np
import os

BASE_DIR = '../'
GLOBAL_OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
OUTPUT_DIR = os.path.join(GLOBAL_OUTPUT_DIR, 'breast_cancer')
h5ad_path = os.path.join(OUTPUT_DIR, 'breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad')
PRIMARY_TISSUES = ['breast']

print(f"Loading {h5ad_path}")
try:
    adata = sc.read_h5ad(h5ad_path)
    adata.obs['site'] = np.where(adata.obs['tissue_general'].isin(PRIMARY_TISSUES), 'Primary', 'Metastasis')

    print("Cells by site:")
    print(adata.obs['site'].value_counts())
    print("Unique tissues:")
    print(adata.obs['tissue_general'].value_counts())
except Exception as e:
    print(f"Failed: {e}")
