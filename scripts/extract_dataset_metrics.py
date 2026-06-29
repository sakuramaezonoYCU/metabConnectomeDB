"""
Purpose: Extracts exactly how many cells were evaluated for the TME and Malignant compartments across primary and metastatic sites.
"""
import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_h5ad_path
import os
import pandas as pd
import scanpy as sc

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')
os.makedirs(META_RESULTS_DIR, exist_ok=True)

# Configuration for cancer datasets.
# Format: (h5ad_prefix, tissue_general_filter, display_name)
# The display_name (3rd value) is the proper capitalized name used in plots and reports.
from pan_cancer_config import CANCER_CAP, CANCER_PRIMARY_TISSUE
cancers_config = []
for cancer in CANCER_CAP.keys():
    tissue = CANCER_PRIMARY_TISSUE.get(cancer, cancer)
    title_name = cancer.capitalize()
    cancers_config.append((cancer, tissue, title_name))

data = []
for prefix, tissue, title_name in cancers_config:
    h5ad_path = get_h5ad_path(prefix)
    try:
        adata = sc.read_h5ad(h5ad_path, backed='r')
    except Exception as e:
        print(f"Could not load {h5ad_path}: {e}")
        continue
        raise
    
    obs = adata.obs
    
    # Identify malignant cells vs TME cells
    tumor_mask = obs['cell_type'].str.contains('malignant|tumor|epithelial|cancer', case=False, na=False)
    
    if isinstance(tissue, list):
        primary_mask = obs['tissue_general'].isin(tissue)
    else:
        primary_mask = obs['tissue_general'] == tissue
        
    met_mask = ~primary_mask
    
    pri_tme = len(obs[primary_mask & ~tumor_mask])
    pri_mal = len(obs[primary_mask & tumor_mask])
    met_tme = len(obs[met_mask & ~tumor_mask])
    met_mal = len(obs[met_mask & tumor_mask])
    
    data.append({
        'Dataset': title_name,
        'Total Primary TME Cells': pri_tme,
        'Primary Malignant Cells': pri_mal,
        'Total Metastatic TME Cells': met_tme,
        'Metastatic Malignant Cells': met_mal
    })

df = pd.DataFrame(data)
out_path = os.path.join(META_RESULTS_DIR, f'cell_type_counts{ANALYSIS_SUFFIX}.csv')
df.to_csv(out_path, index=False)
print(f"Saved dataset overview counts to {out_path}")
