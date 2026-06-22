"""
Purpose: Identifies the total and unique metabolic target candidates implicated in cell-cell communication (CCC) potential across different tissue microenvironments.
"""
import sys
import os
import pandas as pd
import scanpy as sc
from pan_cancer_config import get_h5ad_path, ANALYSIS_SUFFIX

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output', 'pan_cancer_meta_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

from pan_cancer_config import CANCER_CAP, CANCER_PRIMARY_TISSUE
cancers_config = []
for cancer in CANCER_CAP.keys():
    tissue = CANCER_PRIMARY_TISSUE.get(cancer, cancer)
    cancers_config.append((cancer, tissue))

results = []

for prefix, tissue in cancers_config:
    h5ad_path = get_h5ad_path(prefix)
    adata = sc.read_h5ad(h5ad_path, backed='r')
    
    # Use strict masks by default to prevent normal epithelial contamination
    from pan_cancer_config import STRICT_MASK_CANCERS
    if prefix in STRICT_MASK_CANCERS:
        tumor_mask = adata.obs['cell_type'] == 'malignant cell'
    else:
        tumor_mask = adata.obs['cell_type'].str.contains('malignant|tumor|epithelial|cancer', case=False, na=False)
        
    if isinstance(tissue, list):
        total_pri = adata[adata.obs['tissue_general'].isin(tissue)]
        total_met = adata[~adata.obs['tissue_general'].isin(tissue)]
        adata_pri = adata[tumor_mask & (adata.obs['tissue_general'].isin(tissue))]
        adata_met = adata[tumor_mask & (~adata.obs['tissue_general'].isin(tissue))]
    else:
        total_pri = adata[adata.obs['tissue_general'] == tissue]
        total_met = adata[adata.obs['tissue_general'] != tissue]
        adata_pri = adata[tumor_mask & (adata.obs['tissue_general'] == tissue)]
        adata_met = adata[tumor_mask & (adata.obs['tissue_general'] != tissue)]
        
    pri_count = len(adata_pri)
    met_count = len(adata_met)
    t_pri = len(total_pri)
    t_met = len(total_met)
    
    print(f"{prefix.capitalize()}: Primary={pri_count}/{t_pri}, Metastatic={met_count}/{t_met}")
    results.append({
        'Dataset': prefix.capitalize(),
        'Total Primary TME Cells': t_pri,
        'Primary Malignant Cells': pri_count,
        'Total Metastatic TME Cells': t_met,
        'Metastatic Malignant Cells': met_count
    })

df = pd.DataFrame(results)
csv_path = os.path.join(OUTPUT_DIR, f'cell_type_counts{ANALYSIS_SUFFIX}.csv')
df.to_csv(csv_path, index=False)
print(f"Saved counts to {csv_path}")
