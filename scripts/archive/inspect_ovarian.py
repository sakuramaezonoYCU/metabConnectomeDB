import sys
import scanpy as sc
from pan_cancer_config import get_h5ad_path

def inspect_ovarian():
    print("Loading Ovarian...")
    adata_ov = sc.read_h5ad(get_h5ad_path('ovarian'), backed='r')
    print("\nOvarian cell types:")
    print(adata_ov.obs['cell_type'].value_counts())
    print("\nOvarian tissue_general:")
    print(adata_ov.obs['tissue_general'].value_counts())
    print("\nOvarian malignant cell tissue distribution:")
    if 'malignant cell' in adata_ov.obs['cell_type'].values:
        print(adata_ov.obs[adata_ov.obs['cell_type'] == 'malignant cell']['tissue_general'].value_counts())
    else:
        print("NO MALIGNANT CELLS FOUND IN OVARIAN!")

if __name__ == "__main__":
    inspect_ovarian()
