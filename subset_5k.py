import scanpy as sc
import os
import numpy as np

src_path = "output/lung_results/lung_all_whole_transcriptome_2025-11-08.h5ad"
dst_path = "output/lung_results/lung_brain_bone_liver_adrenal-gland_5k_whole_transcriptome_2025-11-08.h5ad"

print(f"Loading {src_path}...")
adata = sc.read_h5ad(src_path)
print(f"Original shape: {adata.shape}")

# Randomly subset 5000 cells
if adata.n_obs > 5000:
    indices = np.random.choice(adata.n_obs, 5000, replace=False)
    adata = adata[indices].copy()
    print(f"Subset shape: {adata.shape}")

print(f"Saving to {dst_path}...")
adata.write_h5ad(dst_path)
print("Done!")
