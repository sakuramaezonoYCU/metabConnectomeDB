import os
import numpy as np
import scanpy as sc
import anndata as ad

print("Starting to prepare balanced lung test cache...")

# Paths
output_dir = 'output/lung_results'
lung_path = os.path.join(output_dir, 'lung_all_whole_transcriptome_2025-11-08.h5ad')
meta_path = os.path.join(output_dir, 'brain_adrenal-gland_skeletal-system_100k_whole_transcriptome_2025-11-08.h5ad')

# Expected output path for TISSUE_FILTER = ['lung', 'brain', 'bone', 'liver', 'adrenal gland']
# when CAP = 5000 is used
# tissue_slug = "lung_brain_bone_liver_adrenal-gland"
# filename: lung_brain_bone_liver_adrenal-gland_5k_whole_transcriptome_2025-11-08.h5ad
target_path_1 = os.path.join(output_dir, 'lung_brain_bone_liver_adrenal-gland_5k_whole_transcriptome_2025-11-08.h5ad')
# Also save to the "lung_lung-brain..." path just in case
target_path_2 = os.path.join(output_dir, 'lung_lung-brain-bone-liver-adrenal-gland_5k_whole_transcriptome_2025-11-08.h5ad')

# Set random seed
np.random.seed(42)

print("Loading primary lung data...")
adata_lung = sc.read_h5ad(lung_path)
print(f"Primary lung data shape: {adata_lung.shape}")

print("Loading metastatic data...")
adata_meta = sc.read_h5ad(meta_path)
print(f"Metastatic data shape: {adata_meta.shape}")

# Sample 2500 primary lung cells
print("Sampling 2500 primary lung cells...")
idx_lung = np.random.choice(adata_lung.obs_names, size=2500, replace=False)
adata_lung_sub = adata_lung[idx_lung].copy()

# Sample 2500 metastatic cells
print("Sampling 2500 metastatic cells (1250 brain, 750 adrenal, 500 skeletal system)...")
brain_cells = adata_meta.obs_names[adata_meta.obs['tissue_general'] == 'brain']
adrenal_cells = adata_meta.obs_names[adata_meta.obs['tissue_general'] == 'adrenal gland']
skeletal_cells = adata_meta.obs_names[adata_meta.obs['tissue_general'] == 'skeletal system']

idx_brain = np.random.choice(brain_cells, size=1250, replace=False)
idx_adrenal = np.random.choice(adrenal_cells, size=750, replace=False)
idx_skeletal = np.random.choice(skeletal_cells, size=500, replace=False)

idx_meta = np.concatenate([idx_brain, idx_adrenal, idx_skeletal])
adata_meta_sub = adata_meta[idx_meta].copy()

print("Concatenating subsets...")
# Concatenate observations
adata_combined = ad.concat([adata_lung_sub, adata_meta_sub], join='inner')

# Copy over var metadata
adata_combined.var = adata_lung.var.copy()

print(f"Combined shape: {adata_combined.shape}")
print("Tissue distribution in combined dataset:")
print(adata_combined.obs['tissue_general'].value_counts())

# Save combined datasets
print(f"Saving to {target_path_1}...")
adata_combined.write_h5ad(target_path_1)

print(f"Saving to {target_path_2}...")
adata_combined.write_h5ad(target_path_2)

print("SUCCESS: Balanced lung test cache is ready!")
