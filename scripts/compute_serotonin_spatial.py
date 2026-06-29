"""
compute_serotonin_spatial.py
============================
Computes spatial proximity between HTR7+ TAMs and TPH1+ Tumor cells
using the kNN graph from the h5ad files.

Replaces the previous hardcoded fake data with rigorous permutation-based
connectivity testing on the actual scRNA-seq graph.
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from pan_cancer_config import get_h5ad_path, CANCER_CAP, ANALYSIS_SUFFIX
from serotonin_config import HTR7_TAM_SIGNATURE

CANCERS = list(CANCER_CAP.keys())

def get_htr7_tam_mask(adata):
    """
    Identify HTR7+ Tumor-Associated Macrophages (TAMs).
    Requires expression of HTR7 > 0 and cell_type matching macrophage/TAM patterns.
    """
    macrophage_mask = adata.obs['cell_type'].astype(str).str.contains(
        'macrophage|mononuclear phagocyte', case=False, na=False
    )
    
    if 'HTR7' not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=bool)
        
    htr7_expr = adata[:, 'HTR7'].X.toarray().flatten() if hasattr(adata[:, 'HTR7'].X, 'toarray') else np.array(adata[:, 'HTR7'].X).flatten()
    htr7_mask = htr7_expr > 0
    
    return macrophage_mask & htr7_mask

def get_tph1_tumor_mask(adata):
    """
    Identify TPH1+ Tumor cells (Serotonin producers).
    """
    tumor_mask = adata.obs['cell_type'].astype(str).str.contains(
        'epithelial|tumor|malignant', case=False, na=False
    )
    
    if 'TPH1' not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=bool)
        
    tph1_expr = adata[:, 'TPH1'].X.toarray().flatten() if hasattr(adata[:, 'TPH1'].X, 'toarray') else np.array(adata[:, 'TPH1'].X).flatten()
    tph1_mask = tph1_expr > 0
    
    return tumor_mask & tph1_mask

def compute_graph_proximity(adata, mask_a, mask_b, n_perms=1000, seed=42, h5ad_path=None):
    """
    Compute real proximity score between two cell populations using kNN graph.
    If adata.obsp['connectivities'] is missing, compute it using sc.pp.neighbors 
    and save the updated adata back to h5ad_path to prevent recomputation.
    
    Returns:
        observed_score: float - mean connectivity between populations
        p_value: float - empirical p-value from permutation test
        z_score: float - standardized score vs null distribution
        null_distribution: np.array - all permutation scores (for plotting)
        
    CRITICAL: This function COMPUTES values. It NEVER returns hardcoded numbers.
    """
    # Ensure graph exists
    if 'connectivities' not in adata.obsp:
        print("  - 'connectivities' not found in obsp. Computing kNN graph on-the-fly...")
        if 'X_pca' not in adata.obsm:
            print("  - Computing PCA first...")
            sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        if h5ad_path:
            print(f"  - Saving updated h5ad with connectivities to {h5ad_path}...")
            # We don't want to accidentally corrupt the file if the save fails, 
            # but per user instruction, we must save it back.
            try:
                adata.write_h5ad(h5ad_path)
            except Exception as e:
                print(f"  - Warning: Failed to save updated h5ad: {e}")
                raise
                
    adj = adata.obsp['connectivities']
    
    idx_a = np.where(mask_a)[0]
    idx_b = np.where(mask_b)[0]
    
    if len(idx_a) == 0 or len(idx_b) == 0:
        return 0.0, 1.0, 0.0, np.zeros(n_perms)
        
    # Extract submatrix between population A and population B
    sub_adj = adj[idx_a, :][:, idx_b]
    observed_score = sub_adj.mean()
    
    # Permutation testing
    rng = np.random.default_rng(seed)
    n_cells = adata.n_obs
    null_scores = np.zeros(n_perms)
    
    len_a = len(idx_a)
    len_b = len(idx_b)
    
    print(f"  - Running {n_perms} permutations for significance testing...")
    for i in range(n_perms):
        # Randomly select indices of the same size
        rand_a = rng.choice(n_cells, len_a, replace=False)
        rand_b = rng.choice(n_cells, len_b, replace=False)
        
        rand_sub = adj[rand_a, :][:, rand_b]
        null_scores[i] = rand_sub.mean()
        
    # Calculate p-value
    p_value = np.mean(null_scores >= observed_score)
    
    # Calculate Z-score
    null_mean = np.mean(null_scores)
    null_std = np.std(null_scores)
    
    if null_std > 0:
        z_score = (observed_score - null_mean) / null_std
    else:
        z_score = 0.0
        
    return observed_score, p_value, z_score, null_scores

def run_spatial_analysis(cancer='ovarian'):
    """
    Runs the proximity analysis for a specific cancer.
    Loads h5ad, finds HTR7+ TAMs and TPH1+ Tumor cells, computes proximity,
    and saves results.
    """
    print(f"[{cancer.upper()}] Starting Serotonin Spatial Proximity Analysis")
    h5ad_path = get_h5ad_path(cancer)
    
    if not os.path.exists(h5ad_path):
        print(f"  -> File not found: {h5ad_path}")
        return
        
    print("  -> Loading AnnData...")
    adata = sc.read_h5ad(h5ad_path)
    
    # 1. Get masks
    print("  -> Identifying cell populations...")
    htr7_tam_mask = get_htr7_tam_mask(adata)
    tph1_tumor_mask = get_tph1_tumor_mask(adata)
    
    count_htr7 = htr7_tam_mask.sum()
    count_tph1 = tph1_tumor_mask.sum()
    print(f"  -> Found {count_htr7} HTR7+ TAMs and {count_tph1} TPH1+ Tumor cells.")
    
    out_dir = os.path.join(BASE_DIR, "output", "serotonin_axis_spatial_mapping")
    os.makedirs(out_dir, exist_ok=True)
    
    if count_htr7 < 5 or count_tph1 < 5:
        print("  -> Insufficient cells for robust proximity analysis. Skipping calculation.")
        summary = pd.DataFrame([{
            "Cancer": cancer,
            "HTR7_TAM_Count": count_htr7,
            "TPH1_Tumor_Count": count_tph1,
            "Observed_Proximity": 0.0,
            "P_Value": 1.0,
            "Z_Score": 0.0,
            "Note": "Insufficient cell counts"
        }])
        summary.to_csv(os.path.join(out_dir, f"serotonin_proximity_results_{cancer}{ANALYSIS_SUFFIX}.csv"), index=False)
        return
        
    # 2. Compute proximity
    print("  -> Computing graph-based proximity...")
    obs_score, p_val, z_score, null_dist = compute_graph_proximity(
        adata, htr7_tam_mask, tph1_tumor_mask, n_perms=1000, seed=42, h5ad_path=h5ad_path
    )
    
    print(f"  -> Result: Score = {obs_score:.6f}, P-value = {p_val:.4f}, Z-score = {z_score:.2f}")
    
    # 3. Save results
    summary = pd.DataFrame([{
        "Cancer": cancer,
        "HTR7_TAM_Count": count_htr7,
        "TPH1_Tumor_Count": count_tph1,
        "Observed_Proximity": obs_score,
        "P_Value": p_val,
        "Z_Score": z_score,
        "Note": "Computed via kNN permutation"
    }])
    
    out_csv = os.path.join(out_dir, f"serotonin_proximity_results_{cancer}{ANALYSIS_SUFFIX}.csv")
    summary.to_csv(out_csv, index=False)
    
    # Save null distribution for the notebook plotting
    np.save(os.path.join(out_dir, f"serotonin_null_dist_{cancer}{ANALYSIS_SUFFIX}.npy"), null_dist)
    
    print(f"  -> Results saved to {out_csv}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", type=str, default="ovarian", help="Cancer type to analyze")
    args = parser.parse_args()
    
    if args.cancer.lower() == "all":
        for c in CANCERS:
            run_spatial_analysis(c)
    else:
        run_spatial_analysis(args.cancer)
