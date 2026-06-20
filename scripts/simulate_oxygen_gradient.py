import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
import os
import glob
import scanpy as sc
import pandas as pd
import numpy as np

def simulate_oxygen_gradient():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(base_dir, "output")
    out_dir = os.path.join(output_dir, f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}", "oxygen_gradient")
    os.makedirs(out_dir, exist_ok=True)
    
    from pan_cancer_config import CANCER_PRIMARY_TISSUE
    
    # Find all h5ad files
    h5ad_files = glob.glob(os.path.join(output_dir, "*_results", "*100k_whole_transcriptome_2025-11-08.h5ad"))
    
    if not h5ad_files:
        print("Error: No h5ad files found.")
        return
        
    for h5ad_path in h5ad_files:
        # Extract cancer name from path
        cancer_prefix = os.path.basename(os.path.dirname(h5ad_path)).replace('_results', '')
        print(f"\n--- Processing {cancer_prefix} ---")
        print(f"Loading {h5ad_path}...")
        
        try:
            adata = sc.read_h5ad(h5ad_path)
        except Exception as e:
            print(f"Error loading {h5ad_path}: {e}")
            continue
        
        # Define signatures
        hypoxia_sig = ["VEGFA", "SLC2A1", "BNIP3"]
        pan_cancer_file = os.path.join(output_dir, "pan_cancer_meta_results", f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv")
        
        if not os.path.exists(pan_cancer_file):
            print(f"Error: {pan_cancer_file} not found.")
            continue
            
        pan_cancer_genes = pd.read_csv(pan_cancer_file)['Strictly_Conserved_Gene'].tolist()
        
        # Filter for genes present in adata
        hypoxia_present = [g for g in hypoxia_sig if g in adata.var_names]
        pan_cancer_present = [g for g in pan_cancer_genes if g in adata.var_names]
        
        if len(hypoxia_present) == 0 or len(pan_cancer_present) == 0:
            print("Missing necessary genes for scoring. Skipping.")
            continue
            
        print(f"Scoring {len(hypoxia_present)} hypoxia genes and {len(pan_cancer_present)} pan-cancer genes...")
        
        # Score cells
        sc.tl.score_genes(adata, gene_list=hypoxia_present, score_name="hypoxia_score", use_raw=False)
        sc.tl.score_genes(adata, gene_list=pan_cancer_present, score_name="metastatic_score", use_raw=False)
        
        df = adata.obs[['tissue_general', 'cell_type', 'hypoxia_score', 'metastatic_score']].copy()
        
        primary_val = CANCER_PRIMARY_TISSUE.get(cancer_prefix, cancer_prefix)
        if isinstance(primary_val, list):
            primary_keyword = '|'.join(primary_val)
        else:
            primary_keyword = primary_val
        
        # Filter for primary tumor cells
        primary_tumor = df[(df['tissue_general'].str.contains(primary_keyword, case=False, na=False)) & 
                           (df['cell_type'].str.contains('malignant|tumor|epithelial', case=False, na=False))].copy()
        
        print(f"Found {len(primary_tumor)} primary tumor cells.")
        if len(primary_tumor) == 0:
            # Fallback to all cells in primary if specific cell type filtering fails
            primary_tumor = df[(df['tissue_general'].str.contains(primary_keyword, case=False, na=False))].copy()
            print(f"Fallback: Found {len(primary_tumor)} primary cells of any type.")
        
        if len(primary_tumor) > 0:
            # Calculate correlation
            correlation = primary_tumor['hypoxia_score'].corr(primary_tumor['metastatic_score'])
            print(f"Correlation between Hypoxia and Metastatic Score in Primary: {correlation:.3f}")
            
            # Save results
            out_file = os.path.join(out_dir, f"{cancer_prefix}_primary_oxygen_gradient_scores.csv")
            primary_tumor.to_csv(out_file)
            print(f"Saved {out_file}")

if __name__ == "__main__":
    simulate_oxygen_gradient()
