import os
import sys
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import pearsonr
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

import cellxgene_census
import serotonin_config as cfg

# Set basic configuration
sc.settings.verbosity = 3
output_dir = os.path.join(BASE_DIR, 'output', 'serotonin_axis_spatial_mapping')
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir

H5AD_CACHE_PATH = os.path.join(BASE_DIR, 'output', 'ovarian_full_uncapped_cellxgene.h5ad')

def fetch_ovarian_census():
    """Fetches the full ovarian cancer dataset from cellxgene census."""
    if os.path.exists(H5AD_CACHE_PATH):
        print(f"Loading cached dataset from {H5AD_CACHE_PATH}...")
        adata = sc.read_h5ad(H5AD_CACHE_PATH)
        if 'feature_name' in adata.var:
            adata.var_names = adata.var['feature_name'].values
            adata.var_names_make_unique()
        return adata
    
    # Build the list of genes we actually need for scoring
    immune_genes = ['CD3E', 'CD4', 'CD8A', 'NCAM1', 'NCR1', 'CD14', 'FCGR3A', 'PTPRC']
    for gene_list in cfg.SIGNATURES_TO_SCORE.values():
        immune_genes.extend(gene_list)
        
    genes_to_fetch = list(set(immune_genes))
    
    print("Downloading full ovarian cancer dataset from CellxGene Census. This may take a while...")
    # Fetch data
    with cellxgene_census.open_soma(census_version='2025-11-08') as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            measurement_name="RNA",
            obs_value_filter="disease in ['ovarian cancer', 'malignant ovarian serous tumor'] and is_primary_data == True and tissue_general in ['ovary', 'abdomen', 'omentum', 'uterus']",
            var_value_filter=f"feature_name in {genes_to_fetch}",
            column_names={"obs": ["tissue_general", "cell_type", "disease"]}
        )
    
    # Preprocessing
    print("Preprocessing data (normalization, log1p)...")
    adata.var_names = adata.var['feature_name'].values
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Save to cache
    print(f"Saving to {H5AD_CACHE_PATH}...")
    adata.write_h5ad(H5AD_CACHE_PATH)
    return adata

def safely_score_genes(adata, gene_list, score_name):
    """Scores genes, strictly checking if they exist in the matrix to prevent random/hallucinated data."""
    valid_genes = [g for g in gene_list if g in adata.var_names]
    missing = [g for g in gene_list if g not in adata.var_names]
    if missing:
        print(f"  WARNING: Missing genes for {score_name}: {missing}")
        
    if not valid_genes:
        print(f"  ERROR: No valid genes found for {score_name}. Returning NaN.")
        adata.obs[score_name] = np.nan
        return False
        
    print(f"  Scoring {score_name} using: {valid_genes}")
    import scipy.sparse
    if scipy.sparse.issparse(adata.X):
        adata.obs[score_name] = adata[:, valid_genes].X.mean(axis=1).A1
    else:
        adata.obs[score_name] = adata[:, valid_genes].X.mean(axis=1)
    return True

def main():
    print("1. Fetching Ovarian scRNA-seq Data")
    adata = fetch_ovarian_census()
    print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes.")

    print("\n2. Scoring Immune Signatures")
    for score_name, gene_list in cfg.SIGNATURES_TO_SCORE.items():
        safely_score_genes(adata, gene_list, score_name)
    
    print("\n3. Defining Niches & Tissues")
    # Tissues
    adata.obs['Is_Metastatic'] = adata.obs['tissue_general'].isin(['omentum', 'abdomen', 'uterus'])
    adata.obs['Is_Primary'] = adata.obs['tissue_general'].isin(['ovary'])
    
    # Cell categories (CellxGene standard cell_type ontologies)
    cell_types = adata.obs['cell_type'].str.lower().fillna('')
    adata.obs['Is_Macrophage'] = cell_types.str.contains('macrophage') | cell_types.str.contains('monocyte')
    adata.obs['Is_CAF'] = cell_types.str.contains('fibroblast')
    adata.obs['Is_T_NK_Cell'] = cell_types.str.contains('t cell') | cell_types.str.contains('natural killer')
    
    print("\n4. Aggregating Results")
    
    # To answer the prompt: we check if HTR7+ TAMs and CAFs are enriched in Metastatic vs Primary
    # and if they correlate with Exhaustion Targets
    
    results = []
    for tissue_class, mask in [('Primary', adata.obs['Is_Primary']), ('Metastatic', adata.obs['Is_Metastatic'])]:
        sub = adata[mask]
        if sub.n_obs == 0:
            continue
            
        # Macrophages
        macs = sub[sub.obs['Is_Macrophage']]
        mean_htr7_tam = macs.obs['HTR7_TAM_Score'].mean() if macs.n_obs > 0 else np.nan
        mean_suppressive_tam = macs.obs['Suppressive_Target_Score'].mean() if macs.n_obs > 0 else np.nan
        
        # CAFs
        cafs = sub[sub.obs['Is_CAF']]
        mean_htr7_caf = cafs.obs['HTR7_TAM_Score'].mean() if cafs.n_obs > 0 else np.nan
        mean_suppressive_caf = cafs.obs['Suppressive_Target_Score'].mean() if cafs.n_obs > 0 else np.nan
        
        # T/NK Cells
        tnk = sub[sub.obs['Is_T_NK_Cell']]
        mean_exhaustion = tnk.obs['Exhaustion_Target_Score'].mean() if tnk.n_obs > 0 else np.nan
        mean_treg = tnk.obs['Treg_Target_Score'].mean() if tnk.n_obs > 0 else np.nan
        
        results.append({
            'Niche': tissue_class,
            'Total_Cells': sub.n_obs,
            'Macrophage_Count': macs.n_obs,
            'CAF_Count': cafs.n_obs,
            'T_NK_Count': tnk.n_obs,
            'Mean_HTR7_TAM_Score': mean_htr7_tam,
            'Mean_HTR7_CAF_Score': mean_htr7_caf,
            'Mean_Suppressive_TAM_Score': mean_suppressive_tam,
            'Mean_Suppressive_CAF_Score': mean_suppressive_caf,
            'Mean_Exhaustion_Score': mean_exhaustion,
            'Mean_Treg_Score': mean_treg
        })
        
    df_res = pd.DataFrame(results)
    print("\nPrimary vs Metastatic Enrichment:")
    print(df_res.to_string())
    
    out_csv = os.path.join(output_dir, 'primary_vs_metastatic_immune_evasion_summary.csv')
    df_res.to_csv(out_csv, index=False)
    print(f"Saved summary to {out_csv}")
    
    # Save the processed metadata (we don't need the whole matrix) for notebook plotting
    out_meta = os.path.join(output_dir, 'immune_evasion_metadata.csv')
    cols_to_save = ['tissue_general', 'cell_type', 'Is_Metastatic', 'Is_Primary', 
                    'Is_Macrophage', 'Is_CAF', 'Is_T_NK_Cell',
                    'HTR7_TAM_Score', 'Exhaustion_Target_Score', 
                    'Suppressive_Target_Score', 'Treg_Target_Score']
    adata.obs[cols_to_save].to_csv(out_meta)
    print(f"Saved full metadata to {out_meta}")

if __name__ == '__main__':
    main()
