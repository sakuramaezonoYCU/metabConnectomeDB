"""
compute_serotonin_htr7_axis.py
==============================
Computes signature scores for the Serotonin-HTR7-TAM-IP4-HR repair axis 
across all cell types in the h5ad files.
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from pan_cancer_config import get_h5ad_path, CANCER_CAP, ANALYSIS_SUFFIX
from serotonin_config import (
    HTR7_TAM_SIGNATURE, SEROTONIN_SYNTHESIS, SEROTONIN_TRANSPORT, 
    SEROTONIN_DEGRADATION, SEROTONIN_RECEPTORS, INOSITOL_PATHWAY, 
    HR_REPAIR_GENES, EV_BIOGENESIS, EV_CARGO_SORTING, SEROTONIN_DEP_SECRETION,
    get_all_genes
)

CANCERS = list(CANCER_CAP.keys())

def check_gene_presence(adata, cancer):
    """
    Checks which genes from the signatures are actually in the dataset.
    Returns a DataFrame.
    """
    results = []
    
    modules = {
        'HTR7_TAM': HTR7_TAM_SIGNATURE,
        'Serotonin_Synthesis': SEROTONIN_SYNTHESIS,
        'Serotonin_Transport': SEROTONIN_TRANSPORT,
        'Serotonin_Degradation': SEROTONIN_DEGRADATION,
        'Serotonin_Receptors': SEROTONIN_RECEPTORS,
        'Inositol_Pathway': INOSITOL_PATHWAY,
        'HR_Repair': HR_REPAIR_GENES,
        'EV_Biogenesis': EV_BIOGENESIS,
        'EV_Cargo': EV_CARGO_SORTING,
        'EV_Secretion': SEROTONIN_DEP_SECRETION
    }
    
    for mod_name, genes in modules.items():
        for g in genes:
            present = g in adata.var_names
            results.append({
                'cancer': cancer,
                'gene': g,
                'gene_set': mod_name,
                'present_in_h5ad': present
            })
            
    return pd.DataFrame(results)

def calculate_module_score(adata, gene_list, score_name):
    """
    Calculates a single cell module score using sc.tl.score_genes.
    Returns the array of scores.
    """
    valid_genes = [g for g in gene_list if g in adata.var_names]
    if len(valid_genes) == 0:
        return np.zeros(adata.n_obs)
        
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=score_name)
    score = adata.obs[score_name].values
    
    # Clean up the object
    del adata.obs[score_name]
    return score

def score_htr7_tam_signature(adata, cancer):
    """
    Scores the HTR7+ TAM signature for every cell_type x tissue.
    """
    # Create the score column
    adata.obs['tam_module_score'] = calculate_module_score(adata, HTR7_TAM_SIGNATURE, 'temp_tam')
    
    # Has HTR7?
    if 'HTR7' in adata.var_names:
        htr7_expr = adata[:, 'HTR7'].X.toarray().flatten() if hasattr(adata[:, 'HTR7'].X, 'toarray') else np.array(adata[:, 'HTR7'].X).flatten()
    else:
        htr7_expr = np.zeros(adata.n_obs)
        
    adata.obs['htr7_expr'] = htr7_expr
    adata.obs['is_htr7_pos'] = htr7_expr > 0
    
    results = []
    
    # Group by cell_type and tissue_general
    if 'tissue_general' not in adata.obs.columns:
        adata.obs['tissue_general'] = 'Unknown'
        
    for (ctype, tissue), group in adata.obs.groupby(['cell_type', 'tissue_general'], observed=False):
        n_cells = len(group)
        if n_cells < 10:
            continue
            
        htr7_frac = group['is_htr7_pos'].mean()
        htr7_mean = group['htr7_expr'].mean()
        tam_score = group['tam_module_score'].mean()
        
        results.append({
            'cancer': cancer,
            'cell_type': str(ctype),
            'tissue': str(tissue),
            'n_cells': n_cells,
            'HTR7_frac_positive': htr7_frac,
            'HTR7_mean_expr': htr7_mean,
            'TAM_module_score': tam_score
        })
        
    # Clean up
    del adata.obs['tam_module_score'], adata.obs['htr7_expr'], adata.obs['is_htr7_pos']
    return pd.DataFrame(results)

def score_serotonin_axis(adata, cancer):
    """
    Full serotonin axis expression: synthesis, transport, degradation, receptors.
    """
    genes_to_track = SEROTONIN_SYNTHESIS + SEROTONIN_TRANSPORT + SEROTONIN_DEGRADATION + SEROTONIN_RECEPTORS
    valid_genes = [g for g in genes_to_track if g in adata.var_names]
    
    results = []
    for (ctype, tissue), group_indices in adata.obs.groupby(['cell_type', 'tissue_general'], observed=False).groups.items():
        if len(group_indices) < 10:
            continue
            
        sub_adata = adata[group_indices, valid_genes]
        expr_matrix = sub_adata.X.toarray() if hasattr(sub_adata.X, 'toarray') else np.array(sub_adata.X)
        
        mean_expr = expr_matrix.mean(axis=0)
        frac_pos = (expr_matrix > 0).mean(axis=0)
        
        for i, g in enumerate(valid_genes):
            if g in SEROTONIN_SYNTHESIS: category = 'Synthesis'
            elif g in SEROTONIN_TRANSPORT: category = 'Transport'
            elif g in SEROTONIN_DEGRADATION: category = 'Degradation'
            else: category = 'Receptor'
            
            results.append({
                'cancer': cancer,
                'cell_type': str(ctype),
                'tissue': str(tissue),
                'gene': g,
                'mean_expr': mean_expr[i],
                'frac_positive': frac_pos[i],
                'gene_category': category
            })
            
    return pd.DataFrame(results)

def score_inositol_hr_repair(adata, cancer):
    """
    Scores Inositol pathway and HR repair pathway, computes cell-level correlation
    for each cell_type x tissue.
    """
    adata.obs['inositol_score'] = calculate_module_score(adata, INOSITOL_PATHWAY, 'temp_inos')
    adata.obs['hr_score'] = calculate_module_score(adata, HR_REPAIR_GENES, 'temp_hr')
    
    results = []
    for (ctype, tissue), group in adata.obs.groupby(['cell_type', 'tissue_general'], observed=False):
        if len(group) < 10:
            continue
            
        in_scores = group['inositol_score'].values
        hr_scores = group['hr_score'].values
        
        mean_in = np.mean(in_scores)
        mean_hr = np.mean(hr_scores)
        
        # Avoid constant array warnings
        if np.std(in_scores) > 0 and np.std(hr_scores) > 0:
            r, p = stats.pearsonr(in_scores, hr_scores)
        else:
            r, p = 0.0, 1.0
            
        results.append({
            'cancer': cancer,
            'cell_type': str(ctype),
            'tissue': str(tissue),
            'inositol_score': mean_in,
            'hr_score': mean_hr,
            'pearson_r': r,
            'p_value': p
        })
        
    del adata.obs['inositol_score'], adata.obs['hr_score']
    return pd.DataFrame(results)

def score_ev_machinery(adata, cancer):
    """
    Scores EV Biogenesis, Cargo, and Serotonin-dependent secretion.
    """
    adata.obs['ev_biogenesis'] = calculate_module_score(adata, EV_BIOGENESIS, 'temp_biogen')
    adata.obs['ev_cargo'] = calculate_module_score(adata, EV_CARGO_SORTING, 'temp_cargo')
    adata.obs['serotonin_sec'] = calculate_module_score(adata, SEROTONIN_DEP_SECRETION, 'temp_sec')
    
    results = []
    for (ctype, tissue), group in adata.obs.groupby(['cell_type', 'tissue_general'], observed=False):
        if len(group) < 10:
            continue
            
        results.append({
            'cancer': cancer,
            'cell_type': str(ctype),
            'tissue': str(tissue),
            'ev_biogenesis_score': group['ev_biogenesis'].mean(),
            'ev_cargo_score': group['ev_cargo'].mean(),
            'serotonin_dep_score': group['serotonin_sec'].mean()
        })
        
    del adata.obs['ev_biogenesis'], adata.obs['ev_cargo'], adata.obs['serotonin_sec']
    return pd.DataFrame(results)

def run_axis_computation(cancer='ovarian'):
    h5ad_path = get_h5ad_path(cancer)
    if not os.path.exists(h5ad_path):
        print(f"File not found: {h5ad_path}")
        return
        
    print(f"[{cancer.upper()}] Loading data for Serotonin Axis Computation...")
    adata = sc.read_h5ad(h5ad_path)
    
    out_dir = os.path.join(BASE_DIR, "output", "serotonin_axis_spatial_mapping")
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Gene Presence
    print("  -> Checking gene presence...")
    df_pres = check_gene_presence(adata, cancer)
    df_pres.to_csv(os.path.join(out_dir, f"gene_presence_{cancer}{ANALYSIS_SUFFIX}.csv"), index=False)
    
    # 2. HTR7+ TAM Scoring
    print("  -> Scoring HTR7+ TAM signature...")
    df_tam = score_htr7_tam_signature(adata, cancer)
    df_tam.to_csv(os.path.join(out_dir, f"htr7_tam_scoring_{cancer}{ANALYSIS_SUFFIX}.csv"), index=False)
    
    # 3. Full Serotonin Axis
    print("  -> Scoring full serotonin axis (receptors/transport)...")
    df_sero = score_serotonin_axis(adata, cancer)
    df_sero.to_csv(os.path.join(out_dir, f"serotonin_full_axis_{cancer}{ANALYSIS_SUFFIX}.csv"), index=False)
    
    # 4. Inositol & HR Repair
    print("  -> Scoring Inositol and HR Repair correlation...")
    df_hr = score_inositol_hr_repair(adata, cancer)
    df_hr.to_csv(os.path.join(out_dir, f"inositol_hr_correlation_{cancer}{ANALYSIS_SUFFIX}.csv"), index=False)
    
    # 5. EV Machinery
    print("  -> Scoring EV Machinery...")
    df_ev = score_ev_machinery(adata, cancer)
    df_ev.to_csv(os.path.join(out_dir, f"ev_machinery_{cancer}{ANALYSIS_SUFFIX}.csv"), index=False)
    
    print(f"[{cancer.upper()}] All module scoring complete.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer", type=str, default="ovarian", help="Cancer type to analyze")
    args = parser.parse_args()
    
    if args.cancer.lower() == "all":
        for c in CANCERS:
            run_axis_computation(c)
    else:
        run_axis_computation(args.cancer)
