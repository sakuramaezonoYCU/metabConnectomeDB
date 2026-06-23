import os
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path: sys.path.append(BASE_DIR)
from pan_cancer_config import ANALYSIS_SUFFIX
try:
    from lifelines import CoxPHFitter
except ImportError:
    print("Please install lifelines (pip install lifelines) to run the permutation test.")
    exit(1)

def compute_permutation_null(signature_csv, n_permutations=100):
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_input = os.path.join(BASE_DIR, "input", "TCGA")
    out_dir = os.path.join(BASE_DIR, "output", f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}", "tcga_validation")
    os.makedirs(out_dir, exist_ok=True)
    
    from pan_cancer_config import TCGA_MAPPING
    cancer_codes = []
    for cohorts in TCGA_MAPPING.values():
        cancer_codes.extend(cohorts)
    cancer_codes = list(set(cancer_codes))
    if not cancer_codes:
        raise RuntimeError("CRITICAL ERROR: TCGA_MAPPING from pan_cancer_config is empty or missing. Cannot validate TCGA signatures without defined cohorts.")
    
    # Determine signature size
    if not os.path.exists(signature_csv):
        print(f"Missing {signature_csv}. Cannot determine signature size.")
        return
    
    sig_df = pd.read_csv(signature_csv)
    signature_size = len(sig_df)
    print(f"Dynamically determined signature size: {signature_size} from {signature_csv}")
    
    # 1. We need the list of all available genes to sample from.
    # Dynamically grab the first cohort to scan for the gene universe
    first_cohort = cancer_codes[0]
    first_exp = os.path.join(dir_input, f"TCGA-{first_cohort}.star_fpkm.tsv.gz")
    if not os.path.exists(first_exp):
        raise RuntimeError(f"Missing {first_exp}. Cannot determine gene universe.")
        
    print("Scanning transcriptome for gene universe...")
    genes_df = pd.read_csv(first_exp, sep='\t', usecols=[0])
    genes_df = genes_df[genes_df.iloc[:, 0].str.startswith('ENSG')]
    gene_col = genes_df.columns[0]
    all_genes = genes_df[gene_col].dropna().unique().tolist()
    print(f"Found {len(all_genes)} unique genes.")
    
    np.random.seed(42)
    random_sets = [np.random.choice(all_genes, size=signature_size, replace=False).tolist() for _ in range(n_permutations)]
    all_needed_genes = set([g for s in random_sets for g in s])
    
    all_results = []
    
    for cancer in cancer_codes:
        print(f"\nRunning {n_permutations} permutations for TCGA-{cancer}...")
        exp_file = os.path.join(dir_input, f"TCGA-{cancer}.star_fpkm.tsv.gz")
        surv_file_gdc = os.path.join(dir_input, f"TCGA-{cancer}.survival.tsv.gz")
        surv_file_gdc_unzipped = os.path.join(dir_input, f"TCGA-{cancer}.survival.tsv")
        surv_file_tcga = os.path.join(dir_input, f"{cancer}_survival.txt")
        
        if not os.path.exists(exp_file):
            continue
            
        surv_file = None
        for sf in [surv_file_gdc, surv_file_gdc_unzipped, surv_file_tcga]:
            if os.path.exists(sf):
                surv_file = sf
                break
                
        if surv_file is None:
            continue
            
        # Load Expression
        try:
            exp_df = pd.read_csv(exp_file, sep='\t')
            gene_col = exp_df.columns[0]
            # Filter only needed genes
            exp_df = exp_df[exp_df[gene_col].isin(all_needed_genes)]
            exp_df = exp_df.set_index(gene_col).T
            exp_df.index.name = 'sample'
            exp_df = exp_df.reset_index()
        except Exception as e:
            print(f"  Error reading {exp_file}: {e}")
            continue
            
        # Load Survival
        try:
            surv_df = pd.read_csv(surv_file, sep='\t')
            col_map = {c: c.lower() for c in surv_df.columns}
            surv_df = surv_df.rename(columns=col_map)
            sample_col = 'sample' if 'sample' in surv_df.columns else ('id' if 'id' in surv_df.columns else surv_df.columns[0])
            surv_df = surv_df.rename(columns={sample_col: 'sample'})
        except Exception as e:
            continue
            
        # Merge
        merged_base = pd.merge(exp_df, surv_df, on='sample', how='inner')
        if merged_base.empty:
            continue
            
        time_cols = [c for c in merged_base.columns if 'time' in c or 'days_to' in c or 'os.time' in c]
        event_cols = [c for c in merged_base.columns if 'event' in c or 'status' in c or c == 'os']
        
        if not time_cols or not event_cols:
            continue
            
        time_col = time_cols[0]
        event_col = event_cols[0]
        
        merged_base[time_col] = pd.to_numeric(merged_base[time_col], errors='coerce')
        merged_base[event_col] = pd.to_numeric(merged_base[event_col], errors='coerce')
        merged_base = merged_base.dropna(subset=[time_col, event_col])
        
        # Run iterations
        for i, random_sig in enumerate(tqdm(random_sets)):
            available_genes = [g for g in random_sig if g in merged_base.columns]
            if not available_genes:
                continue
                
            merged_base['Score'] = merged_base[available_genes].mean(axis=1)
            med = merged_base['Score'].median()
            merged_base['Risk_Binary'] = (merged_base['Score'] >= med).astype(int)
            
            cph = CoxPHFitter()
            fit_df = merged_base[[time_col, event_col, 'Risk_Binary']]
            cph.fit(fit_df, duration_col=time_col, event_col=event_col)
            hr = cph.hazard_ratios_['Risk_Binary']
            p_val = cph.summary['p']['Risk_Binary']
            
            all_results.append({
                'TCGA_Cohort': cancer,
                'Permutation': i,
                'Hazard_Ratio': hr,
                'P_Value': p_val
            })
    if all_results:
        results_df = pd.DataFrame(all_results)
        out_file = os.path.join(out_dir, "null_distribution_metrics.csv")
        results_df.to_csv(out_file, index=False)
        print(f"\nSaved null distribution of {len(results_df)} iterations to {out_file}")
    else:
        print("No null distribution metrics were calculated.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--signature_csv', required=True, help="Path to signature CSV")
    args = parser.parse_args()
    compute_permutation_null(args.signature_csv)
