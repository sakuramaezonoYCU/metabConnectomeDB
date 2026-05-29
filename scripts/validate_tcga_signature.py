import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
import os
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter

def validate_tcga_signature():
    print("Validating strictly-conserved metabolic signature in TCGA cohorts...")
    
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_input = os.path.join(BASE_DIR, "input", "TCGA")
    out_dir = os.path.join(BASE_DIR, "output", f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}", "tcga_validation")
    os.makedirs(out_dir, exist_ok=True)
    
    # Gene signature
    sig_file = os.path.join(BASE_DIR, "output", "pan_cancer_meta_results", f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv")
    if not os.path.exists(sig_file):
        print(f"Error: {sig_file} not found. Please run earlier pipeline steps first.")
        return
        
    sig_df = pd.read_csv(sig_file)
    if 'Strictly_Conserved_Gene' in sig_df.columns:
        metab_genes = sig_df['Strictly_Conserved_Gene'].tolist()
    elif 'gene' in sig_df.columns:
        metab_genes = sig_df['gene'].tolist()
    else:
        metab_genes = sig_df.iloc[:, 0].tolist()
        
    # Map to Ensembl IDs via local HGNC database
    import json
    hgnc_path = os.path.join(BASE_DIR, "input", "hgnc_approved_genes.json")
    valid_ensembl_ids = set()
    
    if os.path.exists(hgnc_path):
        with open(hgnc_path, 'r') as f:
            hgnc_data = json.load(f)
            
        docs = hgnc_data.get('response', {}).get('docs', [])
        for doc in docs:
            if doc.get('symbol') in metab_genes and 'ensembl_gene_id' in doc:
                valid_ensembl_ids.add(doc['ensembl_gene_id'])
                
        print(f"Mapped {len(metab_genes)} symbols to {len(valid_ensembl_ids)} Ensembl IDs via local HGNC.")
    else:
        print(f"Error: {hgnc_path} not found. Cannot map gene symbols to Ensembl IDs.")
        return
        
    cancer_codes = ["BRCA", "COAD", "READ", "LUAD", "LUSC", "SKCM", "OV"]
    all_cancer_metrics = []
    
    for cancer in cancer_codes:
        print(f"Processing TCGA-{cancer}...")
        exp_file = os.path.join(dir_input, f"TCGA-{cancer}.star_fpkm.tsv.gz")
        surv_file_gdc = os.path.join(dir_input, f"TCGA-{cancer}.survival.tsv.gz")
        surv_file_gdc_unzipped = os.path.join(dir_input, f"TCGA-{cancer}.survival.tsv")
        surv_file_tcga = os.path.join(dir_input, f"{cancer}_survival.txt")
        surv_file_tcga_tsv = os.path.join(dir_input, f"{cancer}_survival.tsv")
        
        if not os.path.exists(exp_file):
            print(f"  Missing expression data for {cancer}")
            continue
            
        surv_file = None
        for sf in [surv_file_gdc, surv_file_gdc_unzipped, surv_file_tcga, surv_file_tcga_tsv]:
            if os.path.exists(sf):
                surv_file = sf
                break
                
        if surv_file is None:
            print(f"  Missing survival data for {cancer}")
            continue
            
        # 1. Load Expression Data
        try:
            exp_df = pd.read_csv(exp_file, sep='\t')
            gene_col = exp_df.columns[0]
            
            # Create a temporary column with stripped Ensembl IDs (removing .15 etc)
            stripped_ids = exp_df[gene_col].astype(str).str.split('.').str[0]
            
            # Filter only needed genes using our mapped Ensembl IDs
            mask = stripped_ids.isin(valid_ensembl_ids)
            exp_df = exp_df[mask]
            
            if exp_df.empty:
                print(f"  No signature genes found in {exp_file}. (ID mismatch?)")
                continue
                
            # Replace the Ensembl ID with the stripped version for cleaner matching/scoring if needed
            exp_df.loc[:, gene_col] = stripped_ids[mask]
            
            exp_df = exp_df.set_index(gene_col).T
            exp_df['signature_score'] = exp_df.sum(axis=1)
            exp_df.index.name = 'sample'
            exp_df = exp_df.reset_index()
            
        except Exception as e:
            print(f"  Error reading {exp_file}: {e}")
            continue
            
        # 2. Load Survival Data
        try:
            surv_df = pd.read_csv(surv_file, sep='\t')
            col_map = {c: c.lower() for c in surv_df.columns}
            surv_df = surv_df.rename(columns=col_map)
            
            sample_col = 'sample' if 'sample' in surv_df.columns else ('id' if 'id' in surv_df.columns else surv_df.columns[0])
            surv_df = surv_df.rename(columns={sample_col: 'sample'})
        except Exception as e:
            print(f"  Error reading {surv_file}: {e}")
            continue
            
        # 3. Merge and Compute
        merged = pd.merge(exp_df[['sample', 'signature_score']], surv_df, on='sample', how='inner')
        if merged.empty:
            print(f"  No overlapping samples for {cancer}.")
            continue
            
        time_cols = [c for c in merged.columns if 'time' in c or 'days_to' in c or 'os.time' in c]
        event_cols = [c for c in merged.columns if 'event' in c or 'status' in c or c == 'os']
        
        if not time_cols or not event_cols:
            print(f"  Could not find time/event columns in survival data for {cancer}.")
            continue
            
        time_col = time_cols[0]
        event_col = event_cols[0]
        
        merged[time_col] = pd.to_numeric(merged[time_col], errors='coerce')
        merged[event_col] = pd.to_numeric(merged[event_col], errors='coerce')
        merged = merged.dropna(subset=[time_col, event_col, 'signature_score'])
        
        if len(merged) < 20:
            print(f"  Not enough valid samples for {cancer}.")
            continue
            
        median_score = merged['signature_score'].median()
        merged['risk_group'] = (merged['signature_score'] >= median_score).astype(int)
        
        cph = CoxPHFitter()
        try:
            cox_df = merged[[time_col, event_col, 'risk_group']]
            cph.fit(cox_df, duration_col=time_col, event_col=event_col)
            hr = cph.hazard_ratios_['risk_group']
            p_val = cph.summary['p']['risk_group']
        except Exception as e:
            print(f"  CoxPH failed: {e}")
            hr = 1.0
            p_val = 1.0
            
        all_cancer_metrics.append({
            'TCGA_Cohort': cancer,
            'Hazard_Ratio': hr,
            'P_Value': p_val,
            'N_Samples': len(merged)
        })
        
    if all_cancer_metrics:
        results_df = pd.DataFrame(all_cancer_metrics)
        out_file = os.path.join(out_dir, "true_signature_metrics.csv")
        results_df.to_csv(out_file, index=False)
        print(f"\nSaved TCGA validation metrics to {out_file}")
        print(results_df)
    else:
        print("\nNo metrics were generated.")

if __name__ == "__main__":
    validate_tcga_signature()
