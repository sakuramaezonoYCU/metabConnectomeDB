import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import argparse

import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_de_csv_path, CANCERS_TO_RUN, TCGA_MAPPING

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')

def load_data(meta_results_dir):
    diff_file = os.path.join(meta_results_dir, 'differential_abundance_per_cancer.csv')
    match_file = os.path.join(meta_results_dir, 'metabolite_match_table.csv')
    
    if not os.path.exists(diff_file) or not os.path.exists(match_file):
        print(f"Skipping Cross-Cohort Comparison: Missing input files in {meta_results_dir} (Likely 0 matched metabolites).")
        import sys
        sys.exit(0)

    diff_df = pd.read_csv(diff_file)
    match_df = pd.read_csv(match_file)
    
    # We map MS metabolite to Target Gene list
    # Because one MS metabolite can link to multiple genes, and one gene to multiple metabolites,
    # let's create a flat mapping: Gene -> ms_metabolite
    gene_to_ms = []
    for _, row in match_df.iterrows():
        ms_metab = row['ms_metabolite']
        genes = str(row['linked_genes']).split(';')
        for g in genes:
            gene_to_ms.append({'Gene': g, 'ms_metabolite': ms_metab})
            
    gene_ms_df = pd.DataFrame(gene_to_ms).drop_duplicates()
    
    return diff_df, gene_ms_df

def process_cohort(cohort_name, sc_file, ms_cancer, diff_df, gene_ms_df):
    if not os.path.exists(sc_file):
        print(f"Skipping {cohort_name}: File not found ({sc_file})")
        return None
        
    sc_df = pd.read_csv(sc_file)
    # Filter sc_df to only Significant ones?
    # sc_df has columns: names, logfoldchanges, pvals_adj
    sc_df = sc_df.rename(columns={'names': 'Gene', 'logfoldchanges': 'scRNA_logFC'})
    
    # Get MS differential results for this cancer
    ms_cancer_df = diff_df[diff_df['Cancer'] == ms_cancer].copy()
    if ms_cancer_df.empty:
        print(f"No MS data for {ms_cancer}")
        return None
        
    ms_cancer_df = ms_cancer_df.rename(columns={'Metabolite': 'ms_metabolite', 'Log2FC': 'MS_log2FC'})
    
    # Merge
    merged = pd.merge(sc_df, gene_ms_df, on='Gene', how='inner')
    merged = pd.merge(merged, ms_cancer_df, on='ms_metabolite', how='inner')
    
    # Filter for genes with significant scRNA p-values if any, or just take top
    merged['Cohort'] = cohort_name
    return merged

def plot_cross_cohort(all_merged, meta_results_dir):
    if all_merged.empty:
        print("No intersecting data found.")
        return
        
    # Scatter plot: scRNA_logFC vs MS_log2FC
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=all_merged, 
        x='scRNA_logFC', 
        y='MS_log2FC', 
        hue='Cohort', 
        style='Cohort',
        s=150, 
        alpha=0.8
    )
    
    # Label top points
    for _, row in all_merged.iterrows():
        # Only label if highly changed in either
        if abs(row['scRNA_logFC']) > 1.5 or abs(row['MS_log2FC']) > 1.0:
            label = f"{row['Gene']}\n({row['ms_metabolite']})"
            plt.text(row['scRNA_logFC']+0.05, row['MS_log2FC'], label, fontsize=8)

    plt.axhline(0, color='gray', linestyle='--')
    plt.axvline(0, color='gray', linestyle='--')
    plt.title("Cross-Cohort Comparison: scRNA (Primary vs Met) vs Mass Spec (Tumor vs Normal)")
    plt.xlabel("scRNA-seq logFC (Metastasis vs Primary)")
    plt.ylabel("Mass Spec log2FC (Tumor vs Normal)")
    
    out_file = os.path.join(meta_results_dir, 'cross_cohort_comparison_scatter.png')
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save table
    all_merged.to_csv(os.path.join(meta_results_dir, 'cross_cohort_comparison_table.csv'), index=False)
    print(f"Cross cohort comparison complete. Saved scatter plot and table.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Mass Spec Cross Cohort Comparison")
    parser.add_argument('--signature-name', required=True, help="Name of the signature (e.g. Conserved)")
    args = parser.parse_args()
    
    print(f"Running Cross-Cohort Comparison for signature: {args.signature_name}...")
    meta_results_dir = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', args.signature_name)
    
    diff_df, gene_ms_df = load_data(meta_results_dir)
    
    dfs = []
    
    for cancer in CANCERS_TO_RUN:
        tcga_ids = TCGA_MAPPING.get(cancer)
        if not tcga_ids:
            print(f"Warning: No TCGA ID mapping for {cancer}. Skipping.")
            continue
        
        if isinstance(tcga_ids, str):
            tcga_ids = [tcga_ids]
            
        sc_file = get_de_csv_path(cancer)
        
        # We process each TCGA ID separately and then append
        for tcga_id in tcga_ids:
            display_name = f"{cancer.capitalize()} ({tcga_id})"
            merged = process_cohort(display_name, sc_file, tcga_id, diff_df, gene_ms_df)
            
            if merged is not None and not merged.empty:
                dfs.append(merged)
            
    if dfs:
        all_merged = pd.concat(dfs, ignore_index=True)
        plot_cross_cohort(all_merged, meta_results_dir)
    else:
        print("No data available to plot cross-cohort.")
