import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import fdrcorrection
import argparse
from sklearn.decomposition import PCA

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join(BASE_DIR, 'input')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
MASS_SPEC_FILE = os.path.join(INPUT_DIR, 'massSpecDataMetabolicData_7panCancer_PMID29396322.csv')
HMDB_ANNO_FILE = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')

def clean_name(name):
    """Normalize metabolite name for better exact matching"""
    n = str(name).lower().strip()
    # Remove common stereochemistry/isomers prefixes that might mismatch
    n = re.sub(r'^(l-|d-|dl-|alpha-|beta-|gamma-|n-|\(r\)-|\(s\)-)', '', n)
    # Remove asterisks and trailing spaces
    n = n.replace('*', '').strip()
    return n

# -----------------------------------------------------------------------------
# Module 1: Metabolite Matching
# -----------------------------------------------------------------------------
def match_metabolites(gene_list, meta_results_dir):
    print("Module 1: Loading data and matching metabolites...")
    ms_df = pd.read_csv(MASS_SPEC_FILE, index_col=0)
    # The first column 'X' contains metabolite names
    ms_metabolites = ms_df['X'].tolist()
    ms_clean_map = {clean_name(m): m for m in ms_metabolites}
    
    hmdb_df = pd.read_csv(HMDB_ANNO_FILE, low_memory=False)
    hmdb_filtered = hmdb_df[hmdb_df['Target'].isin(gene_list)].copy()
    
    hmdb_metab_map = {}
    gene_links = {}
    
    for _, row in hmdb_filtered.iterrows():
        orig_name = str(row['Metabolite_Name'])
        target = row['Target']
        
        # Helper to add to maps
        def add_name(n):
            cn = clean_name(n)
            hmdb_metab_map[cn] = orig_name
            if cn not in gene_links:
                gene_links[cn] = set()
            gene_links[cn].add(target)

        add_name(orig_name)
        if pd.notna(row.get('HMDB_Name')):
            add_name(row['HMDB_Name'])
        if pd.notna(row.get('Synonyms')):
            for syn in str(row['Synonyms']).split(';'):
                add_name(syn)
                
    matches = []
    for cn, ms_orig in ms_clean_map.items():
        if cn in hmdb_metab_map:
            genes = list(gene_links[cn])
            matches.append({
                'ms_metabolite': ms_orig,
                'hmdb_metabolite': hmdb_metab_map[cn],
                'linked_genes': ';'.join(genes),
                'match_type': 'Cleaned Exact/Synonym'
            })
            
    if not matches:
        match_df = pd.DataFrame(columns=['ms_metabolite', 'hmdb_metabolite', 'linked_genes', 'match_type'])
    else:
        match_df = pd.DataFrame(matches)
        
    match_df.to_csv(os.path.join(meta_results_dir, 'metabolite_match_table.csv'), index=False)
    print(f"Found {len(match_df)} matching metabolites linked to signature genes.")
    
    return ms_df, match_df, hmdb_filtered

# -----------------------------------------------------------------------------
# Module 2: Differential Abundance
# -----------------------------------------------------------------------------
def run_differential_abundance(ms_df, match_df, meta_results_dir):
    print("Module 2: Running Differential Abundance...")
    if match_df.empty:
        print("No matched metabolites. Skipping Differential Abundance.")
        pd.DataFrame(columns=['Cancer', 'Metabolite', 'Tumor_Samples', 'Normal_Samples', 'Log2FC', 'P_Value', 'FDR']).to_csv(os.path.join(meta_results_dir, 'differential_abundance_per_cancer.csv'), index=False)
        return pd.DataFrame(), pd.DataFrame(), {}
        
    matched_ms_names = match_df['ms_metabolite'].tolist()
    # Filter ms_df to only matched metabolites
    ms_sub = ms_df[ms_df['X'].isin(matched_ms_names)].set_index('X')
    
    cols = ms_sub.columns.tolist()
    
    # Parse cancer types
    cancer_samples = {}
    for c in cols:
        parts = c.split('.')
        if len(parts) >= 3:
            cancer = parts[0]
            stype = parts[-1]
            if cancer not in cancer_samples:
                cancer_samples[cancer] = {'Tumor': [], 'Normal': []}
            if stype in ['Tumor', 'Normal']:
                cancer_samples[cancer][stype].append(c)
                
    results = []
    for cancer, samps in cancer_samples.items():
        tumors = samps['Tumor']
        normals = samps['Normal']
        
        if len(tumors) < 5 or len(normals) < 5:
            continue
            
        tumor_data = ms_sub[tumors]
        normal_data = ms_sub[normals]
        
        for metab in ms_sub.index:
            t_vals = tumor_data.loc[metab].dropna().values
            n_vals = normal_data.loc[metab].dropna().values
            
            if len(t_vals) < 3 or len(n_vals) < 3:
                continue
                
            try:
                stat, pval = mannwhitneyu(t_vals, n_vals, alternative='two-sided')
                # Log2 fold change of medians
                t_med = np.median(t_vals)
                n_med = np.median(n_vals)
                # handle zeros
                t_med = t_med if t_med > 0 else 1e-6
                n_med = n_med if n_med > 0 else 1e-6
                lfc = np.log2(t_med / n_med)
                
                results.append({
                    'Cancer': cancer,
                    'Metabolite': metab,
                    'Tumor_Samples': len(t_vals),
                    'Normal_Samples': len(n_vals),
                    'Log2FC': lfc,
                    'P_Value': pval
                })
            except Exception:
                pass
                
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        # FDR correction per cancer
        res_df['FDR'] = 1.0
        for cancer in res_df['Cancer'].unique():
            mask = res_df['Cancer'] == cancer
            pvals = res_df.loc[mask, 'P_Value']
            _, fdr = fdrcorrection(pvals)
            res_df.loc[mask, 'FDR'] = fdr
            
        res_df.to_csv(os.path.join(meta_results_dir, 'differential_abundance_per_cancer.csv'), index=False)
        
        # Volcano plots
        for cancer in res_df['Cancer'].unique():
            sub = res_df[res_df['Cancer'] == cancer]
            plt.figure(figsize=(8, 6))
            sns.scatterplot(x='Log2FC', y=-np.log10(sub['FDR']+1e-10), data=sub, 
                            hue=(sub['FDR'] < 0.05), palette={True: 'red', False: 'grey'}, s=100)
            plt.title(f'{cancer} Tumor vs Normal Volcano Plot')
            plt.axhline(-np.log10(0.05), ls='--', color='black')
            plt.axvline(0, ls='-', color='black', alpha=0.5)
            
            # Label significant ones
            sig = sub[sub['FDR'] < 0.05]
            for _, row in sig.iterrows():
                plt.text(row['Log2FC'], -np.log10(row['FDR']+1e-10), row['Metabolite'], fontsize=8)
                
            plt.tight_layout()
            plt.savefig(os.path.join(meta_results_dir, f'volcano_{cancer}.png'), dpi=300)
            plt.close()
    return res_df, ms_sub, cancer_samples

# -----------------------------------------------------------------------------
# Module 3: Co-abundance Correlation
# -----------------------------------------------------------------------------
def run_coabundance(ms_sub, meta_results_dir):
    print("Module 3: Co-abundance Correlation Heatmap...")
    if len(ms_sub) < 2:
        print("Not enough metabolites for co-abundance clustering. Skipping.")
        return
        
    # Fill NaN with min value or row median for correlation
    # We will compute spearmanr pairwise, ignoring NaNs
    corr_matrix = pd.DataFrame(index=ms_sub.index, columns=ms_sub.index, dtype=float)
    
    for i, m1 in enumerate(ms_sub.index):
        for j, m2 in enumerate(ms_sub.index):
            if i <= j:
                v1 = ms_sub.loc[m1]
                v2 = ms_sub.loc[m2]
                mask = ~(v1.isna() | v2.isna())
                if mask.sum() > 10:
                    r, _ = spearmanr(v1[mask], v2[mask])
                    corr_matrix.iloc[i, j] = r
                    corr_matrix.iloc[j, i] = r
                else:
                    corr_matrix.iloc[i, j] = 0
                    corr_matrix.iloc[j, i] = 0
                    
    corr_matrix = corr_matrix.fillna(0)
    plt.figure(figsize=(10, 8))
    sns.clustermap(corr_matrix, cmap='coolwarm', center=0, figsize=(10, 10))
    plt.title("Metabolite Co-abundance (Spearman)")
    plt.savefig(os.path.join(meta_results_dir, 'metabolite_coabundance_heatmap.png'), dpi=300)
    plt.close()

# -----------------------------------------------------------------------------
# Module 4: PCA Metabolic Signature
# -----------------------------------------------------------------------------
def run_pca_signature(ms_sub, cancer_samples, meta_results_dir):
    print("Module 4: PCA Metabolic Signature Score...")
    if ms_sub.empty or len(ms_sub) == 0:
        print("Not enough metabolites. Skipping PCA.")
        return
        
    # Impute missing values with row median for PCA
    ms_imputed = ms_sub.apply(lambda row: row.fillna(row.median()), axis=1)
    
    pca = PCA(n_components=1)
    # Standardize rows (metabolites)
    X = ms_imputed.T.values
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0) + 1e-6
    X_scaled = (X - X_mean) / X_std
    
    pca_scores = pca.fit_transform(X_scaled).flatten()
    
    # Check if PC1 correlates positively with the majority of metabolites, if not, flip sign
    loadings = pca.components_[0]
    if np.sum(loadings < 0) > np.sum(loadings > 0):
        pca_scores = -pca_scores
        
    score_df = pd.DataFrame({
        'Sample': ms_imputed.columns,
        'Signature_Score': pca_scores
    })
    
    # Add metadata
    score_df['Cancer'] = 'Unknown'
    score_df['Type'] = 'Unknown'
    for c, samps in cancer_samples.items():
        for t in ['Tumor', 'Normal']:
            score_df.loc[score_df['Sample'].isin(samps[t]), 'Cancer'] = c
            score_df.loc[score_df['Sample'].isin(samps[t]), 'Type'] = t
            
    score_df = score_df[score_df['Type'] != 'Unknown']
    score_df.to_csv(os.path.join(meta_results_dir, 'metabolite_signature_pca_scores.csv'), index=False)
    
    # Plot boxplot
    plt.figure(figsize=(14, 6))
    sns.boxplot(x='Cancer', y='Signature_Score', hue='Type', data=score_df)
    plt.title('Metabolic Signature Score by Cancer and Sample Type')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(meta_results_dir, 'metabolite_pca_boxplot.png'), dpi=300)
    plt.close()

# -----------------------------------------------------------------------------
# Module 5: Bipartite Network
# -----------------------------------------------------------------------------
def run_bipartite_network(hmdb_filtered, match_df, meta_results_dir):
    print("Module 5: Gene-Metabolite Bipartite Network...")
    if hmdb_filtered.empty:
        print("No genes to build network. Skipping.")
        pd.DataFrame(columns=['Gene', 'Metabolite', 'Is_Detected_MassSpec', 'Super_Class']).to_csv(os.path.join(meta_results_dir, 'gene_metabolite_network_edges.csv'), index=False)
        return
        
    matched_hmdb = set(match_df['hmdb_metabolite']) if not match_df.empty else set()
    
    G = nx.Graph()
    edges = []
    
    for _, row in hmdb_filtered.iterrows():
        metab = row['Metabolite_Name']
        target = row['Target']
        is_matched = metab in matched_hmdb
        sc = row.get('Super_Class', 'Unknown')
        
        G.add_node(metab, bipartite=0, matched=is_matched, super_class=sc)
        G.add_node(target, bipartite=1, matched=True, super_class='Gene')
        G.add_edge(metab, target)
        
        edges.append({
            'Gene': target,
            'Metabolite': metab,
            'Is_Detected_MassSpec': is_matched,
            'Super_Class': sc
        })
        
    edges_df = pd.DataFrame(edges)
    edges_df.to_csv(os.path.join(meta_results_dir, 'gene_metabolite_network_edges.csv'), index=False)
    
    # Plot
    plt.figure(figsize=(16, 12))
    pos = nx.spring_layout(G, k=0.5, iterations=50)
    
    # Nodes
    metab_matched = [n for n, d in G.nodes(data=True) if d.get('bipartite')==0 and d.get('matched')]
    metab_unmatched = [n for n, d in G.nodes(data=True) if d.get('bipartite')==0 and not d.get('matched')]
    genes = [n for n, d in G.nodes(data=True) if d.get('bipartite')==1]
    
    nx.draw_networkx_nodes(G, pos, nodelist=metab_matched, node_color='green', node_size=300, alpha=0.8, label='Detected Metab')
    nx.draw_networkx_nodes(G, pos, nodelist=metab_unmatched, node_color='lightgray', node_size=100, alpha=0.5, label='Undetected Metab')
    nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color='lightcoral', node_size=500, alpha=0.9, label='Target Gene')
    
    nx.draw_networkx_edges(G, pos, alpha=0.3)
    
    # Labels for genes and matched metabs
    labels = {n: n for n in genes + metab_matched}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
    
    plt.legend()
    plt.title("Gene vs Metabolite Network (Green = Found in Mass Spec)")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(meta_results_dir, 'gene_metabolite_network.png'), dpi=300)
    plt.close()

# -----------------------------------------------------------------------------
# Module 6: Pan-Cancer Heatmap
# -----------------------------------------------------------------------------
def run_pan_cancer_heatmap(ms_sub, cancer_samples, meta_results_dir):
    print("Module 6: Pan-Cancer Metabolite Heatmap...")
    if ms_sub.empty or not cancer_samples:
        print("No data. Skipping Heatmap.")
        return
        
    tumor_medians = {}
    
    for c, samps in cancer_samples.items():
        tumors = samps['Tumor']
        if tumors:
            tumor_medians[c] = ms_sub[tumors].median(axis=1)
            
    if tumor_medians:
        df_med = pd.DataFrame(tumor_medians)
        # Z-score rows
        df_z = df_med.apply(lambda x: (x - x.mean()) / (x.std() + 1e-6), axis=1)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(df_z, cmap='viridis', center=0)
        plt.title('Pan-Cancer Median Tumor Metabolite Abundance (Z-score)')
        plt.tight_layout()
        plt.savefig(os.path.join(meta_results_dir, 'pan_cancer_metabolite_heatmap.png'), dpi=300)
        plt.close()

# -----------------------------------------------------------------------------
# Module 7: Per-Gene Profile
# -----------------------------------------------------------------------------
def run_per_gene_profile(edges_df, res_df, meta_results_dir):
    print("Module 7: Per-Gene Metabolite Profile...")
    if edges_df.empty:
        print("No edges. Skipping Profile.")
        pd.DataFrame(columns=['Gene', 'Total_Linked_Metabolites', 'Detected_in_MassSpec', 'Detected_Metabolites']).to_csv(os.path.join(meta_results_dir, 'per_gene_metabolite_profile.csv'), index=False)
        return
        
    profile = []
    
    for gene in edges_df['Gene'].unique():
        sub = edges_df[edges_df['Gene'] == gene]
        all_m = sub['Metabolite'].tolist()
        det_m = sub[sub['Is_Detected_MassSpec']]['Metabolite'].tolist()
        
        # Summarize differential abundance
        sig_cancers = set()
        if res_df is not None and not res_df.empty:
            for m in det_m:
                # We need to map HMDB name to MS name, or just search loosely
                # simpler: just list the detection
                pass
                
        profile.append({
            'Gene': gene,
            'Total_Linked_Metabolites': len(all_m),
            'Detected_in_MassSpec': len(det_m),
            'Detected_Metabolites': ';'.join(det_m)
        })
        
    df_prof = pd.DataFrame(profile)
    df_prof.to_csv(os.path.join(meta_results_dir, 'per_gene_metabolite_profile.csv'), index=False)

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Mass Spec Metabolomics Integration")
    parser.add_argument('--signature_csv', required=True, help="Path to the signature CSV file containing a 'Gene' or 'Strictly_Conserved_Gene' column")
    args = parser.parse_args()

    if not os.path.exists(args.signature_csv):
        raise FileNotFoundError(f"CRITICAL ERROR: Signature CSV {args.signature_csv} does not exist.")

    df_sig = pd.read_csv(args.signature_csv)
    if 'Gene' in df_sig.columns:
        genes = df_sig['Gene'].dropna().unique().tolist()
    elif 'Target' in df_sig.columns:
        genes = df_sig['Target'].dropna().unique().tolist()
    elif 'Strictly_Conserved_Gene' in df_sig.columns:
        genes = df_sig['Strictly_Conserved_Gene'].dropna().unique().tolist()
    else:
        raise ValueError(f"CRITICAL ERROR: Could not find 'Gene', 'Target', or 'Strictly_Conserved_Gene' column in {args.signature_csv}")

    sig_name = os.path.basename(args.signature_csv).replace('.csv', '')

    meta_results_dir = os.path.join(OUTPUT_DIR, 'massspec_metabolomics', sig_name)
    os.makedirs(meta_results_dir, exist_ok=True)
    
    ms_df, match_df, hmdb_filtered = match_metabolites(genes, meta_results_dir)
    res_df, ms_sub, cancer_samples = run_differential_abundance(ms_df, match_df, meta_results_dir)
    run_coabundance(ms_sub, meta_results_dir)
    run_pca_signature(ms_sub, cancer_samples, meta_results_dir)
    run_bipartite_network(hmdb_filtered, match_df, meta_results_dir)
    run_pan_cancer_heatmap(ms_sub, cancer_samples, meta_results_dir)
    
    edges_df = pd.read_csv(os.path.join(meta_results_dir, 'gene_metabolite_network_edges.csv'))
    run_per_gene_profile(edges_df, res_df, meta_results_dir)
    print(f"Analysis complete. Check {meta_results_dir} for results.")
