import pandas as pd
import os
import json
import requests
import matplotlib.pyplot as plt
import seaborn as sns
from nr1d2_config import OUTPUT_DIR, CANCERS, ENRICHR_LIBRARIES, NR1D2_RESULTS_DIR

def get_pan_cancer_genes():
    """Finds metabolic genes consistently upregulated in metastasis across all 5 cancers."""
    cancer_genes = []
    
    for cancer in CANCERS:
        res_file = os.path.join(OUTPUT_DIR, f"{cancer}_results", f"primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv")
        if os.path.exists(res_file):
            df = pd.read_csv(res_file)
            up_genes = set(df[df['Significance'] == 'Up in Metastasis']['names'].tolist())
            cancer_genes.append(up_genes)
        else:
            print(f"Warning: {res_file} not found.")
            
    if not cancer_genes:
        return []
        
    pan_cancer = set.intersection(*cancer_genes)
    print(f"Found {len(pan_cancer)} pan-cancer upregulated metabolic genes.")
    return list(pan_cancer)

def run_enrichment_analysis():
    pan_genes = get_pan_cancer_genes()
    if not pan_genes:
        print("No pan-cancer genes found.")
        return None
        
    print(f"Found {len(pan_genes)} pan-cancer upregulated metabolic genes.")
    
    # Upload list once
    add_list_url = 'https://maayanlab.cloud/Enrichr/addList'
    payload = {
        'list': (None, '\n'.join(pan_genes)),
        'description': (None, 'Pan-Cancer Metabolic Targets')
    }
    response = requests.post(add_list_url, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    
    data = json.loads(response.text)
    user_list_id = data['userListId']
    
    all_results = []
    
    for lib in ENRICHR_LIBRARIES:
        print(f"Querying Enrichr library: {lib}...")
        enrich_url = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={lib}'
        res = requests.get(enrich_url)
        if not res.ok:
            print(f"Error fetching {lib}")
            continue
            
        results = json.loads(res.text)[lib]
        # Enrichr format: [Rank, Term name, P-value, Z-score, Combined score, Overlapping genes, Adjusted p-value, Old p-value, Old adjusted p-value]
        for row in results:
            term = row[1]
            pval = row[2]
            adj_pval = row[6]
            combined_score = row[4]
            overlap = row[5]
            
            # Enrichr term for TFs often contains the TF name first
            tf_name = term.split('_')[0].split(' ')[0].upper()
            
            all_results.append({
                'Library': lib,
                'Term': term,
                'Transcription Factor': tf_name,
                'P-value': pval,
                'Adjusted P-value': adj_pval,
                'Combined Score': combined_score,
                'Overlapping Genes': ','.join(overlap)
            })
            
    df_enrich = pd.DataFrame(all_results)
    
    # Save raw results
    csv_path = os.path.join(NR1D2_RESULTS_DIR, 'tf_enrichment_results.csv')
    df_enrich.to_csv(csv_path, index=False)
    print(f"Saved TF enrichment results to {csv_path}")
    
    return df_enrich

def plot_nr1d2_enrichment(df_enrich):
    if df_enrich is None or df_enrich.empty:
        return
        
    # Get top 15 TFs by combined score across all libraries
    top_tfs = df_enrich.sort_values('Combined Score', ascending=False).drop_duplicates(subset=['Transcription Factor']).head(15)
    
    plt.figure(figsize=(10, 8))
    sns.barplot(
        data=top_tfs,
        y='Transcription Factor',
        x='Combined Score',
        palette=['crimson' if tf == 'NR1D2' or tf == 'REV-ERBB' else 'lightgray' for tf in top_tfs['Transcription Factor']]
    )
    plt.title('Top Transcription Factors Regulating Pan-Cancer Metastatic Metabolism\n(Enrichr Combined Score)', pad=20)
    plt.xlabel('Enrichment Score (Combined Score)')
    plt.ylabel('Transcription Factor')
    
    # Add highlighting for NR1D2
    plt.tight_layout()
    plot_path = os.path.join(NR1D2_RESULTS_DIR, 'top_tfs_barplot.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {plot_path}")

if __name__ == "__main__":
    df = run_enrichment_analysis()
    plot_nr1d2_enrichment(df)
