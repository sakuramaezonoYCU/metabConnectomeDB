import os
import json
import urllib.request
import pandas as pd

# Define paths
META_RESULTS_DIR = 'output/pan_cancer_meta_results'
OUTPUT_DIR = 'output'
ANALYSIS_SUFFIX = '_Br500k_Co100k_Lu500k_Me100k_Ov100k'

def fetch_uniprot_roles():
    conserved_csv = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    if not os.path.exists(conserved_csv):
        print(f"File not found: {conserved_csv}")
        return
        
    sig_df = pd.read_csv(conserved_csv)
    genes = sig_df['Strictly_Conserved_Gene'].tolist()
    
    results = []
    
    import requests
    import re
    
    print("Fetching biological roles from UniProt API...")
    for gene in genes:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606&fields=accession,id,gene_names,protein_name,cc_function&format=json"
        
        try:
            response = requests.get(url, headers={'User-Agent': 'python-requests/2.25.1'})
            response.raise_for_status()
            
            import gzip
            if response.content.startswith(b'\x1f\x8b'):
                data = json.loads(gzip.decompress(response.content).decode('utf-8'))
            else:
                data = response.json()
                
            function_text = ""
            if 'results' in data and len(data['results']) > 0:
                first_result = data['results'][0]
                
                # Search for the FUNCTION comment
                if 'comments' in first_result:
                    for comment in first_result['comments']:
                        if comment.get('commentType') == 'FUNCTION':
                            if 'texts' in comment and len(comment['texts']) > 0:
                                function_text = comment['texts'][0].get('value', '')
                                break
                                
            # Clean up the text a bit (remove PubMed citations to keep it concise in the table)
            function_text = re.sub(r'\(PubMed:\d+(?:, PubMed:\d+)*\)', '', function_text).strip()
            
            if not function_text:
                function_text = "No curated function found in UniProt."
                
            results.append({
                "Gene": gene,
                "UniProt_Biological_Role": function_text
            })
            print(f"Fetched: {gene}")
            
        except Exception as e:
            print(f"Error fetching {gene}: {e}")
            results.append({
                "Gene": gene,
                "UniProt_Biological_Role": "Error fetching from API."
            })
            
    # Save to CSV
    df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, 'uniprot_biological_roles.csv')
    df.to_csv(out_path, index=False)
    print(f"\nSaved UniProt biological roles to {out_path}")

if __name__ == "__main__":
    fetch_uniprot_roles()
