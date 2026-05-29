import pandas as pd
import os
import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_de_csv_path

import requests

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output')
DRUGGABILITY_DIR = os.path.join(OUTPUT_DIR, 'druggability')
os.makedirs(DRUGGABILITY_DIR, exist_ok=True)

CANCERS = ['breast', 'colorectal', 'lung', 'melanoma', 'ovarian']
DGIDB_GQL = "https://dgidb.org/api/graphql"

def get_conserved_genes():
    cancer_genes = []
    gene_counts = {}
    
    for cancer in CANCERS:
        res_file = get_de_csv_path(cancer)
        if os.path.exists(res_file):
            df = pd.read_csv(res_file)
            up_genes = df[df['Significance'] == 'Up in Metastasis']['names'].tolist()
            cancer_genes.append(set(up_genes))
            for g in up_genes:
                gene_counts[g] = gene_counts.get(g, 0) + 1
                
    genes_23 = [g for g, c in gene_counts.items() if c == 5]
    genes_181 = [g for g, c in gene_counts.items() if c >= 4]
    
    print(f"Found {len(genes_23)} all-5 conserved genes.")
    print(f"Found {len(genes_181)} >=4 conserved genes.")
    return genes_23, genes_181

def query_dgidb_graphql(genes):
    print(f"[DGIdb] Querying {len(genes)} genes...")
    results = []
    
    chunk_size = 20
    for i in range(0, len(genes), chunk_size):
        chunk = genes[i:i+chunk_size]
        genes_str = '", "'.join(chunk)
        q = f"""
        query {{
          genes(names: ["{genes_str}"]) {{
            nodes {{
              name
              interactions {{
                drug {{
                  name
                }}
              }}
            }}
          }}
        }}
        """
        try:
            r = requests.post(DGIDB_GQL, json={"query": q})
            data = r.json()
            nodes = data.get('data', {}).get('genes', {}).get('nodes', [])
            for node in nodes:
                gene_name = node.get('name')
                interactions = node.get('interactions', [])
                for interaction in interactions:
                    drug_name = interaction.get('drug', {}).get('name')
                    if drug_name:
                        results.append({
                            'Gene': gene_name,
                            'Drug': drug_name,
                            'Database': 'DGIdb'
                        })
        except Exception as e:
            print(f"Error querying {chunk}: {e}")
            
    df = pd.DataFrame(results)
    if not df.empty:
        df = df.drop_duplicates(subset=['Gene', 'Drug'])
    return df

if __name__ == '__main__':
    genes_23, genes_181 = get_conserved_genes()
    
    print("Querying DGIdb for the strictly conserved pan-cancer genes...")
    df_23 = query_dgidb_graphql(genes_23)
    path_23 = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv')
    df_23.to_csv(path_23, index=False)
    print(f"Saved {path_23} with {len(df_23)} interactions.")
    
    print("Querying DGIdb for the broadly conserved (>=4) cancer genes...")
    df_181 = query_dgidb_graphql(genes_181)
    path_181 = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv')
    df_181.to_csv(path_181, index=False)
    print(f"Saved {path_181} with {len(df_181)} interactions.")
