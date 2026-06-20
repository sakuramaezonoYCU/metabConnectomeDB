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

import sys
if ".." not in sys.path:
    sys.path.append("..")
from pan_cancer_config import CANCERS_TO_RUN as CANCERS
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
                
    max_cancers = len(CANCERS)
    genes_strict = [g for g, c in gene_counts.items() if c == max_cancers]
    genes_broad = [g for g, c in gene_counts.items() if c >= max_cancers - 1]
    
    if len(genes_strict) == 0:
        print(f"Fail-safe triggered: 0 strictly conserved genes across all {max_cancers} cancers. Falling back to {max_cancers - 1} cancers.")
        genes_strict = genes_broad
    
    print(f"Found {len(genes_strict)} strictly conserved genes (after potential fallback).")
    print(f"Found {len(genes_broad)} broadly conserved (>={max_cancers - 1}) genes.")
    return genes_strict, genes_broad

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
            
    df = pd.DataFrame(results, columns=['Gene', 'Drug', 'Database'])
    if not df.empty:
        df = df.drop_duplicates(subset=['Gene', 'Drug'])
    return df

if __name__ == '__main__':
    genes_strict, genes_broad = get_conserved_genes()
    
    print("Querying DGIdb for the strictly conserved pan-cancer genes...")
    df_strict = query_dgidb_graphql(genes_strict)
    path_strict = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv')
    df_strict.to_csv(path_strict, index=False)
    print(f"Saved {path_strict} with {len(df_strict)} interactions.")
    
    print("Querying DGIdb for the broadly conserved (>=4) cancer genes...")
    df_broad = query_dgidb_graphql(genes_broad)
    path_broad = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv')
    df_broad.to_csv(path_broad, index=False)
    print(f"Saved {path_broad} with {len(df_broad)} interactions.")
