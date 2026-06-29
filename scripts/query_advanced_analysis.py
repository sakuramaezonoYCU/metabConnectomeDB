import requests
import pandas as pd
from IPython.display import display
from pan_cancer_config import OPENTARGETS_API_URL

import os
import json
import hashlib
from datetime import datetime

def _get_api_cache_dir():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cache_dir = os.path.join(root, "input", "api_cache")
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir

def query_string_ppi(genes):
    """
    Query STRING Database for Protein-Protein Interactions among the target genes.
    This helps establish if the 'axis' has physical or functional evidence of interaction.
    """
    print(f"[STRING PPI] Querying network for: {genes}")
    cache_dir = _get_api_cache_dir()
    gene_hash = hashlib.md5("_".join(sorted(genes)).encode()).hexdigest()
    cache_file = os.path.join(cache_dir, f"string_ppi_{gene_hash}.csv")
    failures_file = os.path.join(cache_dir, "api_failures.log")

    if os.path.exists(cache_file):
        print(f"[STRING PPI] Loading cached network from version-controlled cache: {cache_file}")
        df = pd.read_csv(cache_file)
        if not df.empty:
            df['Source'] = 'Cache'
        return df

    url = "https://string-db.org/api/json/network"
    
    params = {
        "identifiers": "%0D".join(genes),
        "species": 9606, # Human
        "caller_identity": "metabConnectomeDB"
    }
    
    try:
        r = requests.get(url, params=params)
        r.raise_for_status()
        data = r.json()
        
        if not data:
            import warnings
            warnings.warn("[STRING PPI] No interactions found among these genes in STRING.")
            df = pd.DataFrame(columns=[])
            df.to_csv(cache_file, index=False)
            return df
            
        results = []
        for interaction in data:
            results.append({
                'Gene_A': interaction.get('preferredName_A'),
                'Gene_B': interaction.get('preferredName_B'),
                'Combined_Score': interaction.get('score'),
                'Experimental_Score': interaction.get('escore'),
                'Database_Score': interaction.get('dscore')
            })
            
        df = pd.DataFrame(results)
        df.to_csv(cache_file, index=False)
        print(f"[STRING PPI] Found {len(df)} interaction(s). Cached to version-controlled file: {cache_file}.")
        if not df.empty:
            df['Source'] = 'API'
        return df
    except Exception as e:
        msg = f"{datetime.now().isoformat()} - [STRING PPI] Failed to query STRING: {e}"
        print(f"{msg}. Logged to {failures_file}.")
        with open(failures_file, "a") as f:
            f.write(msg + "\n")
        return pd.DataFrame(columns=[])
        raise

def query_tractability(genes):
    """
    Query Open Targets for Target Tractability.
    Even if no drugs currently exist, tractability predicts if the target 
    is chemically druggable (small molecule) or targetable by antibodies.
    """
    print(f"[Tractability] Querying Open Targets for tractability of {genes}...")
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) metabConnectomeDB'}
    
    cache_dir = _get_api_cache_dir()
    cache_file = os.path.join(cache_dir, "tractability_cache.json")
    failures_file = os.path.join(cache_dir, "api_failures.log")
    
    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            cache = json.load(f)
    else:
        cache = {}
        
    results = []
    
    for gene in genes:
        if gene in cache:
            cached_result = dict(cache[gene])
            cached_result['Source'] = 'Cache'
            results.append(cached_result)
            continue
            
        # First resolve symbol to Ensembl ID
        query_search = """
        query search($queryString: String!) {
          search(queryString: $queryString, entityNames: ["target"]) {
            hits {
              id
              object {
                ... on Target {
                  approvedSymbol
                }
              }
            }
          }
        }
        """
        try:
            r = requests.post(OPENTARGETS_API_URL, json={'query': query_search, 'variables': {'queryString': gene}}, headers=headers)
            hits = r.json().get('data', {}).get('search', {}).get('hits', [])
            ensembl_id = None
            for hit in hits:
                if hit.get('object', {}).get('approvedSymbol') == gene:
                    ensembl_id = hit.get('id')
                    break
                    
            if not ensembl_id:
                continue
                
            # Now query tractability
            query_trac = """
            query target($ensemblId: String!){
              target(ensemblId: $ensemblId){
                tractability {
                  id
                  modality
                  value
                }
              }
            }
            """
            r = requests.post(OPENTARGETS_API_URL, json={'query': query_trac, 'variables': {'ensemblId': ensembl_id}}, headers=headers)
            data = r.json()
            tractability = data.get('data', {}).get('target', {}).get('tractability', [])
            
            # We filter for modalities that have 'value' == True (i.e. considered tractable)
            tractable_modalities = set([t['modality'] for t in tractability if t.get('value') is True])
            
            result = {
                'Target_Gene': gene,
                'Small_Molecule_Tractable': 'SM' in tractable_modalities,
                'Antibody_Tractable': 'AB' in tractable_modalities,
                'Other_Modalities_Tractable': 'PR' in tractable_modalities or 'Other' in tractable_modalities,
            }
            cache[gene] = result
            
            # For the current run's dataframe, explicitly mark it as API
            current_run_result = dict(result)
            current_run_result['Source'] = 'API'
            results.append(current_run_result)
            
        except Exception as e:
            msg = f"{datetime.now().isoformat()} - [Tractability] Failed for {gene}: {e}"
            print(f"[Tractability Warning] Failed for {gene}: {e}. Logged to failures.")
            with open(failures_file, "a") as f:
                f.write(msg + "\n")
            continue
            raise
            
    # Save cache
    with open(cache_file, 'w') as f:
        json.dump(cache, f, indent=4)
        
    df = pd.DataFrame(results)
    if not df.empty:
        sources = df['Source'].value_counts().to_dict()
        print(f"[Tractability] Found tractability data for {len(df)} genes.")
        print(f"[Tractability] Data Sources -> {sources}")
    return df
