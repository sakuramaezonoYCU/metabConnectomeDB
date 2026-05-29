import requests
import pandas as pd
from IPython.display import display
from druggability_config import OPENTARGETS_API_URL

def query_string_ppi(genes):
    """
    Query STRING Database for Protein-Protein Interactions among the target genes.
    This helps establish if the 'axis' has physical or functional evidence of interaction.
    """
    print(f"[STRING PPI] Querying network for: {genes}")
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
            print("[STRING PPI] No interactions found among these genes in STRING.")
            return pd.DataFrame()
            
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
        print(f"[STRING PPI] Found {len(df)} interaction(s).")
        return df
    except Exception as e:
        print(f"[STRING PPI Error] Failed to query STRING: {e}")
        return pd.DataFrame()

def query_tractability(genes):
    """
    Query Open Targets for Target Tractability.
    Even if no drugs currently exist, tractability predicts if the target 
    is chemically druggable (small molecule) or targetable by antibodies.
    """
    print(f"[Tractability] Querying Open Targets for tractability of {genes}...")
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) metabConnectomeDB'}
    
    results = []
    
    for gene in genes:
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
            
            results.append({
                'Target_Gene': gene,
                'Small_Molecule_Tractable': 'SM' in tractable_modalities,
                'Antibody_Tractable': 'AB' in tractable_modalities,
                'Other_Modalities_Tractable': 'PR' in tractable_modalities or 'Other' in tractable_modalities
            })
            
        except Exception as e:
            print(f"[Tractability Error] Failed for {gene}: {e}")
            
    df = pd.DataFrame(results)
    if not df.empty:
        print(f"[Tractability] Found tractability data for {len(df)} genes.")
    return df
