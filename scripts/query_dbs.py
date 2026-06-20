import requests
import pandas as pd
import time
from pan_cancer_config import DGIDB_API_URL, OPENTARGETS_API_URL

headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) metabConnectomeDB'}

def query_dgidb(genes):
    """
    Query the Drug Gene Interaction Database (DGIdb) for a list of genes using GraphQL.
    """
    print(f"[DGIdb] Querying interactions for {len(genes)} genes...")
    results = []
    
    DGIDB_GQL = "https://dgidb.org/api/graphql"
    genes_str = '", "'.join(genes)
    q = f"""
    query {{
      genes(names: ["{genes_str}"]) {{
        nodes {{
          name
          interactions {{
            drug {{
              name
            }}
            interactionTypes {{
              type
            }}
            interactionAttributes {{
              name
              value
            }}
          }}
        }}
      }}
    }}
    """
    
    try:
        r = requests.post(DGIDB_GQL, json={"query": q}, headers=headers)
        r.raise_for_status()
        data = r.json()
        
        nodes = data.get('data', {}).get('genes', {}).get('nodes', [])
        for node in nodes:
            gene_name = node.get('name')
            for interaction in node.get('interactions', []):
                drug_name = interaction.get('drug', {}).get('name')
                if drug_name:
                    types = [t.get('type') for t in interaction.get('interactionTypes', []) if t.get('type')]
                    results.append({
                        'Database': 'DGIdb',
                        'Target_Gene': gene_name,
                        'Drug_Name': drug_name,
                        'Interaction_Type': ', '.join(types) if types else 'Targeted',
                        'Sources': 'DGIdb GraphQL',
                        'Approval_Status': 'N/A'
                    })
        print(f"[DGIdb] Found {len(results)} interactions.")
        return pd.DataFrame(results)
    except Exception as e:
        import warnings
        warnings.warn(f"[DGIdb Error] Failed to query DGIdb: {e}")
        return pd.DataFrame()

def query_open_targets(genes):
    """
    Query Open Targets Platform GraphQL API for known drugs for specific genes.
    Open Targets aggregates ChEMBL, clinicaltrials.gov, and OMIM/MalaCards diseases.
    """
    print(f"[OpenTargets] Resolving gene symbols to Ensembl IDs for {len(genes)} genes...")
    
    # Step 1: Resolve symbols to Ensembl IDs
    symbol_to_id = {}
    for gene in genes:
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
            data = r.json()
            hits = data.get('data', {}).get('search', {}).get('hits', [])
            for hit in hits:
                if hit.get('object', {}).get('approvedSymbol') == gene:
                    symbol_to_id[gene] = hit.get('id')
                    break
        except Exception as e:
            import warnings
            warnings.warn(f"Error resolving {gene}: {e}")
            continue
            
    print(f"[OpenTargets] Resolved IDs: {symbol_to_id}")
    
    results = []
    # Step 2: Query Known Drugs for each Ensembl ID
    for gene, ensembl_id in symbol_to_id.items():
        print(f"[OpenTargets] Fetching known drugs for {gene} ({ensembl_id})...")
        query_drugs = """
        query target($ensemblId: String!){
          target(ensemblId: $ensemblId){
            drugAndClinicalCandidates {
              rows {
                drug {
                  name
                }
                maxClinicalStage
                diseases {
                  disease {
                    name
                  }
                }
              }
            }
          }
        }
        """
        try:
            r = requests.post(OPENTARGETS_API_URL, json={'query': query_drugs, 'variables': {'ensemblId': ensembl_id}}, headers=headers)
            data = r.json()
            
            if 'errors' in data:
                import warnings
                warnings.warn(f"[OpenTargets GraphQL Error for {gene}]: {data['errors'][0]['message']}")
                continue
                
            drugs = data.get('data', {}).get('target', {})
            if drugs:
                drugs = drugs.get('drugAndClinicalCandidates', {})
                if drugs:
                    drugs = drugs.get('rows', [])
                else:
                    drugs = []
            else:
                drugs = []
            
            for row in drugs:
                drug_info = row.get('drug', {})
                diseases = row.get('diseases', [])
                disease_names = [d.get('disease', {}).get('name') for d in diseases if d.get('disease')]
                
                results.append({
                    'Database': 'OpenTargets',
                    'Target_Gene': gene,
                    'Drug_Name': drug_info.get('name'),
                    'Interaction_Type': 'Targeted (ChEMBL/CT.gov)',
                    'Approval_Status': f"Max Phase {row.get('maxClinicalStage', 'Unknown')}",
                    'Indication_Disease': ', '.join(disease_names)
                })
        except Exception as e:
            import warnings
            warnings.warn(f"[OpenTargets Error] Failed fetching drugs for {gene}: {e}")
            continue
            
            
    print(f"[OpenTargets] Found {len(results)} known drug indications.")
    return pd.DataFrame(results)

def query_diseases(genes):
    """
    Query disease associations for the target genes.
    We use the OpenTargets GraphQL API, which natively integrates JensenLab DISEASES, 
    Orphanet, and other disease-gene association databases.
    """
    print(f"[DISEASES] Querying disease associations for {len(genes)} genes via OpenTargets...")
    results = []
    headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) metabConnectomeDB'}
    from pan_cancer_config import OPENTARGETS_API_URL
    import requests
    
    symbol_to_id = {}
    for gene in genes:
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
            for hit in hits:
                if hit.get('object', {}).get('approvedSymbol') == gene:
                    symbol_to_id[gene] = hit.get('id')
                    break
        except Exception as e:
            import warnings
            warnings.warn(f"[DISEASES Error] Error resolving {gene}: {e}")
            continue
            
    for gene, ensembl_id in symbol_to_id.items():
        query_assoc = """
        query targetDiseaseAssociations($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            associatedDiseases(page: { size: 10, index: 0 }) {
              rows {
                score
                disease {
                  name
                }
              }
            }
          }
        }
        """
        try:
            r = requests.post(OPENTARGETS_API_URL, json={'query': query_assoc, 'variables': {'ensemblId': ensembl_id}}, headers=headers)
            data = r.json()
            if 'errors' in data:
                import warnings
                warnings.warn(f"[DISEASES GraphQL Error for {gene}]: {data['errors'][0]['message']}")
                continue
                
            diseases = data.get('data', {}).get('target', {}).get('associatedDiseases', {}).get('rows', [])
            disease_names = [d.get('disease', {}).get('name') for d in diseases if d.get('disease')]
            
            if disease_names:
                results.append({
                    'Database': 'OpenTargets_DISEASES',
                    'Target_Gene': gene,
                    'Drug_Name': 'N/A',
                    'Interaction_Type': 'Disease Association',
                    'Approval_Status': 'N/A',
                    'Indication_Disease': ', '.join(disease_names)
                })
        except Exception as e:
            import warnings
            warnings.warn(f"[DISEASES Error] Failed fetching diseases for {gene}: {e}")
            continue

    df = pd.DataFrame(results)
    print(f"[DISEASES] Found disease associations for {len(df)} genes.")
    return df

def compile_drug_databases(genes):
    """
    Run all database queries and merge results.
    Note: MalaCards and OMIM are effectively queried through OpenTargets which 
    integrates orphan disease and OMIM data under the hood. DISEASES is also 
    integrated in OpenTargets.
    """
    df_dgi = query_dgidb(genes)
    df_ot = query_open_targets(genes)
    df_dis = query_diseases(genes)
    
    combined = pd.concat([df_dgi, df_ot, df_dis], ignore_index=True)
    if combined.empty:
        print("No interactions found across queried databases.")
        return combined
        
    combined.drop_duplicates(subset=['Target_Gene', 'Drug_Name', 'Database'], inplace=True)
    return combined

if __name__ == "__main__":
    from dynamic_genes import get_dynamic_genes
    TARGET_GENES = get_dynamic_genes('.')
    df = compile_drug_databases(TARGET_GENES)
    print(df.head())
