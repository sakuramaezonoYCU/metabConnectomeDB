import os
import requests
import json
import pandas as pd

# Define paths
META_RESULTS_DIR = 'output/pan_cancer_meta_results'
OUTPUT_DIR = 'output'
ANALYSIS_SUFFIX = '_Br500k_Co100k_Lu500k_Me100k_Ov100k'

def fetch_opentargets_diseases():
    conserved_csv = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
    if not os.path.exists(conserved_csv):
        print(f"File not found: {conserved_csv}")
        return
        
    sig_df = pd.read_csv(conserved_csv)
    genes = sig_df['Strictly_Conserved_Gene'].tolist()
    
    results = []
    
    print("Fetching disease associations from Open Targets Platform GraphQL API...")
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    
    # GraphQL Query to search for gene symbol, get its Ensembl ID, and fetch top 3 associated diseases
    query = """
    query targetSearch($queryString: String!) {
      search(queryString: $queryString, entityNames: ["target"]) {
        hits {
          id
          name
          object {
            ... on Target {
              associatedDiseases(page: {index: 0, size: 3}) {
                count
                rows {
                  disease {
                    name
                  }
                  score
                }
              }
            }
          }
        }
      }
    }
    """
    
    for gene in genes:
        try:
            response = requests.post(url, json={"query": query, "variables": {"queryString": gene}}, headers={'User-Agent': 'python-requests/2.25.1'})
            response.raise_for_status()
            data = response.json()
            
            diseases_text = "No diseases found"
            if 'data' in data and data['data']['search']['hits']:
                # Get the first hit (best match for gene symbol)
                hit = data['data']['search']['hits'][0]
                target_obj = hit.get('object', {})
                if target_obj and 'associatedDiseases' in target_obj:
                    disease_rows = target_obj['associatedDiseases']['rows']
                    if disease_rows:
                        disease_names = [r['disease']['name'] for r in disease_rows if 'disease' in r and 'name' in r['disease']]
                        if disease_names:
                            diseases_text = ", ".join(disease_names)
            
            results.append({
                "Gene": gene,
                "OpenTargets_Diseases": diseases_text
            })
            print(f"Fetched {gene}: {diseases_text}")
            
        except Exception as e:
            print(f"Error fetching {gene}: {e}")
            results.append({
                "Gene": gene,
                "OpenTargets_Diseases": "Error fetching from API"
            })
            
    # Save to CSV
    df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, 'opentargets_diseases.csv')
    df.to_csv(out_path, index=False)
    print(f"\nSaved OpenTargets disease associations to {out_path}")

if __name__ == "__main__":
    fetch_opentargets_diseases()
