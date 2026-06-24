import os
import json
import subprocess

BASE_DIR = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
INPUT_DIR = os.path.join(BASE_DIR, 'input')
CONFIG_PATH = os.path.join(INPUT_DIR, 'pipeline.config.json')

def fetch_reactome_pathway(pathway_id, key):
    url = f"https://reactome.org/ContentService/data/participants/{pathway_id}/referenceEntities"
    genes = []
    
    try:
        # Using curl via subprocess to bypass Python HTTP/2 protocol errors in urllib/requests
        result = subprocess.run(["curl", "-s", "-H", "Accept: application/json", url], capture_output=True, text=True, check=True)
        if not result.stdout.strip():
            raise RuntimeError(f"Empty response from {url}")
            
        data = json.loads(result.stdout)
        for item in data:
            if 'geneName' in item and item['geneName']:
                genes.extend(item['geneName'])
            elif 'name' in item and item['name']:
                genes.append(item['name'][0])
        # Filter out Non-HGNC / None values and unique
        genes = list(set([str(g).upper() for g in genes if g]))
        
        filename = f"reactome_{pathway_id}_{key}.json"
        filepath = os.path.join(INPUT_DIR, filename)
        with open(filepath, 'w') as f:
            json.dump(genes, f, indent=4)
        print(f"Successfully saved {len(genes)} genes for {key} ({pathway_id}) to {filepath}")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error fetching {pathway_id} via curl: {e}")
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse JSON for {pathway_id}: {e}\nResponse preview: {result.stdout[:200]}")
        raise RuntimeError(f"Failed to fetch {pathway_id}: {e}")

if __name__ == "__main__":
    if not os.path.exists(CONFIG_PATH):
        raise FileNotFoundError(f"Configuration file not found: {CONFIG_PATH}")
        
    with open(CONFIG_PATH, 'r') as f:
        config = json.load(f)
        
    serotonin_axis = config.get("SEROTONIN_AXIS", {})
    reactome_pathways = serotonin_axis.get("REACTOME_PATHWAYS", {})
    
    for key, info in reactome_pathways.items():
        pathway_id = info["id"]
        print(f"Fetching Reactome Pathway: {key} ({pathway_id})")
        fetch_reactome_pathway(pathway_id, key)
