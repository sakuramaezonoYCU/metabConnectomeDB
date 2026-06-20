import os
import json
import urllib.request

BASE_DIR = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
INPUT_DIR = os.path.join(BASE_DIR, 'input')
CONFIG_PATH = os.path.join(INPUT_DIR, 'pipeline.config.json')

def fetch_reactome_pathway(pathway_id, key):
    url = f"https://reactome.org/ContentService/data/participants/{pathway_id}/referenceEntities"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36', 'Accept': 'application/json'})
    genes = []
    
    # Bypass environment proxies that might be causing the 403 Tunnel Forbidden
    proxy_handler = urllib.request.ProxyHandler({})
    opener = urllib.request.build_opener(proxy_handler)
    
    try:
        with opener.open(req) as response:
            if response.status == 200:
                data = json.loads(response.read().decode('utf-8'))
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
            else:
                print(f"Error fetching {pathway_id}: HTTP {response.status}")
    except Exception as e:
        print(f"Failed to fetch {pathway_id}: {e}")

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
