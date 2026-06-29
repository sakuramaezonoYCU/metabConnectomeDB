import os
import sys
import json
import requests

def fetch_kegg_genes(pathway_id):
    """
    Fetches the list of genes associated with a specific KEGG pathway using the KEGG REST API.
    """
    url = f"http://rest.kegg.jp/get/{pathway_id}"
    print(f"Fetching {pathway_id} from {url}...")
    
    response = requests.get(url)
    response.raise_for_status()
    
    genes = []
    in_gene_section = False
    
    for line in response.text.split('\n'):
        if line.startswith('GENE'):
            in_gene_section = True
            line = line[12:] # Strip "GENE        "
        elif in_gene_section:
            if not line.startswith(' '):
                break # End of gene section
            line = line.strip()
        else:
            continue
            
        # Parse the HGNC symbol from the KEGG GENE format
        # Format example: "2821  GPI; glucose-6-phosphate isomerase [KO:K01810] [EC:5.3.1.9]"
        parts = line.split(';')
        if parts:
            first_part = parts[0].strip()
            tokens = first_part.split()
            if len(tokens) > 1:
                # The primary symbol is the second token (first is numeric ID)
                symbol = tokens[1].strip(',')
                if not symbol.startswith('[KO'):
                    genes.append(symbol)

    unique_genes = sorted(list(set(genes)))
    print(f"  -> Found {len(unique_genes)} unique genes for {pathway_id}")
    return unique_genes

if __name__ == "__main__":
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    INPUT_DIR = os.path.join(BASE_DIR, 'input')
    os.makedirs(INPUT_DIR, exist_ok=True)
    
    # Import KEGG_PATHWAYS from our centralized config
    from pan_cancer_config import KEGG_PATHWAYS
    
    for pathway_key, pathway_info in KEGG_PATHWAYS.items():
        kegg_id = pathway_info["id"]
        # E.g. hsa00010 -> kegg_hsa00010_GLYCOLYSIS.json
        filename = f"kegg_{kegg_id}_{pathway_key}.json"
        out_path = os.path.join(INPUT_DIR, filename)
        
        if os.path.exists(out_path):
            print(f"Skipping {pathway_key} - {out_path} already exists.")
            continue
            
        genes = fetch_kegg_genes(kegg_id)
        with open(out_path, 'w') as f:
            json.dump(genes, f, indent=4)
        print(f"Saved {pathway_key} to {out_path}")

    print("\n✅ Successfully fetched and updated KEGG pathway gene sets.")
