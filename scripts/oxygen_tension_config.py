import os
import sys

# Ensure the scripts directory is in path to import our modular scripts
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import CANCERS_TO_RUN, CANCER_PO2_CSV_MAPPING, normalize_cancer_name, KEGG_PATHWAYS

import json

# Paths
BASE_DIR = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
INPUT_DIR = os.path.join(BASE_DIR, 'input')

def load_kegg_genes(pathway_key):
    pathway_info = KEGG_PATHWAYS[pathway_key]
    kegg_id = pathway_info["id"]
    filename = f"kegg_{kegg_id}_{pathway_key}.json"
    filepath = os.path.join(INPUT_DIR, filename)
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"Missing {filename}. Please run 'python scripts/fetch_kegg_pathways.py' "
            f"to reproducibly fetch the {pathway_key} gene set from the KEGG REST API."
        )
    with open(filepath, 'r') as f:
        return json.load(f)

# Core OXPHOS and Glycolysis metabolic gene sets (Reproducibly derived from KEGG REST API)
GLYCOLYSIS_GENES = load_kegg_genes("GLYCOLYSIS")
OXPHOS_GENES = load_kegg_genes("OXPHOS")
HIF1_GENES = load_kegg_genes("HIF1")

OUTPUT_DIR = os.path.join(BASE_DIR, "output", "oxygen_tension")
os.makedirs(OUTPUT_DIR, exist_ok=True)
