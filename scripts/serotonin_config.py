"""
serotonin_config.py
===================
Centralized configurations and gene lists for the Serotonin-HTR7-TAM-IP4-HR repair axis.
All lists are dynamically loaded from curated JSON files or fetched Reactome APIs.
"""
import os
import json

BASE_DIR = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
INPUT_DIR = os.path.join(BASE_DIR, 'input')
CONFIG_PATH = os.path.join(INPUT_DIR, 'pipeline.config.json')

def load_json(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Missing required configuration file: {filepath}")
    with open(filepath, 'r') as f:
        return json.load(f)

# 1. Load the main pipeline config
config = load_json(CONFIG_PATH)
serotonin_axis = config.get("SEROTONIN_AXIS", {})

# 2. Load the HTR7 TAM Signature
htr7_file = os.path.join(BASE_DIR, serotonin_axis.get("HTR7_TAM_SIGNATURE_FILE", "input/htr7_tam_signature.json"))
HTR7_TAM_SIGNATURE = load_json(htr7_file)

# 3. Load Curated Serotonin Targets
curated_file = os.path.join(BASE_DIR, serotonin_axis.get("CURATED_SEROTONIN_TARGETS_FILE", "input/serotonin_curated_targets.json"))
curated_targets = load_json(curated_file)

SEROTONIN_SYNTHESIS = curated_targets.get("SEROTONIN_SYNTHESIS", [])
SEROTONIN_TRANSPORT = curated_targets.get("SEROTONIN_TRANSPORT", [])
SEROTONIN_DEGRADATION = curated_targets.get("SEROTONIN_DEGRADATION", [])
SEROTONIN_RECEPTORS = curated_targets.get("SEROTONIN_RECEPTORS", [])
PARACRINE_PAIRS = curated_targets.get("PARACRINE_PAIRS", [])
EV_CARGO_SORTING = curated_targets.get("EV_CARGO_SORTING", [])
SEROTONIN_DEP_SECRETION = curated_targets.get("SEROTONIN_DEP_SECRETION", [])

# 4. Load Reactome Pathways
def load_reactome_pathway(key):
    reactome_pathways = serotonin_axis.get("REACTOME_PATHWAYS", {})
    if key not in reactome_pathways:
        return []
    pathway_id = reactome_pathways[key]["id"]
    filepath = os.path.join(INPUT_DIR, f"reactome_{pathway_id}_{key}.json")
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"Missing {filepath}. Please run 'python scripts/fetch_reactome_pathways.py' "
            f"to reproducibly fetch the {key} gene set from the Reactome REST API."
        )
    return load_json(filepath)

INOSITOL_PATHWAY = load_reactome_pathway("INOSITOL_PATHWAY")
HR_REPAIR_GENES = load_reactome_pathway("HR_REPAIR")
EV_BIOGENESIS = load_reactome_pathway("EV_BIOGENESIS")
EXHAUSTION_TARGETS = load_reactome_pathway("EXHAUSTION_TARGETS")
SUPPRESSIVE_LIGAND_TARGETS = load_reactome_pathway("SUPPRESSIVE_LIGAND_TARGETS")
TREG_TARGETS = load_reactome_pathway("TREG_TARGETS")

def get_all_genes():
    """Returns a flat list of all genes tracked in this module."""
    all_genes = (
        HTR7_TAM_SIGNATURE + 
        SEROTONIN_SYNTHESIS + 
        SEROTONIN_TRANSPORT + 
        SEROTONIN_DEGRADATION + 
        SEROTONIN_RECEPTORS + 
        INOSITOL_PATHWAY + 
        HR_REPAIR_GENES + 
        EXHAUSTION_TARGETS +
        SUPPRESSIVE_LIGAND_TARGETS +
        TREG_TARGETS +
        EV_BIOGENESIS + 
        EV_CARGO_SORTING + 
        SEROTONIN_DEP_SECRETION
    )
    return list(set(all_genes))
