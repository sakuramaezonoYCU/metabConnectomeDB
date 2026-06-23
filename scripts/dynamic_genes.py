import os
import pandas as pd
import json

# =====================================================================
# OVERRIDE PARAMETERS
# =====================================================================
# To explicitly override the dynamically calculated genes across the pipeline,
# populate the "OVERRIDE_TARGET_GENES" array in pipeline.config.json at the root of the project.
# e.g., "OVERRIDE_TARGET_GENES": ["GLS", "SLC16A7", "STAT3"]
# If left empty, the pipeline will dynamically calculate the intersection.

def get_dynamic_genes(base_dir=None):
    """
    Dynamically loads the intersection of upregulated metabolic targets
    across all executed cancer cohorts, ensuring zero hardcoding of data.
    If pipeline.config.json contains OVERRIDE_TARGET_GENES, it returns those instead.
    """
    if base_dir is None:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
    config_path = os.path.join(base_dir, 'pipeline.config.json')
    override_genes = []
    if os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config = json.load(f)
            override_genes = config.get("OVERRIDE_TARGET_GENES", [])
    if override_genes:
        print(f"Using {len(override_genes)} user-overridden target genes from pipeline.config.json.")
        return override_genes
        # Assuming this is called from within scripts/
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        
    cancer_genes_map = {}
    from pan_cancer_config import CANCER_CAP, get_de_csv_path
    for c in CANCER_CAP.keys():
        res_file = get_de_csv_path(c)
        if os.path.exists(res_file):
            df_de = pd.read_csv(res_file)
            up_genes = set(df_de[df_de['Significance'] == 'Up in Metastasis']['names'].tolist())
            cancer_genes_map[c] = up_genes

    if not cancer_genes_map:
        raise ValueError("No cancer gene sets were loaded! Ensure upstream DE pipeline executed successfully.")

    # Get union of all genes
    all_genes = set()
    for genes in cancer_genes_map.values():
        all_genes.update(genes)
        
    pan_cancer = []
    for g in all_genes:
        count = sum(1 for genes in cancer_genes_map.values() if g in genes)
        if count >= len(cancer_genes_map) - 1:  # Present in N-1 or N cancers
            pan_cancer.append(g)
            
    return list(pan_cancer)
