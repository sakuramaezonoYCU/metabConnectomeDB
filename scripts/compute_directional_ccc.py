import sys
import os
import pandas as pd
import numpy as np

# Compute BASE_DIR first
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path: 
    sys.path.append(BASE_DIR)

from pan_cancer_config import ANALYSIS_SUFFIX

def compute_directional_ccc():
    # 1. Load strictly conserved pan-cancer genes
    meta_path = os.path.join(BASE_DIR, "output", "pan_cancer_meta_results", f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv")
    if not os.path.exists(meta_path):
        print("Error: strictly conserved genes list not found.")
        return
    pan_cancer = pd.read_csv(meta_path)
    
    # 2. Assign classes
    def assign_class(gene):
        # Known producers/enzymes from the strictly-conserved list
        producers = ['GLS', 'SGMS1', 'SPTLC1']
        if gene in producers:
            return 'producer'
        else:
            return 'consumer'
            
    pan_cancer['Direction_Class'] = pan_cancer['Strictly_Conserved_Gene'].apply(assign_class)
    
    out_dir = os.path.join(BASE_DIR, "output", f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}", "directional_ccc")
    os.makedirs(out_dir, exist_ok=True)
    
    # Save the classification
    classification_path = os.path.join(out_dir, "metalinks_direction_classes.csv")
    pan_cancer.to_csv(classification_path, index=False)
    print(f"Saved direction classification to {classification_path}")
    
if __name__ == "__main__":
    compute_directional_ccc()
