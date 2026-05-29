import pandas as pd
import numpy as np
import os
import json
from oxygen_tension_config import CANCERS, GLYCOLYSIS_GENES, OXPHOS_GENES, OXYGEN_TENSION_MAP, BASE_DIR

def get_de_file_path(cancer):
    # Depending on how the output folder is structured, try standard locations
    paths_to_try = [
        os.path.join(BASE_DIR, "output", f"{cancer}_results", f"primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv"),
        os.path.join(BASE_DIR, "output", f"primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv")
    ]
    for path in paths_to_try:
        if os.path.exists(path):
            return path
    return None

def compute_enrichment_ratios():
    results = []
    
    for cancer in CANCERS:
        file_path = get_de_file_path(cancer)
        
        if file_path is None:
            print(f"[Warning] Data for {cancer} not found. Skipping...")
            continue
            
        try:
            df = pd.read_csv(file_path)
        except Exception as e:
            print(f"[Error] Failed to read {file_path}: {e}")
            continue
            
        # Ensure gene names are uppercase for string matching
        df['names'] = df['names'].str.upper()
        
        # Filter sets
        glyco_df = df[df['names'].isin(GLYCOLYSIS_GENES)]
        oxphos_df = df[df['names'].isin(OXPHOS_GENES)]
        
        # Calculate mean Log2 Fold Change (Metastasis vs Primary)
        # Positive LFC means upregulated in Metastasis
        mean_glyco_lfc = glyco_df['logfoldchanges'].mean() if not glyco_df.empty else 0
        mean_oxphos_lfc = oxphos_df['logfoldchanges'].mean() if not oxphos_df.empty else 0
        
        # Compute OXPHOS / Glycolysis Enrichment Ratio in Linear Scale
        # LFC difference: LFC_oxphos - LFC_glyco
        # Then linear ratio: 2^(LFC_oxphos - LFC_glyco)
        lfc_diff = mean_oxphos_lfc - mean_glyco_lfc
        oxphos_glyco_ratio = 2 ** lfc_diff
        
        cancer_name = cancer.capitalize()
        o2_tension = OXYGEN_TENSION_MAP.get(cancer_name, np.nan)
        
        results.append({
            "Cancer": cancer_name,
            "Mean_Glycolysis_LFC": mean_glyco_lfc,
            "Mean_OXPHOS_LFC": mean_oxphos_lfc,
            "OXPHOS_Glycolysis_Ratio": oxphos_glyco_ratio,
            "O2_Tension_Pct": o2_tension,
            "Glycolysis_Genes_Found": len(glyco_df),
            "OXPHOS_Genes_Found": len(oxphos_df)
        })
        
    return pd.DataFrame(results)

if __name__ == "__main__":
    df_res = compute_enrichment_ratios()
    print(df_res)
