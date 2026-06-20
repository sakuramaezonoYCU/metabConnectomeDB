import pandas as pd
import numpy as np
import os
import sys

# Ensure the scripts directory is in path to import our modular scripts
if '..' not in sys.path: sys.path.append('..')
from oxygen_tension_config import GLYCOLYSIS_GENES, OXPHOS_GENES, HIF1_GENES, BASE_DIR
from pan_cancer_config import CANCERS_TO_RUN, CANCER_PO2_CSV_MAPPING, normalize_cancer_name

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
    
    # Load the physiological oxygen tension reference CSV
    csv_path = os.path.join(BASE_DIR, "input", "pO2_guide_24588669.csv")
    try:
        po2_df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"[Error] Failed to read {csv_path}: {e}")
        po2_df = pd.DataFrame()

    for raw_cancer in CANCERS_TO_RUN:
        cancer = normalize_cancer_name(raw_cancer)
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
        hif1_df = df[df['names'].isin(HIF1_GENES)]
        
        # Calculate mean Log2 Fold Change (Metastasis vs Primary)
        # Positive LFC means upregulated in Metastasis
        mean_glyco_lfc = glyco_df['logfoldchanges'].mean() if not glyco_df.empty else 0
        mean_oxphos_lfc = oxphos_df['logfoldchanges'].mean() if not oxphos_df.empty else 0
        mean_hif1_lfc = hif1_df['logfoldchanges'].mean() if not hif1_df.empty else 0
        
        # Compute OXPHOS / Glycolysis Enrichment Ratio in Linear Scale
        # LFC difference: LFC_oxphos - LFC_glyco
        # Then linear ratio: 2^(LFC_oxphos - LFC_glyco)
        lfc_diff = mean_oxphos_lfc - mean_glyco_lfc
        oxphos_glyco_ratio = 2 ** lfc_diff
        
        # Fetch oxygen tension data
        csv_name = CANCER_PO2_CSV_MAPPING.get(cancer, "")
        o2_tumour = np.nan
        o2_normal = np.nan
        pmid = ""
        
        if csv_name and not po2_df.empty:
            row = po2_df[po2_df['Tumour type'] == csv_name]
            if not row.empty:
                try:
                    o2_tumour = float(row['Median % oxygen_tumour'].values[0])
                except Exception:
                    pass
                try:
                    o2_normal = float(row['Median % oxygen_normal'].values[0])
                except Exception:
                    pass
                pmid = str(row['Reference'].values[0]).replace('"', '')
                
                # Default empty normal tissue to 6.0 average
                if pd.isna(o2_normal) and pd.notna(o2_tumour):
                    o2_normal = 6.0
        
        results.append({
            "Cancer": cancer.capitalize(),
            "Mean_Glycolysis_LFC": mean_glyco_lfc,
            "Mean_OXPHOS_LFC": mean_oxphos_lfc,
            "Mean_HIF1_LFC": mean_hif1_lfc,
            "OXPHOS_Glycolysis_Ratio": oxphos_glyco_ratio,
            "O2_Tension_Tumour_Pct": o2_tumour,
            "O2_Tension_Normal_Pct": o2_normal,
            "PMID_Reference": pmid,
            "Glycolysis_Genes_Found": len(glyco_df),
            "OXPHOS_Genes_Found": len(oxphos_df),
            "HIF1_Genes_Found": len(hif1_df)
        })
        
    return pd.DataFrame(results)

if __name__ == "__main__":
    df_res = compute_enrichment_ratios()
    print(df_res)
