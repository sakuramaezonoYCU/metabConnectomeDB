import os
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import json
from pan_cancer_config import DEPMAP_DATA_PATH, ANALYSIS_SUFFIX

# ---------------------------------------------------------
# PATH CONFIGURATION
# ---------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "druggability")
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
OUTPUT_BASENAME = f"druggability_axis{ANALYSIS_SUFFIX}"

def analyze_depmap_synergy(genes):
    """
    Analyzes DepMap CRISPR effect data for a set of genes to identify co-dependencies.
    High positive correlation indicates they operate in the same pathway.
    Negative correlation may indicate synthetic lethality.
    """
    if not os.path.exists(DEPMAP_DATA_PATH):
        print(f"[DepMap] Data not found at {DEPMAP_DATA_PATH}.")
        print("[DepMap] To enable synergy analysis, download 'CRISPRGeneEffect.csv' from the DepMap portal (https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2026Q1&filename=CRISPRGeneEffect.csv) and place it in the input/ folder.")
        return None
        
    print(f"[DepMap] Loading DepMap CRISPR data from {DEPMAP_DATA_PATH}...")
    try:
        # We only need the columns for our specific genes.
        # Let's read the first few rows to find the exact column names.
        df_head = pd.read_csv(DEPMAP_DATA_PATH, nrows=5)
        
        # The first column is the index (ModelID, DepMap_ID, or Unnamed: 0)
        index_col_name = df_head.columns[0]
        
        # Map our gene symbols to the actual column names
        col_mapping = {}
        for col in df_head.columns:
            symbol = col.split(' ')[0]
            if symbol in genes:
                col_mapping[symbol] = col
                
        if not col_mapping:
            print("[DepMap] None of the specified genes were found in the DepMap dataset.")
            return None
            
        print(f"[DepMap] Found columns for: {list(col_mapping.keys())}")
        
        # Read only the necessary columns
        cols_to_use = [index_col_name] + list(col_mapping.values())
        df = pd.read_csv(DEPMAP_DATA_PATH, usecols=cols_to_use)
        df.set_index(index_col_name, inplace=True)
        
        # Rename columns to just the gene symbol for easier plotting
        df.rename(columns={v: k for k, v in col_mapping.items()}, inplace=True)
        
        # Calculate correlation matrix
        corr_matrix = df.corr()
        
        # Plot correlation heatmap
        plt.figure(figsize=(8, 6))
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
        plt.title('DepMap CRISPR Essentiality Correlation (Co-dependency)')
        plt.tight_layout()
        
        plot_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_depmap_correlation.png")
        plt.savefig(plot_path, dpi=300)
        print(f"[DepMap] Saved correlation plot to {plot_path}")
        plt.show()  # Display inline for the notebook/HTML
        plt.close()
        
        return corr_matrix
        
    except Exception as e:
        print(f"[DepMap Error] Failed to analyze DepMap data: {e}")
        return None

if __name__ == "__main__":
    from dynamic_genes import get_dynamic_genes
    TARGET_GENES = get_dynamic_genes('.')
    print(f"Target genes loaded: {TARGET_GENES}")
    analyze_depmap_synergy(TARGET_GENES)
