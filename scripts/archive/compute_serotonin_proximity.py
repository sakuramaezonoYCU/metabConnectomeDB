"""
compute_serotonin_proximity.py
==============================
Wrapper to run the Serotonin spatial proximity analysis 
(from compute_serotonin_spatial.py) across all pan-cancer datasets.
"""
import os
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from compute_serotonin_spatial import run_spatial_analysis
from pan_cancer_config import CANCER_CAP

def compute_all_proximity():
    """Run spatial proximity for all cancers."""
    cancers = list(CANCER_CAP.keys())
    
    print(f"Starting Pan-Cancer Serotonin Proximity Analysis for {len(cancers)} cancers...")
    for cancer in cancers:
        print(f"\n{'='*50}")
        try:
            run_spatial_analysis(cancer)
        except Exception as e:
            print(f"❌ Error processing {cancer}: {e}")
            
    print(f"\n{'='*50}\nPan-Cancer Analysis Complete.")

if __name__ == "__main__":
    compute_all_proximity()
