import sys
import os
import builtins
import pandas as pd
import numpy as np

# Override globals
builtins.OUTPUT_DIR = os.path.abspath('output/breast_results')
builtins.h5ad_path = os.path.join(builtins.OUTPUT_DIR, 'breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad')
builtins.PRIMARY_TISSUES = ['breast']
builtins.CANCER_TYPE_NAME = 'breast_cancer'
builtins.GLOBAL_OUTPUT_DIR = os.path.abspath('output')

print(f"Testing with: {builtins.h5ad_path}")
print(f"File exists: {os.path.exists(builtins.h5ad_path)}")

# Run the script
with open('scripts/primary_vs_metastasis_comparison.py', 'r') as f:
    script_content = f.read()
    
exec(script_content, globals())
