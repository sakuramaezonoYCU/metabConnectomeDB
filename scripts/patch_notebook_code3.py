import nbformat as nbf
import os
import re

filepath = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/pan_cancer_meta_analysis.ipynb"

try:
    with open(filepath, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    modified = False
    for cell in nb.cells:
        if cell.cell_type == 'code':
            old_source = cell.source
            
            # Fix signature score pngs
            import sys
            if '..' not in sys.path: sys.path.append('..')
            try:
                from pan_cancer_config import CANCERS_TO_RUN
                cancers = CANCERS_TO_RUN
            except ImportError:
                cancers = ['breast', 'lung', 'colorectal', 'melanoma', 'ovarian']
                
            for cancer in cancers:
                target = f"'{cancer}_primary_signature_score.png'"
                if target in cell.source:
                    if "ANALYSIS_SUFFIX" not in cell.source:
                        cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                    cell.source = cell.source.replace(target, f"f'{cancer}_primary_signature_score{{ANALYSIS_SUFFIX}}.png'")
                    
            if cell.source != old_source:
                modified = True
                
    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched signature score pngs in pan_cancer_meta_analysis.ipynb")
    else:
        print("No changes made.")

except Exception as e:
    print(f"Error: {e}")
