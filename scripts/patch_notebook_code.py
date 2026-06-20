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
            
            if "upset_plot.png" in cell.source:
                if "ANALYSIS_SUFFIX" not in cell.source:
                    cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                cell.source = cell.source.replace("'output/pan_cancer_meta_results/upset_plot.png'", "f'output/pan_cancer_meta_results/upset_plot{ANALYSIS_SUFFIX}.png'")
                
            if "tissue_specific_correlation.png" in cell.source:
                if "ANALYSIS_SUFFIX" not in cell.source:
                    cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                cell.source = cell.source.replace("'output/pan_cancer_meta_results/tissue_specific_correlation.png'", "f'output/pan_cancer_meta_results/tissue_specific_correlation{ANALYSIS_SUFFIX}.png'")

            if "pan_cancer_23_genes.csv" in cell.source:
                if "ANALYSIS_SUFFIX" not in cell.source:
                    cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                cell.source = cell.source.replace("'output/pan_cancer_meta_results/pan_cancer_23_genes.csv'", "f'output/pan_cancer_meta_results/pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv'")

            if cell.source != old_source:
                modified = True
                
    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched pan_cancer_meta_analysis.ipynb code cells")

except Exception as e:
    print(f"Error: {e}")
