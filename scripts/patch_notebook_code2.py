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
            
            # Already patched upset_plot and tissue_specific_correlation
            if "metabolite_target_network.png" in cell.source:
                if "ANALYSIS_SUFFIX" not in cell.source:
                    cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                cell.source = cell.source.replace("'output/pan_cancer_meta_results/metabolite_target_network.png'", "f'output/pan_cancer_meta_results/metabolite_target_network{ANALYSIS_SUFFIX}.png'")
                
            if cell.source != old_source:
                modified = True
                
        if cell.cell_type == 'markdown':
            old_source = cell.source
            
            if "metabolite_target_network_edges.csv" in cell.source:
                cell.source = cell.source.replace("metabolite_target_network_edges.csv", "metabolite_target_network_edges{ANALYSIS_SUFFIX}.csv")
                
            cell.source = re.sub(r"([a-z]+)_primary_signature_scores\.csv", r"\1_primary_signature_scores{ANALYSIS_SUFFIX}.csv", cell.source)
            
            if cell.source != old_source:
                modified = True
                
    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched remaining paths in pan_cancer_meta_analysis.ipynb")

except Exception as e:
    print(f"Error: {e}")

# also fix generate_combined_pan_cancer_notebook.py
filepath_py = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/generate_combined_pan_cancer_notebook.py"
try:
    with open(filepath_py, 'r', encoding='utf-8') as f:
        content = f.read()
    
    old_content = content
    content = content.replace("metabolite_target_network_edges.csv", "metabolite_target_network_edges{ANALYSIS_SUFFIX}.csv")
    content = content.replace("metabolite_target_network.png", "metabolite_target_network{ANALYSIS_SUFFIX}.png")
    
    if content != old_content:
        with open(filepath_py, 'w', encoding='utf-8') as f:
            f.write(content)
        print("Patched generate_combined_pan_cancer_notebook.py")
except Exception as e:
    print(f"Error: {e}")
