import nbformat as nbf
import os
import re

filepath = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/pan_cancer_meta_analysis.ipynb"

try:
    with open(filepath, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    modified = False
    for cell in nb.cells:
        if cell.cell_type == 'markdown':
            # Fix hardcoded cancer CSV paths
            old_source = cell.source
            
            # Replace lines like: - **Breast:** `output/breast_results/primary_vs_metastasis_breast_results_DE_metabolic_targets.csv`
            # with: - **Breast:** `output/breast_results/primary_vs_metastasis_breast_results_DE_metabolic_targets_{cancer_suffix}.csv`
            
            cell.source = re.sub(r"(output/[a-z]+_results/primary_vs_metastasis_[a-z]+_results_DE_metabolic_targets)\.csv", r"\1_{cancer_suffix}.csv", cell.source)
            
            # Fix "23 strictly conserved" text
            cell.source = cell.source.replace("23 strictly conserved pan-cancer metabolic targets", "strictly conserved pan-cancer metabolic targets")
            
            # Fix upset_plot_data.csv
            cell.source = cell.source.replace("upset_plot_data.csv", "upset_plot_data{ANALYSIS_SUFFIX}.csv")
            
            if cell.source != old_source:
                modified = True
                
    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched pan_cancer_meta_analysis.ipynb markdown cells")

except Exception as e:
    print(f"Error: {e}")

# Now fix generate_combined_pan_cancer_notebook.py
filepath_py = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/generate_combined_pan_cancer_notebook.py"
try:
    with open(filepath_py, 'r', encoding='utf-8') as f:
        content = f.read()
        
    old_content = content
    content = re.sub(r"(output/[a-z]+_results/primary_vs_metastasis_[a-z]+_results_DE_metabolic_targets)\.csv", r"\1_{cancer_suffix}.csv", content)
    content = content.replace("upset_plot_data.csv", "upset_plot_data{ANALYSIS_SUFFIX}.csv")
    content = content.replace("23 strictly conserved pan-cancer", "strictly conserved pan-cancer")
    
    if content != old_content:
        with open(filepath_py, 'w', encoding='utf-8') as f:
            f.write(content)
        print("Patched generate_combined_pan_cancer_notebook.py markdown generation strings")
except Exception as e:
    print(f"Error: {e}")
