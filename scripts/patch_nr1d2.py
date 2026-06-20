import nbformat as nbf
import os

nb_path = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/nr1d2_master_regulator_analysis.ipynb"

try:
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    modified = False
    for cell in nb.cells:
        old_source = cell.source
        
        # Replace 23s
        cell.source = cell.source.replace("the 23 pan-cancer conserved", "the strictly conserved pan-cancer")
        cell.source = cell.source.replace("these 23 genes", "these strictly conserved genes")
        
        if cell.source != old_source:
            modified = True
            
        # Clear outputs
        if cell.cell_type == 'code':
            cell.outputs = []
            
    if modified:
        with open(nb_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched nr1d2_master_regulator_analysis.ipynb")

except Exception as e:
    print(e)
