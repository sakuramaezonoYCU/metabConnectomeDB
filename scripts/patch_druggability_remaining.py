import nbformat as nbf
import os

nb_path = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/druggability_axis_analysis.ipynb"

try:
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    modified = False
    for cell in nb.cells:
        old_source = cell.source
        
        # Replace 23s and 181s
        cell.source = cell.source.replace("the 181 and 23 core pan-cancer", "the broadly and strictly conserved pan-cancer")
        cell.source = cell.source.replace("**23 pan-cancer genes**", "**strictly conserved pan-cancer genes**")
        cell.source = cell.source.replace("**181 genes**", "**broadly conserved genes**")
        cell.source = cell.source.replace("for 23 conserved genes", "for strictly conserved genes")
        cell.source = cell.source.replace("for 181 broadly conserved genes", "for broadly conserved genes")
        
        if cell.source != old_source:
            modified = True
            
        # Clear outputs
        if cell.cell_type == 'code':
            cell.outputs = []
            
    if modified:
        with open(nb_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched druggability_axis_analysis.ipynb")

except Exception as e:
    print(e)
