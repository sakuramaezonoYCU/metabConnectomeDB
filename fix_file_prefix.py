import nbformat as nbf
import os
import re

# 1. Update cancer_cellxgene_integration.ipynb to remove 'cancer_disease-name_' prefix from outputs
notebook_path = 'scripts/cancer_cellxgene_integration.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

modified = False
for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'def get_cache_filename' in cell.source:
            # We want to remove the disease_slug prefix from the filename returned
            new_source = re.sub(
                r'return f"\{disease_slug\}_\{tissue_slug\}_\{cells_str\}_\{download_mode\}_\{census_version\}\.h5ad"',
                r'return f"{tissue_slug}_{cells_str}_{download_mode}_{census_version}.h5ad"',
                cell.source
            )
            if new_source != cell.source:
                cell.source = new_source
                modified = True

if modified:
    with open(notebook_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print("Fixed cancer_cellxgene_integration.ipynb prefix")
