import os
import nbformat as nbf
import glob

# 1. Modify the notebook to remove the prefix
notebook_path = 'scripts/cancer_cellxgene_integration.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

modified = False
for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'return f"cancer_{disease_slug}_{tissue_slug}_{cells_str}_{download_mode}_{census_version}.h5ad"' in cell.source:
            cell.source = cell.source.replace(
                'return f"cancer_{disease_slug}_{tissue_slug}_{cells_str}_{download_mode}_{census_version}.h5ad"',
                'return f"{tissue_slug}_{cells_str}_{download_mode}_{census_version}.h5ad"'
            )
            modified = True

if modified:
    with open(notebook_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print("Updated notebook h5ad cache filename logic.")
else:
    print("Notebook cell not found or already modified.")

# 2. Rename existing files in output/breast_cancer/ to remove prefix
prefix = "cancer_breast-cancer_"
directory = "output/breast_cancer/"

if os.path.exists(directory):
    for filename in os.listdir(directory):
        if filename.startswith(prefix):
            old_path = os.path.join(directory, filename)
            new_filename = filename[len(prefix):]
            new_path = os.path.join(directory, new_filename)
            os.rename(old_path, new_path)
            print(f"Renamed: {filename} -> {new_filename}")
else:
    print(f"Directory {directory} not found.")

print("All done!")
