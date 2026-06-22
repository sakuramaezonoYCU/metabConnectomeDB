import nbformat as nbf
import os
import re

notebooks_to_patch = [
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/pan_cancer_meta_analysis.ipynb"
]

scripts_to_patch = [
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/generate_combined_pan_cancer_notebook.py",
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/generate_pan_cancer_notebook.py"
]

replacements = [
    ("The 23 Pan-Cancer Genes", "The Conserved Pan-Cancer Genes"),
    ("the 23 pan-cancer conserved genes", "the pan-cancer conserved genes"),
    ("The 23 strictly conserved genes", "The strictly conserved genes"),
    ("the 23 strictly conserved genes", "the strictly conserved genes"),
    ("23-Gene Signature", "Conserved Gene Signature"),
    ("23-gene signature", "conserved gene signature"),
    ("the 23 genes", "the conserved genes"),
    ("Pan-Cancer 23-Gene Conserved", "Pan-Cancer Conserved"),
]

# Patch the notebooks
for nb_path in notebooks_to_patch:
    if not os.path.exists(nb_path): continue
    try:
        with open(nb_path, 'r', encoding='utf-8') as f:
            nb = nbf.read(f, as_version=4)
        modified = False
        for cell in nb.cells:
            old_source = cell.source
            for old_str, new_str in replacements:
                cell.source = cell.source.replace(old_str, new_str)
            if cell.source != old_source:
                modified = True
        if modified:
            with open(nb_path, 'w', encoding='utf-8') as f:
                nbf.write(nb, f)
            print(f"Removed hardcoded 23s from {os.path.basename(nb_path)}")
    except Exception as e:
        print(e)

# Patch the generator scripts
for py_path in scripts_to_patch:
    if not os.path.exists(py_path): continue
    try:
        with open(py_path, 'r', encoding='utf-8') as f:
            content = f.read()
        old_content = content
        for old_str, new_str in replacements:
            content = content.replace(old_str, new_str)
        if content != old_content:
            with open(py_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"Removed hardcoded 23s from {os.path.basename(py_path)}")
    except Exception as e:
        print(e)
