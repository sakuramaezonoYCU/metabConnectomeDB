import os

base_dir = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts"
scripts_to_patch = [
    "compute_directional_ccc.py",
    "compute_permutation_null.py",
    "compute_stat3_network.py",
    "validate_tcga_signature.py",
    "simulate_oxygen_gradient.py"
]

for script in scripts_to_patch:
    path = os.path.join(base_dir, script)
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # We need to make sure the scripts use f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}"
        # They currently have "deepdive_conserved_metabGeneSig" hardcoded.
        content = content.replace('"deepdive_conserved_metabGeneSig"', 'f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}"')
        content = content.replace("'deepdive_conserved_metabGeneSig'", 'f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}"')
        
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Patched {script}")

# Now patch the notebook
import nbformat as nbf
nb_path = os.path.join(base_dir, "deepdive_conserved_metabGeneSig.ipynb")
if os.path.exists(nb_path):
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    modified = False
    for cell in nb.cells:
        if cell.cell_type == 'code':
            old_source = cell.source
            if 'output_base="deepdive_conserved_metabGeneSig"' in cell.source:
                cell.source = cell.source.replace('output_base="deepdive_conserved_metabGeneSig"', 'output_base=f"deepdive_conserved_metabGeneSig{ANALYSIS_SUFFIX}"')
                modified = True
            
            # Clear outputs so the notebook is fresh
            cell.outputs = []
            
    if modified:
        with open(nb_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched deepdive_conserved_metabGeneSig.ipynb")
