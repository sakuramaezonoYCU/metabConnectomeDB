import os
import glob

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
        
        # Replace the hardcoded directory name
        content = content.replace("deepdive_23_metabGeneSig", "deepdive_conserved_metabGeneSig")
        # Replace the stat3 output file name
        content = content.replace("stat3_u87_targets_23genes.csv", "stat3_u87_targets_strictly_conserved.csv")
        # Replace generic "23 genes" text
        content = content.replace("23 pan-cancer conserved genes", "strictly conserved pan-cancer genes")
        content = content.replace("23-gene", "strictly-conserved")
        content = content.replace("23-gene metabolic signature", "strictly conserved metabolic signature")
        content = content.replace("23 pan-cancer genes", "strictly conserved pan-cancer genes")
        
        # In simulate_oxygen_gradient.py: pan_cancer_23
        content = content.replace("pan_cancer_23", "pan_cancer_genes")
        
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Patched {script}")
