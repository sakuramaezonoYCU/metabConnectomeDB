import os
import nbformat as nbf
import glob
import shutil

# 1. Fix primary_vs_metastasis_comparison.ipynb
notebook_path1 = 'scripts/primary_vs_metastasis_comparison.ipynb'
with open(notebook_path1, 'r', encoding='utf-8') as f:
    nb1 = nbf.read(f, as_version=4)

modified1 = False
for cell in nb1.cells:
    if cell.cell_type == 'code':
        if "OUTPUT_DIR = os.path.join(BASE_DIR, 'output')" in cell.source:
            cell.source = cell.source.replace(
                "OUTPUT_DIR = os.path.join(BASE_DIR, 'output')",
                "OUTPUT_DIR = globals().get('OUTPUT_DIR', os.path.join(BASE_DIR, 'output'))"
            )
            modified1 = True

if modified1:
    with open(notebook_path1, 'w', encoding='utf-8') as f:
        nbf.write(nb1, f)
    print("Fixed OUTPUT_DIR in primary_vs_metastasis_comparison.ipynb")

# 2. Fix orphan_metabolic_immune_evasion.ipynb
notebook_path2 = 'scripts/orphan_metabolic_immune_evasion.ipynb'
with open(notebook_path2, 'r', encoding='utf-8') as f:
    nb2 = nbf.read(f, as_version=4)

modified2 = False
for cell in nb2.cells:
    if cell.cell_type == 'code':
        if "output_dir = os.path.join(workspace_dir, 'output')" in cell.source:
            cell.source = cell.source.replace(
                "output_dir = os.path.join(workspace_dir, 'output')",
                "output_dir = globals().get('OUTPUT_DIR', os.path.join(workspace_dir, 'output'))"
            )
            modified2 = True

if modified2:
    with open(notebook_path2, 'w', encoding='utf-8') as f:
        nbf.write(nb2, f)
    print("Fixed output_dir in orphan_metabolic_immune_evasion.ipynb")

# 3. Clean up stray files in output/
base_output_dir = 'output'
breast_dir = 'output/breast_cancer'

os.makedirs(breast_dir, exist_ok=True)

# Move breast cancer specific csv files that landed in output/
csvs_to_move = [
    'primary_vs_metastasis_breast_cancer_DE_metabolic_targets.csv'
]
for csv_file in csvs_to_move:
    src = os.path.join(base_output_dir, csv_file)
    dst = os.path.join(breast_dir, csv_file)
    if os.path.exists(src):
        shutil.move(src, dst)
        print(f"Moved {csv_file} to breast_cancer/")

# Delete stray HTMLs and PNGs that are duplicates or orphans
files_to_delete = [
    'immune_evasion_orphan_metabolic_candidates.html',
    'primary_vs_metastasis_DE_metabolic_targets.html',
    'primary_vs_metastasis.html',
    'immune_cell_specificity_dotplot.png',
    'orphan_metabolic_immune_connectome_network.png',
    'robust_orphan_interactions_bar.png',
    'primary_vs_metastasis_DE_metabolic_targets_immune_evasion_orphan_metabolic_candidates.csv',
    'immune_evasion_orphan_metabolic_candidates.csv',
    'immune_evasion_dark_matter_candidates.csv',
    'immune_evasion_orphan_metabolic_candidates_immune_evasion_orphan_metabolic_candidates.csv'
]

for f in files_to_delete:
    p = os.path.join(base_output_dir, f)
    if os.path.exists(p):
        os.remove(p)
        print(f"Deleted stray file: {f}")

print("Cleanup complete!")
