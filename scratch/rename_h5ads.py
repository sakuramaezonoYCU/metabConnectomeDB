import os
import re

cancer_output_dir = 'output'

mappings = [
    {
        'folder': 'breast_results',
        'old': 'breast_mammary-gland_liver_axilla_chest-wall_500k_whole_transcriptome_2025-11-08.h5ad',
        'new_disease_slug': 'breast-cancer'
    },
    {
        'folder': 'lung_results',
        'old': 'lung_lymph-node_brain_pleural-fluid_500k_whole_transcriptome_2025-11-08.h5ad',
        'new_disease_slug': 'lung-cancer_lung-adenocarcinoma'
    },
    {
        'folder': 'ovarian_results',
        'old': 'ovary_abdomen_omentum_uterus_500k_whole_transcriptome_2025-11-08.h5ad',
        'new_disease_slug': 'ovarian-cancer_malignant-ovarian-serous-tumor'
    }
]

for m in mappings:
    old_path = os.path.join(cancer_output_dir, m['folder'], m['old'])
    new_name = f"{m['new_disease_slug']}_{m['old']}"
    new_path = os.path.join(cancer_output_dir, m['folder'], new_name)
    
    if os.path.exists(old_path):
        os.rename(old_path, new_path)
        print(f"Renamed: {m['old']} -> {new_name}")
    elif os.path.exists(new_path):
        print(f"Already named: {new_name}")
    else:
        print(f"NOT FOUND: {old_path}")
