import json
import os

ipynb_file = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/pan_cancer_meta_analysis.ipynb'

with open(ipynb_file, 'r', encoding='utf-8') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = cell['source']
        for i, line in enumerate(source):
            if "We have relaxed the strict threshold to generate stratified combinations of 4 out of 5 cancers" in line:
                source[i] = '    "print(\\"METHODOLOGY NOTE (MAX CANCER - 1 Rule): Because the strict intersection across all 5 cancers yielded 0 genes, the pipeline automatically falls back to utilizing the union of 4-cancer combinations to ensure a robust meta-signature is evaluated downstream.\\")\\n"\n'

with open(ipynb_file, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Updated pan_cancer_meta_analysis.ipynb with MAX CANCER - 1 rule.")
