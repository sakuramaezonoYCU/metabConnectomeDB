import json
import os

notebook_path = 'scripts/orphan_metabolic_immune_evasion.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

modified = False
for idx, cell in enumerate(nb.get('cells', [])):
    if cell.get('cell_type') == 'code':
        source = "".join(cell.get('source', []))
        
        # 1. Modify cell 1's output_dir assignment to check globals()
        if "output_dir = os.path.join(workspace_dir, 'output')" in source:
            print(f"Modifying cell {idx} for output_dir")
            source = source.replace(
                "output_dir = os.path.join(workspace_dir, 'output')",
                "output_dir = globals().get('OUTPUT_DIR', os.path.join(workspace_dir, 'output'))"
            )
            cell['source'] = [source]
            modified = True
            
        # 2. Modify cell 3's HTML glob, cancer_name resolver, and fallback block
        if "html_files = glob.glob(os.path.join(output_dir, 'cancer_*.html'))" in source:
            print(f"Modifying cell {idx} for html_files, cancer_name, and fallback block")
            source = source.replace(
                "html_files = glob.glob(os.path.join(output_dir, 'cancer_*.html'))",
                "html_files = glob.glob(os.path.join(output_dir, '*.html'))"
            )
            source = source.replace(
                "cancer_name = os.path.basename(f).replace('cancer_', '').split('_')[0]",
                "cancer_name = globals().get('CANCER_TYPE_NAME', os.path.basename(f)).replace('cancer_', '').split('_')[0]\n    cancer_name = cancer_name.replace('_results', '').replace('_cancer', '')"
            )
            source = source.replace(
                "'Cancer': ['lung'] * len(fallback_genes)",
                "'Cancer': [globals().get('CANCER_TYPE_NAME', 'breast_cancer').replace('_results', '').replace('_cancer', '')] * len(fallback_genes)"
            )
            cell['source'] = [source]
            modified = True

if modified:
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(nb, f, indent=1)
    print("🎉 Evasion notebook updated successfully!")
else:
    print("⚠️ No cells matched the modification target text.")
