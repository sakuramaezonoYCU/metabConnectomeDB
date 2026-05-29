import os
import re

files_to_patch = [
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/generate_combined_pan_cancer_notebook.py",
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/generate_pan_cancer_notebook.py"
]

for filepath in files_to_patch:
    if not os.path.exists(filepath):
        continue
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
        
    old_content = content
    
    # 1. Add import ANALYSIS_SUFFIX
    if "from pan_cancer_config import ANALYSIS_SUFFIX" not in content:
        content = content.replace("import os", "import os\nimport sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n")

    # 2. Fix UpSet plot
    content = content.replace("'upset_plot.png'", "f'upset_plot{ANALYSIS_SUFFIX}.png'")
    content = content.replace("'../output/pan_cancer_meta_results/upset_plot.png'", "f'../output/pan_cancer_meta_results/upset_plot{ANALYSIS_SUFFIX}.png'")
    
    # 3. Fix Network plot
    content = content.replace("'metabolite_target_network.png'", "f'metabolite_target_network{ANALYSIS_SUFFIX}.png'")
    content = content.replace("'../output/pan_cancer_meta_results/metabolite_target_network.png'", "f'../output/pan_cancer_meta_results/metabolite_target_network{ANALYSIS_SUFFIX}.png'")
    
    # 4. Fix Signature Scores
    cancers = ['breast', 'lung', 'colorectal', 'melanoma', 'ovarian']
    for c in cancers:
        content = content.replace(f"'{c}_primary_signature_score.png'", f"f'{c}_primary_signature_score{{ANALYSIS_SUFFIX}}.png'")
        content = content.replace(f"'{c}_primary_signature_scores.csv'", f"f'{c}_primary_signature_scores{{ANALYSIS_SUFFIX}}.csv'")

    if content != old_content:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Patched generator script: {os.path.basename(filepath)}")
    else:
        print(f"No changes needed for: {os.path.basename(filepath)}")

# Also let's just clear the outputs in the current pan_cancer_meta_analysis.ipynb so it doesn't show old errors
import nbformat as nbf
nb_path = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/pan_cancer_meta_analysis.ipynb"
try:
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
    for cell in nb.cells:
        if cell.cell_type == 'code':
            cell.outputs = []
    with open(nb_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print("Cleared stale outputs in pan_cancer_meta_analysis.ipynb")
except Exception as e:
    print(e)
