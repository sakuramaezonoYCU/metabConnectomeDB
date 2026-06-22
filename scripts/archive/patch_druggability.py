import os
import re

scripts_to_patch = [
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/compute_druggability_targets.py",
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/update_druggability_notebook.py"
]
notebooks_to_patch = [
    "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/druggability_axis_analysis.ipynb"
]

for py_path in scripts_to_patch:
    if not os.path.exists(py_path): continue
    try:
        with open(py_path, 'r', encoding='utf-8') as f:
            content = f.read()
        old_content = content
        
        # Add imports
        if "from pan_cancer_config" not in content:
            content = content.replace("import os", "import os\nimport sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX, get_de_csv_path\n")
            
        # Fix compute_druggability_targets.py reading logic
        if "f\"primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv\"" in content:
            content = content.replace("res_file = os.path.join(OUTPUT_DIR, f\"{cancer}_results\", f\"primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv\")", "res_file = get_de_csv_path(cancer)")

        # Fix string literals
        content = content.replace("'druggable_targets_23_genes.csv'", "f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv'")
        content = content.replace("'druggable_targets_181_genes.csv'", "f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv'")
        
        # Output printing fixes in compute_druggability_targets.py
        content = content.replace("the 23 pan-cancer genes", "the strictly conserved pan-cancer genes")
        content = content.replace("the 181 >=4 cancer genes", "the broadly conserved (>=4) cancer genes")
        
        if content != old_content:
            with open(py_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"Patched {os.path.basename(py_path)}")
    except Exception as e:
        print(e)
        
import nbformat as nbf
for nb_path in notebooks_to_patch:
    if not os.path.exists(nb_path): continue
    try:
        with open(nb_path, 'r', encoding='utf-8') as f:
            nb = nbf.read(f, as_version=4)
        modified = False
        for cell in nb.cells:
            old_source = cell.source
            if cell.cell_type == 'code':
                if "'druggable_targets_23_genes.csv'" in cell.source:
                    if "ANALYSIS_SUFFIX" not in cell.source:
                        cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                    cell.source = cell.source.replace("'druggable_targets_23_genes.csv'", "f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv'")
                if "'druggable_targets_181_genes.csv'" in cell.source:
                    if "ANALYSIS_SUFFIX" not in cell.source:
                        cell.source = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n\n" + cell.source
                    cell.source = cell.source.replace("'druggable_targets_181_genes.csv'", "f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv'")
            
            # Markdown scrub
            if cell.cell_type == 'markdown':
                cell.source = cell.source.replace("23 Strictly Conserved", "Strictly Conserved")
                cell.source = cell.source.replace("181 Broadly Conserved", "Broadly Conserved")
                cell.source = cell.source.replace("23 strictly conserved", "strictly conserved")
                cell.source = cell.source.replace("181 broadly conserved", "broadly conserved")
                cell.source = cell.source.replace("druggable_targets_23_genes", "druggable_targets_strictly_conserved")
                cell.source = cell.source.replace("druggable_targets_181_genes", "druggable_targets_broadly_conserved")
                
            if cell.source != old_source:
                modified = True
        if modified:
            with open(nb_path, 'w', encoding='utf-8') as f:
                nbf.write(nb, f)
            print(f"Patched {os.path.basename(nb_path)}")
    except Exception as e:
        print(e)
