import os
import glob
import re

scripts_dir = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts"

def patch_python_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()

    # If the file hardcodes pan_cancer_23_genes.csv without using ANALYSIS_SUFFIX,
    # we might need to make sure ANALYSIS_SUFFIX is imported
    # But since many files do different things, let's just do targeted replaces.
    
    modified = False

    # Replace hardcoded 5MetCan_100k assignments with an import
    if "import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX" in content or 'import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX' in content:
        content = re.sub(r"ANALYSIS_SUFFIX\s*=\s*['\"]_5MetCan_100k['\"]", 
                         "import sys\\nif '..' not in sys.path: sys.path.append('..')\\nfrom pan_cancer_config import ANALYSIS_SUFFIX", content)
        modified = True

    if f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv" in content:
        # replace literal strings with f-strings using ANALYSIS_SUFFIX
        # e.g. f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv" -> f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv"
        # We also need to ensure ANALYSIS_SUFFIX is imported
        if "from pan_cancer_config import ANALYSIS_SUFFIX" not in content:
            content = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n" + content
        content = content.replace("f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv'", "f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv'")
        content = content.replace('f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv"', 'f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv"')
        modified = True

    if f"tissue_specific_correlation{ANALYSIS_SUFFIX}.png" in content:
        if "from pan_cancer_config import ANALYSIS_SUFFIX" not in content:
            content = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n" + content
        content = content.replace("f'tissue_specific_correlation{ANALYSIS_SUFFIX}.png'", "f'tissue_specific_correlation{ANALYSIS_SUFFIX}.png'")
        content = content.replace('f"tissue_specific_correlation{ANALYSIS_SUFFIX}.png"', 'f"tissue_specific_correlation{ANALYSIS_SUFFIX}.png"')
        modified = True

    if f"serotonin_proximity_results{ANALYSIS_SUFFIX}.csv" in content:
        if "from pan_cancer_config import ANALYSIS_SUFFIX" not in content:
            content = "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX\n" + content
        content = content.replace('f"serotonin_proximity_results{ANALYSIS_SUFFIX}.csv"', 'f"serotonin_proximity_results{ANALYSIS_SUFFIX}.csv"')
        modified = True

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Patched Python file: {os.path.basename(filepath)}")

def patch_ipynb_file(filepath):
    import nbformat as nbf
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            nb = nbf.read(f, as_version=4)
    except Exception:
        return

    modified = False
    for cell in nb.cells:
        if cell.cell_type == 'code':
            if "import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX" in cell.source or 'import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX' in cell.source:
                cell.source = re.sub(r"ANALYSIS_SUFFIX\s*=\s*['\"]_5MetCan_100k['\"]", 
                         "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX", cell.source)
                modified = True

            # If there's hardcoded output_base like 'nr1d2_master_regulator_analysis_5MetCan_100k'
            if "_5MetCan_100k" in cell.source and "output_base" in cell.source:
                cell.source = re.sub(r"(output_base\s*=\s*['\"][a-zA-Z0-9_]+?)_5MetCan_100k(['\"])", r"\1' + ANALYSIS_SUFFIX", cell.source)
                modified = True

        if cell.cell_type == 'markdown':
            if f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv" in cell.source:
                cell.source = cell.source.replace(f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv", "pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv")
                modified = True

    if modified:
        with open(filepath, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print(f"Patched Notebook: {os.path.basename(filepath)}")

for root, _, files in os.walk(scripts_dir):
    for file in files:
        if file == "pan_cancer_config.py" or file == "rename_outputs.py":
            continue
        filepath = os.path.join(root, file)
        if file.endswith('.py'):
            patch_python_file(filepath)
        elif file.endswith('.ipynb'):
            patch_ipynb_file(filepath)
