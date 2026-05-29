import nbformat as nbf
import os

ipynb_path = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/pan_cancer_meta_analysis.ipynb'
if os.path.exists(ipynb_path):
    with open(ipynb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    for cell in nb.cells:
        if cell.cell_type == 'code':
            if "import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX" in cell.source:
                cell.source = cell.source.replace(
                    "import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX",
                    "import sys\nif '..' not in sys.path: sys.path.append('..')\nfrom pan_cancer_config import ANALYSIS_SUFFIX"
                )
    
    with open(ipynb_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print("Updated pan_cancer_meta_analysis.ipynb")
