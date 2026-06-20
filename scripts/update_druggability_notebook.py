import nbformat as nbf
import os
import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX, get_de_csv_path

import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DRUGGABILITY_DIR = os.path.join(BASE_DIR, 'output', 'druggability')

def update_notebook():
    nb_path = 'druggability_axis_analysis.ipynb'
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    # Find the export cell (usually at the very end)
    export_cell_index = -1
    for i, cell in enumerate(nb.cells):
        src = "".join(cell.get("source", []))
        if "SAVE_AS_HTML" in src and "nbconvert" in src:
            export_cell_index = i
            break
            
    # New cells to add at the end
    new_cells = []
    
    # Also add a config print cell at the very beginning
    print_config_code = """
print('--- INJECTED PIPELINE CONFIGURATION ---')
from pan_cancer_config import CANCER_CAP, CANCERS_TO_RUN, ANALYSIS_SUFFIX, DGIDB_API_URL
print(f'CANCERS_TO_RUN: {CANCERS_TO_RUN}')
print(f'CANCER_CAP: {CANCER_CAP}')
print(f'ANALYSIS_SUFFIX: {ANALYSIS_SUFFIX}')
print(f'DGIDB_API_URL: {DGIDB_API_URL}')
"""
    nb.cells.insert(0, nbf.v4.new_code_cell(print_config_code))
    
    # Markdown
    new_cells.append(nbf.v4.new_markdown_cell("""---
## Druggability of Pan-Cancer Conserved Genes
In addition to the specific GLS axis, we also query the DGIdb database for the strictly conserved **pan-cancer genes** (upregulated in all cancers) and the broadly conserved **genes** (upregulated in $\ge$ max-1 cancers).
"""))

    # Code for strictly conserved genes
    code_strict = """df_strict = pd.read_csv(os.path.join(OUTPUT_DIR, f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv'))
print(f"Total drug interactions for strictly conserved genes: {len(df_strict)}")
display(df_strict.head(10))"""
    cell_strict = nbf.v4.new_code_cell(code_strict)
    
    csv_strict_path = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(csv_strict_path):
        df_strict = pd.read_csv(csv_strict_path)
        output_txt = f"Total drug interactions for strictly conserved genes: {len(df_strict)}\n"
        html_tbl = df_strict.head(10).to_html()
        cell_strict.outputs.append(nbf.v4.new_output("stream", name="stdout", text=output_txt))
        cell_strict.outputs.append(nbf.v4.new_output("execute_result", data={"text/html": html_tbl}, execution_count=2))
        
    new_cells.append(cell_strict)
    
    # Code for broadly conserved genes
    code_broad = """df_broad = pd.read_csv(os.path.join(OUTPUT_DIR, f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv'))
print(f"Total drug interactions for broadly conserved genes: {len(df_broad)}")
display(df_broad.head(10))"""
    cell_broad = nbf.v4.new_code_cell(code_broad)
    
    csv_broad_path = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(csv_broad_path):
        df_broad = pd.read_csv(csv_broad_path)
        output_txt = f"Total drug interactions for broadly conserved genes: {len(df_broad)}\n"
        html_tbl = df_broad.head(10).to_html()
        cell_broad.outputs.append(nbf.v4.new_output("stream", name="stdout", text=output_txt))
        cell_broad.outputs.append(nbf.v4.new_output("execute_result", data={"text/html": html_tbl}, execution_count=3))
        
    new_cells.append(cell_broad)

    if export_cell_index != -1:
        # Insert before the export cell
        nb.cells[export_cell_index:export_cell_index] = new_cells
    else:
        nb.cells.extend(new_cells)
        
    with open(nb_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
        
    print("Notebook updated successfully.")

if __name__ == '__main__':
    update_notebook()
