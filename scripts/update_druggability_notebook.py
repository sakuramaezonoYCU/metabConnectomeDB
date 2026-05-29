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
            
    # New cells to add
    new_cells = []
    
    # Markdown
    new_cells.append(nbf.v4.new_markdown_cell("""---
## Druggability of Pan-Cancer Conserved Genes
In addition to the specific GLS axis, we also query the DGIdb database for the strictly conserved **23 pan-cancer genes** (upregulated in all 5 cancers) and the broadly conserved **181 genes** (upregulated in $\ge$ 4 cancers).
"""))

    # Code for 23 genes
    code_23 = """df_23 = pd.read_csv(os.path.join(OUTPUT_DIR, f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv'))
print(f"Total drug interactions for 23 conserved genes: {len(df_23)}")
display(df_23.head(10))"""
    cell_23 = nbf.v4.new_code_cell(code_23)
    
    csv_23_path = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_strictly_conserved{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(csv_23_path):
        df_23 = pd.read_csv(csv_23_path)
        output_txt = f"Total drug interactions for 23 conserved genes: {len(df_23)}\n"
        html_tbl = df_23.head(10).to_html()
        cell_23.outputs.append(nbf.v4.new_output("stream", name="stdout", text=output_txt))
        cell_23.outputs.append(nbf.v4.new_output("execute_result", data={"text/html": html_tbl}, execution_count=2))
        
    new_cells.append(cell_23)
    
    # Code for 181 genes
    code_181 = """df_181 = pd.read_csv(os.path.join(OUTPUT_DIR, f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv'))
print(f"Total drug interactions for 181 conserved genes: {len(df_181)}")
display(df_181.head(10))"""
    cell_181 = nbf.v4.new_code_cell(code_181)
    
    csv_181_path = os.path.join(DRUGGABILITY_DIR, f'druggable_targets_broadly_conserved{ANALYSIS_SUFFIX}.csv')
    if os.path.exists(csv_181_path):
        df_181 = pd.read_csv(csv_181_path)
        output_txt = f"Total drug interactions for 181 conserved genes: {len(df_181)}\n"
        html_tbl = df_181.head(10).to_html()
        cell_181.outputs.append(nbf.v4.new_output("stream", name="stdout", text=output_txt))
        cell_181.outputs.append(nbf.v4.new_output("execute_result", data={"text/html": html_tbl}, execution_count=3))
        
    new_cells.append(cell_181)

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
