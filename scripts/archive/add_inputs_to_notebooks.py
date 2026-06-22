import os
import glob
import nbformat
import re

def extract_inputs(nb):
    inputs = set()
    configs = set()
    for cell in nb.cells:
        if cell.cell_type == 'code':
            source = cell.source
            # find strings ending in .csv, .h5ad, .json, .txt, .tsv
            matches = re.findall(r'[\'"]([^\'"]+\.(?:csv|h5ad|json|txt|tsv|gz))[\'"]', source)
            for m in matches:
                if 'input/' in m or 'output/' in m or m.startswith('.'):
                    inputs.add(m)
            # find config imports
            config_matches = re.findall(r'from\s+([a-zA-Z0-9_]+_config)\s+import', source)
            for c in config_matches:
                configs.add(c + '.py')
    return list(inputs), list(configs)

def inject_inputs_cell(nb, inputs, configs):
    # Check if already has Inputs section
    for cell in nb.cells:
        if cell.cell_type == 'markdown' and '### Inputs' in cell.source:
            return False
            
    if not inputs and not configs:
        return False
        
    md_content = "### Inputs / Parameters\n*Explicitly documented for traceability and reproducibility.*\n"
    if inputs:
        md_content += "\n**Input Data Files:**\n"
        for i in sorted(inputs):
            md_content += f"- `{i}`\n"
            
    if configs:
        md_content += "\n**Configuration Files:**\n"
        for c in sorted(configs):
            md_content += f"- `{c}`\n"
            
    # Insert after the first markdown cell if it exists
    insert_idx = 0
    if len(nb.cells) > 0 and nb.cells[0].cell_type == 'markdown':
        insert_idx = 1
        
    new_cell = nbformat.v4.new_markdown_cell(md_content)
    nb.cells.insert(insert_idx, new_cell)
    return True

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    notebooks = glob.glob(os.path.join(base_dir, "*.ipynb"))
    
    for nb_path in notebooks:
        with open(nb_path, 'r', encoding='utf-8') as f:
            nb = nbformat.read(f, as_version=4)
            
        inputs, configs = extract_inputs(nb)
        
        if inject_inputs_cell(nb, inputs, configs):
            with open(nb_path, 'w', encoding='utf-8') as f:
                nbformat.write(nb, f)
            print(f"Injected Inputs/Parameters into {os.path.basename(nb_path)}")
        else:
            print(f"Skipped {os.path.basename(nb_path)} (no inputs found or already documented)")

if __name__ == '__main__':
    main()
