import nbformat as nbf

def fix_notebook(nb_path):
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    for cell in nb.cells:
        if cell.cell_type == 'code' and 'nbconvert' in cell.source and '--to html' in cell.source:
            if 'if globals().get("SAVE_AS_HTML", True):' not in cell.source:
                cell.source = 'if globals().get("SAVE_AS_HTML", True):\n' + '\n'.join(['    ' + line for line in cell.source.split('\n')])
                
        if cell.cell_type == 'code' and "output_csv = os.path.join(OUTPUT_DIR, 'primary_vs_metastasis_DE_metabolic_targets.csv')" in cell.source:
            # Fix hardcoded output_csv
            cell.source = cell.source.replace("'primary_vs_metastasis_DE_metabolic_targets.csv'", "f'primary_vs_metastasis_{CANCER_TYPE_NAME}_DE_metabolic_targets.csv'")
            
    with open(nb_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)

fix_notebook('scripts/primary_vs_metastasis_comparison.ipynb')
fix_notebook('scripts/orphan_metabolic_immune_evasion.ipynb')
print("Notebooks fixed!")
