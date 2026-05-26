import nbformat as nbf

notebook_path = 'scripts/primary_vs_metastasis_comparison.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

modified = False
for cell in nb.cells:
    if cell.cell_type == 'code':
        source = cell.source
        
        # We find the cell with Main DE
        if "sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')" in source:
            replacement = """try:
    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
except Exception as e:
    print(f"Warning: Cannot perform Main DE. Check if dataset contains both Primary and Metastasis cells for these PRIMARY_TISSUES. Error: {e}")
    # Create empty dummy results to prevent subsequent cells from crashing
    import numpy as np
    adata.uns['rank_genes_groups'] = {
        'names': np.rec.fromarrays([np.array([], dtype=object)], names=['Metastasis']),
        'scores': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
        'logfoldchanges': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
        'pvals': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
        'pvals_adj': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis'])
    }"""
            
            new_source = source.replace("sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')", replacement)
            cell.source = new_source
            modified = True

if modified:
    with open(notebook_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print("Notebook successfully patched for missing group errors.")
else:
    print("Code pattern not found. Nothing modified.")
