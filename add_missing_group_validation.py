import nbformat as nbf

notebook_path = 'scripts/primary_vs_metastasis_comparison.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

modified = False
for cell in nb.cells:
    if cell.cell_type == 'code':
        source = cell.source
        
        # Check Main DE
        if "sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')" in source:
            new_source = source.replace(
                "sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')",
                "if 'Metastasis' not in adata.obs['site'].values or 'Primary' not in adata.obs['site'].values:\n    print('Warning: Cannot perform Main DE. Dataset must contain both Primary and Metastasis cells.')\nelse:\n    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')"
            )
            # We also need to indent the subsequent code that depends on the DE results, OR just let it throw an error later, OR we can check if DE was run.
            # Actually, it's safer to just wrap it in a try-except, because the rest of the notebook depends on it.
            # Wait, if we wrap it, the next cell will fail when trying to access `adata.uns['rank_genes_groups']`.
            # Let's wrap the entire DE execution block in a try-except so we see a clean error, or just do the if-else and then the next cells will fail cleanly.
            modified = True
            
        # Or better yet, we can fix the cell by wrapping the DE call in try-except and handling the IndexError cleanly.
        
with open(notebook_path, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)
