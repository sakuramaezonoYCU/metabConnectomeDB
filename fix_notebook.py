import json

with open('scripts/primary_vs_metastasis_comparison.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        new_source = []
        for line in cell['source']:
            # replace the bad code
            if line == "            if len(adata_site.obs['comparison_group'].unique()) < 2:\n":
                pass
            elif line == "            pass\n":
                pass
            elif line == "        else:\n":
                pass
            else:
                new_source.append(line)
                
        source_str = "".join(cell['source'])
        bad_str = """        try:
            if len(adata_site.obs['comparison_group'].unique()) < 2:
            pass
        else:
            sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')"""
        
        good_str = """        try:
            if len(adata_site.obs['comparison_group'].unique()) < 2:
                print(f"Warning: Not enough groups for site {site}")
                continue
            sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')"""
            
        if bad_str in source_str:
            new_source_str = source_str.replace(bad_str, good_str)
            
            bad_str2 = """    try:
        if len(adata.obs['site'].unique()) < 2:
        pass
    else:
        sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')"""
        
            good_str2 = """    try:
        if len(adata.obs['site'].unique()) < 2:
            print(f"Warning: Cannot perform Main DE. Check if dataset contains both Primary and Metastasis cells for these PRIMARY_TISSUES. Error: Not enough groups.")
            target_genes_de = pd.DataFrame()
        else:
            sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')"""
            
            if bad_str2 in new_source_str:
                new_source_str = new_source_str.replace(bad_str2, good_str2)
                
            cell['source'] = [s + '\n' for s in new_source_str.split('\n')]
            if not new_source_str.endswith('\n'):
                cell['source'][-1] = cell['source'][-1][:-1]
            elif cell['source']:
                cell['source'].pop()

with open('scripts/primary_vs_metastasis_comparison.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)
