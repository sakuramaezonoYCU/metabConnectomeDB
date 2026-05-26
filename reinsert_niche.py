import nbformat as nbf
import os

notebook_path = 'scripts/primary_vs_metastasis_comparison.ipynb'
with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

code_to_insert = """# Check if we have multiple metastatic sites
meta_sites = [t for t in adata[adata.obs['site'] == 'Metastasis'].obs['tissue_general'].unique() if pd.notnull(t)]

if len(meta_sites) > 1:
    print(f"Found multiple metastatic niches: {meta_sites}\\nRunning site-specific DE...")
    
    # We need UpSetPlot for visualization
    try:
        from upsetplot import from_contents, plot
    except ImportError:
        import subprocess
        print("Installing upsetplot...")
        subprocess.check_call(['pip', 'install', 'upsetplot'])
        from upsetplot import from_contents, plot
        
    site_significant_genes = {}
    
    for site in meta_sites:
        # Create a subset with just the primary cells and THIS specific metastatic site
        adata_site = adata[(adata.obs['tissue_general'] == site) | (adata.obs['site'] == 'Primary')].copy()
        adata_site.obs['comparison_group'] = adata_site.obs.apply(lambda x: site if x['site'] == 'Metastasis' else 'Primary', axis=1)
        
        try:
            sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')
            result = adata_site.uns['rank_genes_groups']
            
            # Filter for significant metabolic targets
            site_df = pd.DataFrame({
                'names': result['names'][site],
                'logfoldchanges': result['logfoldchanges'][site],
                'pvals_adj': result['pvals_adj'][site]
            })
            
            # Significant up-regulated in metastasis
            sig_up = site_df[(site_df['names'].isin(target_genes)) & 
                             (site_df['pvals_adj'] < 0.05) & 
                             (site_df['logfoldchanges'] > 0.5)]['names'].tolist()
            
            if len(sig_up) > 0:
                site_significant_genes[site] = sig_up
                print(f" - {site}: {len(sig_up)} upregulated metabolic targets")
            else:
                print(f" - {site}: 0 upregulated metabolic targets found. Skipping for plot.")
            
        except Exception as e:
            print(f" - {site}: Could not run DE (maybe not enough cells). Error: {e}")
            
    if len(site_significant_genes) > 0:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(10, 6))
        upset_data = from_contents(site_significant_genes)
        axes = plot(upset_data, show_counts=True, sort_by='cardinality', fig=fig)
        plt.suptitle(f"Metastatic Convergence: Upregulated Targets across Niches", fontsize=14, y=1.05)
        plt.show()
        
        # Output the intersection (the pan-metastatic signature)
        if len(site_significant_genes) >= 2:
            common_genes = set.intersection(*[set(g) for g in site_significant_genes.values()])
            print(f"\\n🔥 Pan-Metastatic Conserved Targets ({len(common_genes)} genes):")
            print(list(common_genes))
else:
    print(f"Only one or zero metastatic sites found ({meta_sites}). Site-specific comparison skipped.")
"""

# Find the cell that says "# upsetplot removed due to library bugs" and replace it
found = False
for cell in nb.cells:
    if cell.cell_type == 'code' and 'upsetplot removed due to library bugs' in cell.source:
        cell.source = code_to_insert
        found = True

if found:
    with open(notebook_path, 'w', encoding='utf-8') as f:
        nbf.write(nb, f)
    print("Successfully re-inserted Niche Comparison code!")
else:
    print("Could not find the target cell to replace.")
