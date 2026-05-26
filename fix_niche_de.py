import nbformat as nbf
import re

pvsm_nb = "scripts/primary_vs_metastasis_comparison.ipynb"
with open(pvsm_nb, "r") as f:
    nb = nbf.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == "code" and "sc.tl.rank_genes_groups(adata_site," in cell.source:
        if "if len(adata_site.obs['comparison_group'].unique()) < 2:" not in cell.source:
            # We want to wrap the rank_genes_groups call
            replacement = """        if len(adata_site.obs['comparison_group'].unique()) < 2:
            print(f"Warning: Cannot perform Niche DE for {site}. Missing Primary or Metastasis cells.")
            import numpy as np
            adata_site.uns['rank_genes_groups'] = {
                'names': np.rec.fromarrays([np.array([], dtype=object)], names=[site]),
                'scores': np.rec.fromarrays([np.array([], dtype=float)], names=[site]),
                'logfoldchanges': np.rec.fromarrays([np.array([], dtype=float)], names=[site]),
                'pvals': np.rec.fromarrays([np.array([], dtype=float)], names=[site]),
                'pvals_adj': np.rec.fromarrays([np.array([], dtype=float)], names=[site])
            }
        else:
            try:
                sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')
            except Exception as e:
                print(f"Warning: Niche DE failed for {site}. Error: {e}")
                import numpy as np
                adata_site.uns['rank_genes_groups'] = {
                    'names': np.rec.fromarrays([np.array([], dtype=object)], names=[site]),
                    'scores': np.rec.fromarrays([np.array([], dtype=float)], names=[site]),
                    'logfoldchanges': np.rec.fromarrays([np.array([], dtype=float)], names=[site]),
                    'pvals': np.rec.fromarrays([np.array([], dtype=float)], names=[site]),
                    'pvals_adj': np.rec.fromarrays([np.array([], dtype=float)], names=[site])
                }"""
            
            # regex replace the old try/except
            cell.source = re.sub(
                r'        try:\n\s+sc\.tl\.rank_genes_groups\(adata_site.*?(?=\n        # Extract)',
                replacement + "\n",
                cell.source,
                flags=re.DOTALL
            )
            
            # fallback
            if "if len(adata_site.obs['comparison_group'].unique()) < 2:" not in cell.source:
                print("Fallback replace for Niche DE")
                cell.source = cell.source.replace(
                    "sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')",
                    "if len(adata_site.obs['comparison_group'].unique()) < 2:\n            pass\n        else:\n            sc.tl.rank_genes_groups(adata_site, groupby='comparison_group', groups=[site], reference='Primary', method='wilcoxon')"
                )

with open(pvsm_nb, "w") as f:
    nbf.write(nb, f)
print("Fixed Niche DE logic.")
