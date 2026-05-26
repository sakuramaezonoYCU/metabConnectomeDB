import os
import nbformat as nbf
import re
import shutil

# 1. Update run_cancer_pipeline.py
pipeline_file = "scripts/run_cancer_pipeline.py"
with open(pipeline_file, "r") as f:
    content = f.read()

# Change cancer_name_safe to use _results
old_line = 'cancer_name_safe = disease_filter_str.replace(" ", "_").replace(",", "")'
new_line = 'cancer_name_safe = f"{disease_filter_str.split()[0].lower()}_results"'
content = content.replace(old_line, new_line)

with open(pipeline_file, "w") as f:
    f.write(content)
print("Fixed run_cancer_pipeline.py to use <cancer>_results directory.")

# 2. Update cancer_cellxgene_integration.ipynb to remove prefix
cellxgene_nb = "scripts/cancer_cellxgene_integration.ipynb"
with open(cellxgene_nb, "r") as f:
    nb1 = nbf.read(f, as_version=4)

for cell in nb1.cells:
    if cell.cell_type == "code" and "def get_cache_filename" in cell.source:
        cell.source = re.sub(
            r'return f"\{disease_slug\}_\{tissue_slug\}_\{cells_str\}_\{download_mode\}_\{census_version\}\.h5ad"',
            r'return f"{tissue_slug}_{cells_str}_{download_mode}_{census_version}.h5ad"',
            cell.source
        )

with open(cellxgene_nb, "w") as f:
    nbf.write(nb1, f)
print("Fixed prefix in cancer_cellxgene_integration.ipynb.")

# 3. Update primary_vs_metastasis_comparison.ipynb
pvsm_nb = "scripts/primary_vs_metastasis_comparison.ipynb"
with open(pvsm_nb, "r") as f:
    nb2 = nbf.read(f, as_version=4)

cells_to_keep = []
for cell in nb2.cells:
    # Remove internal nbconvert cell
    if cell.cell_type == "code" and "jupyter nbconvert" in cell.source and "output_base" in cell.source:
        print("Removed internal nbconvert cell from primary_vs_metastasis_comparison.ipynb")
        continue
    
    # Fix DE block to check for both groups
    if cell.cell_type == "code" and "sc.tl.rank_genes_groups(" in cell.source and "groupby='site'" in cell.source:
        if "if len(adata.obs['site'].unique()) < 2:" not in cell.source:
            replacement = """print("Running Differential Expression...")
if len(adata.obs['site'].unique()) < 2:
    print(f"Warning: Cannot perform Main DE. Dataset only contains {list(adata.obs['site'].unique())} cells for these PRIMARY_TISSUES.")
    import numpy as np
    adata.uns['rank_genes_groups'] = {
        'names': np.rec.fromarrays([np.array([], dtype=object)], names=['Metastasis']),
        'scores': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
        'logfoldchanges': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
        'pvals': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
        'pvals_adj': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis'])
    }
else:
    try:
        sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')
    except Exception as e:
        print(f"Warning: Cannot perform Main DE. Error: {e}")
        import numpy as np
        adata.uns['rank_genes_groups'] = {
            'names': np.rec.fromarrays([np.array([], dtype=object)], names=['Metastasis']),
            'scores': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
            'logfoldchanges': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
            'pvals': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis']),
            'pvals_adj': np.rec.fromarrays([np.array([], dtype=float)], names=['Metastasis'])
        }"""
            
            # Use regex to replace the entire try-except or the original call
            # This handles both the original and our previous try-except wrapper
            cell.source = re.sub(r'print\("Running Differential Expression..."\)\n(?:try:\n\s+sc\.tl\.rank_genes_groups.*?(?=\n# Extract DE results))', 
                   replacement + "\n", cell.source, flags=re.DOTALL)
            
            if "if len(adata.obs['site'].unique()) < 2:" not in cell.source:
                print("Regex replace failed for Main DE, using fallback.")
                # fallback if regex didn't match perfectly
                cell.source = cell.source.replace("sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')", 
                                                  "if len(adata.obs['site'].unique()) < 2:\n        pass # Handled by try-except logic replaced previously\n    sc.tl.rank_genes_groups(adata, groupby='site', groups=['Metastasis'], reference='Primary', method='wilcoxon')")
    
    cells_to_keep.append(cell)

nb2.cells = cells_to_keep
with open(pvsm_nb, "w") as f:
    nbf.write(nb2, f)
print("Fixed primary_vs_metastasis_comparison.ipynb DE logic and removed duplicate HTML export.")

# 4. Clean up output directories
output_dir = "output"
breast_results_dir = "output/breast_results"
os.makedirs(breast_results_dir, exist_ok=True)

# Delete stray HTML/CSV files in output/
for f in os.listdir(output_dir):
    p = os.path.join(output_dir, f)
    if os.path.isfile(p):
        if f.endswith('.html') or f.endswith('.csv') or f.endswith('.png') or f.endswith('.txt'):
            # Only delete if it looks like one of our generated stray files
            if "primary_vs_metastasis" in f or "orphan" in f or "immune_evasion" in f or "cellxgene" in f or "cancer_" in f:
                os.remove(p)
                print(f"Deleted stray file: {f}")

# Rename existing files in breast_cancer/ to move to breast_results/ and remove prefix
breast_cancer_dir = "output/breast_cancer"
if os.path.exists(breast_cancer_dir):
    for f in os.listdir(breast_cancer_dir):
        old_p = os.path.join(breast_cancer_dir, f)
        if os.path.isfile(old_p):
            new_f = f
            if f.startswith("cancer_breast-cancer_"):
                new_f = f.replace("cancer_breast-cancer_", "")
            new_p = os.path.join(breast_results_dir, new_f)
            shutil.move(old_p, new_p)
            print(f"Moved and renamed {f} to breast_results/{new_f}")
    
    # Try to remove empty directory
    try:
        os.rmdir(breast_cancer_dir)
        print(f"Removed empty directory {breast_cancer_dir}")
    except OSError:
        pass

print("Cleanup complete.")
