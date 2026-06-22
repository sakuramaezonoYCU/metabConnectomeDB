import json

with open('scripts/oxygen_tension_analysis.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb.get('cells', []):
    if cell.get('cell_type') == 'code':
        source = cell.get('source', [])
        for i, line in enumerate(source):
            if "# Hypotheses 1: Normal Tissue Physioxia" in line:
                source.insert(i, "df_plot_tumour['Cancer'] = df_plot_tumour['Cancer'].astype(str)\n")
                source.insert(i, "df_plot_normal['Cancer'] = df_plot_normal['Cancer'].astype(str)\n")
                break

# Insert KEGG provenance markdown cell at index 1 (right after the first markdown cell)
provenance_cell = {
 "cell_type": "markdown",
 "id": "kegg_provenance_block",
 "metadata": {},
 "source": [
  "### Data Provenance: KEGG Pathway Definition\n",
  "\n",
  "To ensure strict scientific reproducibility, the gene signatures for **Glycolysis** (hsa00010), **Oxidative Phosphorylation** (hsa00190), and **HIF-1 Signaling** (hsa04066) are programmatically fetched from the KEGG REST API rather than being hardcoded.\n",
  "\n",
  "The extraction script (`scripts/fetch_kegg_pathways.py`) retrieves the HGNC symbols corresponding to each KEGG pathway to precisely define our analytical signatures."
 ]
}

# Only insert if not already present
if not any(c.get('id') == 'kegg_provenance_block' for c in nb.get('cells', [])):
    nb['cells'].insert(1, provenance_cell)

with open('scripts/oxygen_tension_analysis.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

print("Notebook patched!")
