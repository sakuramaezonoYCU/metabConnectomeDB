import json

notebook_path = "scripts/cancer_cellxgene_integration.ipynb"

with open(notebook_path, 'r', encoding='utf-8') as f:
    nb = json.load(f)

for cell in nb.get('cells', []):
    if cell['cell_type'] == 'markdown':
        source = cell.get('source', [])
        # Find the cell with Section 9
        if any("## 9. Cancer Pathway-Level Analysis" in line for line in source):
            # Replace the table rows
            new_source = []
            for line in source:
                if "| Pathway | Key Metabolite | Key Genes | Cancer Relevance |" in line:
                    new_source.append("| Pathway | Key Metabolite | Key Genes | Cancer Relevance | References (PMID) |\n")
                elif "|---------|---------------|-----------|-----------------|" in line:
                    new_source.append("|---------|---------------|-----------|-----------------|-------------------|\n")
                elif "| IDO1/Kynurenine |" in line:
                    new_source.append("| IDO1/Kynurenine | Kynurenine, Tryptophan | IDO1, TDO2, AHR, KMO | Immune suppression via T cell exhaustion | [29551271](https://pubmed.ncbi.nlm.nih.gov/29551271/), [34322129](https://pubmed.ncbi.nlm.nih.gov/34322129/) |\n")
                elif "| xCT/Glutamate |" in line:
                    new_source.append("| xCT/Glutamate | Glutamate, Cystine | SLC7A11, SLC3A2 | Ferroptosis resistance, oxidative stress | [36496662](https://pubmed.ncbi.nlm.nih.gov/36496662/), [30872535](https://pubmed.ncbi.nlm.nih.gov/30872535/) |\n")
                elif "| CD73/Adenosine |" in line:
                    new_source.append("| CD73/Adenosine | Adenosine, ATP | NT5E (CD73), ADORA2A, ENTPD1 (CD39) | Purinergic immune checkpoint | [26880461](https://pubmed.ncbi.nlm.nih.gov/26880461/), [29758241](https://pubmed.ncbi.nlm.nih.gov/29758241/) |\n")
                elif "| COX-2/PGE" in line:
                    new_source.append("| COX-2/PGE₂ | PGE₂, Arachidonic acid | PTGS2 (COX-2), PTGER2, PTGER4 | Inflammation, immune evasion | [35121582](https://pubmed.ncbi.nlm.nih.gov/35121582/), [36776289](https://pubmed.ncbi.nlm.nih.gov/36776289/) |\n")
                elif "| SPHK1/S1P |" in line:
                    new_source.append("| SPHK1/S1P | S1P | SPHK1, S1PR1-5 | Tumor angiogenesis, lymphocyte trafficking | [22298596](https://pubmed.ncbi.nlm.nih.gov/22298596/) |\n")
                else:
                    new_source.append(line)
            cell['source'] = new_source

with open(notebook_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1)

print("Notebook successfully updated with PMIDs in Section 9!")
