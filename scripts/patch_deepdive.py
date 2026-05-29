import nbformat as nbf

nb_path = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/deepdive_conserved_metabGeneSig.ipynb"

try:
    with open(nb_path, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
        
    modified = False
    for cell in nb.cells:
        old_source = cell.source
        
        # Scrub 23s
        cell.source = cell.source.replace("23-Gene Metastatic", "Strictly Conserved Metastatic")
        cell.source = cell.source.replace("23-Gene Score", "Metastatic Signature Score")
        cell.source = cell.source.replace("all 23 pan-cancer", "all strictly conserved pan-cancer")
        cell.source = cell.source.replace("stat3_u87_targets_23genes.csv", "stat3_u87_targets_strictly_conserved.csv")
        cell.source = cell.source.replace("23-Gene Metastatic Score", "Metastatic Signature Score")
        cell.source = cell.source.replace("Random 23-Gene Signatures", "Random Gene Signatures")
        
        if cell.source != old_source:
            modified = True
            
        # Clear outputs
        if cell.cell_type == 'code':
            cell.outputs = []
            
    if modified:
        with open(nb_path, 'w', encoding='utf-8') as f:
            nbf.write(nb, f)
        print("Patched deepdive_conserved_metabGeneSig.ipynb")
        
except Exception as e:
    print(e)
