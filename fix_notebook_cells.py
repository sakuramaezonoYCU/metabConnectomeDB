import json
import glob
import os

for path in glob.glob("scripts/*.ipynb"):
    with open(path, "r") as f:
        nb = json.load(f)
    
    new_cells = []
    for cell in nb.get("cells", []):
        if cell["cell_type"] == "markdown":
            source = "".join(cell.get("source", []))
            if "### ⚙️ Pipeline Injected Configuration" in source:
                continue
        new_cells.append(cell)
        
    nb["cells"] = new_cells
    with open(path, "w") as f:
        json.dump(nb, f, indent=1)
    print(f"Cleaned {path}")
