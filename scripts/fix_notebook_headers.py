import json
import os

def fix_notebook(path):
    print(f"Fixing notebook: {path}")
    with open(path, "r", encoding="utf-8") as f:
        nb = json.load(f)
    
    modified = False
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            source_lines = cell.get("source", [])
            new_lines = []
            for line in source_lines:
                if "df['SuperClass']" in line:
                    line = line.replace("df['SuperClass']", "df['Super_Class']")
                    modified = True
                new_lines.append(line)
            cell["source"] = new_lines
            
    if modified:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
        print(f"🎉 Successfully updated 'SuperClass' to 'Super_Class' in {os.path.basename(path)}")
    else:
        print(f"No changes needed in {os.path.basename(path)}")

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    nb_path = os.path.join(script_dir, "metab_targetPair_analysis.ipynb")
    fix_notebook(nb_path)
