import json

path = "scripts/cancer_cellxgene_integration.ipynb"
with open(path, "r") as f:
    nb = json.load(f)

for cell in nb.get("cells", []):
    if cell["cell_type"] == "code":
        source = cell.get("source", [])
        for i, line in enumerate(source):
            if line.startswith("DISEASE_FILTER = ["):
                # Replace the block
                end_idx = i
                while end_idx < len(source) and not source[end_idx].startswith("]"):
                    end_idx += 1
                
                new_lines = [
                    'DISEASE_FILTER = [                     # Cancer types to query\n',
                    '    "lung adenocarcinoma",\n',
                    '    "breast cancer",\n',
                    '    "colorectal cancer",\n',
                    '    "ovarian cancer",\n',
                    '    "brain cancer",\n',
                    '    "melanoma"\n',
                    ']\n'
                ]
                source = source[:i] + new_lines + source[end_idx+1:]
                cell["source"] = source
                break

for cell in nb.get("cells", []):
    if cell["cell_type"] == "code":
        source = cell.get("source", [])
        for i, line in enumerate(source):
            if line.startswith('TISSUE_FILTER  = ["large intestine",'):
                end_idx = i
                while end_idx < len(source) and not source[end_idx].startswith('                 "intestine"]'):
                    end_idx += 1
                source = source[:i] + ['TISSUE_FILTER  = [] # None = all tissues\n'] + source[end_idx+1:]
                cell["source"] = source
                break

for cell in nb.get("cells", []):
    if cell["cell_type"] == "code":
        source = cell.get("source", [])
        for i, line in enumerate(source):
            if line.startswith("CAP = 100000"):
                source[i] = 'CAP = None #os.environ.get("CELLXGENE_CAP", None)\n'
                cell["source"] = source
                break

with open(path, "w") as f:
    json.dump(nb, f, indent=1)

print("Notebook updated!")
