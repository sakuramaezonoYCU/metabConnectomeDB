import json

with open('scripts/metab_targetPair_analysis.ipynb', 'r') as f:
    nb = json.load(f)

for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'markdown':
        src = ''.join(cell['source'])
        if '5.2' in src or '5.3' in src:
            print(f"--- MARKDOWN CELL {i} ---")
            print(src[:100])
    elif cell['cell_type'] == 'code':
        src = ''.join(cell['source'])
        if 'enzyme' in src.lower() or 'inter_clean' in src:
            print(f"--- CODE CELL {i} ---")
            print(src)
