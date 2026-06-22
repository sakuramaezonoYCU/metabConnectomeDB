import json
import ast

def extract_function_source(filename, function_name):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    source = []
    in_function = False
    for line in lines:
        if line.startswith(f"def {function_name}("):
            in_function = True
        if in_function:
            source.append(line)
            # A bit simplistic, but we know the next def or main block ends it
            if len(source) > 1 and (line.startswith("def ") or line.startswith("# ===================")):
                source.pop()
                break
    return "".join(source)

# 1. Extract the function source code
func_source = extract_function_source('scripts/annotate_with_databases.py', 'classify_by_rules')

# 2. Prepare the code cell content
cell_content = [
    "# Note: The `classify_by_rules` function is presented here as a transparent guide for users ",
    "to see the exact naming heuristics used to classify targets.\n",
    "# Per project integrity rules, this function is actually executed globally across the dataset ",
    "during the data enrichment pipeline (`scripts/annotate_with_databases.py`) before this notebook is run.\n\n"
]
cell_content.extend([line + "\n" if not line.endswith("\n") else line for line in func_source.splitlines()])

# 3. Load the notebook
notebook_path = 'scripts/metab_targetPair_analysis.ipynb'
with open(notebook_path, 'r') as f:
    nb = json.load(f)

# 4. Find where to insert it (after the markdown cell)
insert_idx = -1
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'markdown' and 'Local Rule Classifier' in "".join(cell['source']):
        insert_idx = i + 1
        break

if insert_idx == -1:
    print("Error: Could not find the markdown cell.")
    exit(1)

# Check if the next cell is already classify_by_rules
next_cell = nb['cells'][insert_idx] if insert_idx < len(nb['cells']) else None
if next_cell and next_cell['cell_type'] == 'code' and 'classify_by_rules' in "".join(next_cell['source']):
    print("Cell already exists!")
    exit(0)

# 5. Create and insert the cell
new_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": cell_content
}

nb['cells'].insert(insert_idx, new_cell)

# 6. Save the notebook
with open(notebook_path, 'w') as f:
    json.dump(nb, f, indent=1)

print("Successfully restored classify_by_rules code cell to the notebook.")
