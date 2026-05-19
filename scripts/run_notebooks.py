import json
import os
import sys

# Configure Matplotlib for headless non-blocking execution before any other import
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.show = lambda *args, **kwargs: None

def run_notebook(path):
    print(f"Executing notebook: {os.path.basename(path)}...")
    
    # Prepend venv/bin to PATH so subprocess calls inside the notebook use the venv binaries
    venv_bin = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(path))), 'venv', 'bin')
    if os.path.exists(venv_bin):
        os.environ['PATH'] = venv_bin + os.pathsep + os.environ.get('PATH', '')
        
    with open(path, "r", encoding="utf-8") as f:
        nb = json.load(f)
    
    # We maintain a shared global dict for the execution of cells
    global_dict = {
        '__name__': '__main__',
        'plt': plt,  # Force overridden plt into cell namespace
    }
    
    # Change directory to the directory of the notebook
    notebook_dir = os.path.dirname(os.path.abspath(path))
    original_dir = os.getcwd()
    os.chdir(notebook_dir)
    
    try:
        for idx, cell in enumerate(nb.get("cells", [])):
            if cell.get("cell_type") == "code":
                source_lines = cell.get("source", [])
                # Clean source lines from magic commands
                cleaned_lines = []
                for line in source_lines:
                    # Ignore lines starting with % or !
                    stripped = line.strip()
                    if stripped.startswith("%") or stripped.startswith("!"):
                        continue
                    cleaned_lines.append(line)
                
                cell_code = "".join(cleaned_lines)
                if not cell_code.strip():
                    continue
                
                try:
                    # Execute the cell
                    exec(cell_code, global_dict)
                except Exception as e:
                    print(f"❌ Error executing cell {idx} in {os.path.basename(path)}:")
                    print("--- CELL CODE START ---")
                    print(cell_code)
                    print("--- CELL CODE END ---")
                    print(f"Error: {e}")
                    raise e
        print(f"🎉 Successfully executed notebook: {os.path.basename(path)}\n")
    finally:
        os.chdir(original_dir)

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    notebooks = [
        os.path.join(script_dir, "unique_metab_data_exploration.ipynb"),
        os.path.join(script_dir, "metab_targetPair_analysis.ipynb")
    ]
    for nb in notebooks:
        run_notebook(nb)
