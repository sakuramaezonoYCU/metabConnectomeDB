import json
import re

nb_path = "scripts/cancer_cellxgene_integration.ipynb"
with open(nb_path, "r") as f:
    nb = json.load(f)

for cell in nb.get("cells", []):
    if cell.get("cell_type") != "code":
        continue
    
    source = cell["source"]
    source_str = "".join(source)
    
    # 1. Inject LOADED log
    if "LOAD_FROM_CACHE = os.path.exists(h5ad_path)" in source_str:
        if "h5ad_version_log.txt" not in source_str:
            new_lines = []
            for line in source:
                new_lines.append(line)
                if "print(f\"✅ Local H5AD cache found at: {h5ad_path}\")" in line:
                    indent = re.match(r"^(\s*)", line).group(1)
                    new_lines.append(f"{indent}import datetime\n")
                    new_lines.append(f"{indent}with open(os.path.join(output_dir, 'h5ad_version_log.txt'), 'a') as f:\n")
                    new_lines.append(f"{indent}    f.write(f\"[{{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}}] LOADED CACHE: {{h5ad_filename}} (Census Version: {{CENSUS_VERSION}})\\n\")\n")
            cell["source"] = new_lines

    # 2. Inject CREATED log
    if "adata.write_h5ad(h5ad_path, compression='gzip')" in source_str:
        if "h5ad_version_log.txt" not in source_str:
            new_lines = []
            for line in source:
                new_lines.append(line)
                if "print(\"✅ AnnData successfully saved to disk!\")" in line:
                    indent = re.match(r"^(\s*)", line).group(1)
                    new_lines.append(f"{indent}import datetime\n")
                    new_lines.append(f"{indent}with open(os.path.join(output_dir, 'h5ad_version_log.txt'), 'a') as f:\n")
                    new_lines.append(f"{indent}    f.write(f\"[{{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}}] CREATED NEW: {{h5ad_filename}} (Census Version: {{CENSUS_VERSION}}, Cells: {{adata.n_obs}})\\n\")\n")
            cell["source"] = new_lines

with open(nb_path, "w") as f:
    json.dump(nb, f, indent=1)
    
print("Successfully patched notebook to add h5ad tracking log.")
