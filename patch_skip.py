import re
import os

with open("scripts/run_cancer_pipeline.py", "r", encoding="utf-8") as f:
    code = f.read()

# Add skip_if_exists parameter to execute_and_export
old_def = "def execute_and_export(notebook_path, html_path, title_text, global_dict):"
new_def = "def execute_and_export(notebook_path, html_path, title_text, global_dict, skip_if_exists=False):"
code = code.replace(old_def, new_def)

# Add skip logic at the beginning of execute_and_export
skip_logic = """
    if skip_if_exists and os.path.exists(html_path):
        print(f"⏭️ Skipping {os.path.basename(notebook_path)} as {os.path.basename(html_path)} already exists.")
        
        # If this is the cellxgene notebook, we need to recover the h5ad_path
        if 'cancer_cellxgene_integration' in notebook_path:
            # First try parsing the HTML
            try:
                with open(html_path, 'r', encoding='utf-8') as hf:
                    html_content = hf.read()
                    import re
                    match = re.search(r'Saved combined dataset to: ([^<\\n]+)', html_content)
                    if match:
                        return match.group(1).strip()
            except Exception as e:
                pass
                
            # Fallback to finding the .h5ad file in the output directory
            import glob
            output_dir = os.path.dirname(html_path)
            h5ads = glob.glob(os.path.join(output_dir, '*.h5ad'))
            if h5ads:
                return h5ads[0]
                
        # For other notebooks, just return whatever h5ad_path was passed in
        return global_dict.get('h5ad_path', None)
"""

# Insert skip logic after docstring/initialization
code = code.replace(
    "    print(f\"\\n{'='*50}\\nExecuting: {title_text}\\n{'='*50}\")",
    skip_logic + "\n    print(f\"\\n{'='*50}\\nExecuting: {title_text}\\n{'='*50}\")"
)

# Modify calls to pass skip_if_exists=True
code = code.replace(
    "generated_h5ad = execute_and_export(cellxgene_nb, cellxgene_html, f\"CellxGene Integration: {disease_filter_str}\", inject_globals)",
    "generated_h5ad = execute_and_export(cellxgene_nb, cellxgene_html, f\"CellxGene Integration: {disease_filter_str}\", inject_globals, skip_if_exists=True)"
)

code = code.replace(
    "execute_and_export(pvsm_nb, pvsm_html, f\"Primary vs Metastasis: {disease_filter_str}\", pvsm_globals)",
    "execute_and_export(pvsm_nb, pvsm_html, f\"Primary vs Metastasis: {disease_filter_str}\", pvsm_globals, skip_if_exists=True)"
)

code = code.replace(
    "execute_and_export(orphan_nb, orphan_html, f\"Orphan Metabolic Targets: {disease_filter_str}\", orphan_globals)",
    "execute_and_export(orphan_nb, orphan_html, f\"Orphan Metabolic Targets: {disease_filter_str}\", orphan_globals, skip_if_exists=True)"
)

with open("scripts/run_cancer_pipeline.py", "w", encoding="utf-8") as f:
    f.write(code)
print("SUCCESS")
