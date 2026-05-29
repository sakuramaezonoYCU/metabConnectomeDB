import json

nb_path = 'scripts/cancer_cellxgene_integration.ipynb'
with open(nb_path, 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "obs_df_sorted.head(CAP)" in source:
            new_source = source.replace(
                "# 1. Subsample cell IDs FIRST to select a CONTIGUOUS range of cells (prevents remote seek lag)\n    obs_df_sorted = obs_df.sort_values('soma_joinid')\n    if CAP is not None and len(obs_df_sorted) > CAP:\n        print(f\"⚠️  Downsampling metadata from {len(obs_df_sorted):,} to {CAP:,} cells (contiguous slice)...\")\n        obs_df_sub = obs_df_sorted.head(CAP)\n    else:\n        print(f\"✅ Downloading all available cells ({len(obs_df_sorted):,})...\")\n        obs_df_sub = obs_df_sorted.copy()",
                "# 1. Subsample cell IDs FIRST (use random sample to keep tissues proportional)\n    if CAP is not None and len(obs_df) > CAP:\n        print(f\"⚠️  Downsampling metadata from {len(obs_df):,} to {CAP:,} cells (stratified random sample)...\")\n        obs_df_sub = obs_df.sample(n=CAP, random_state=42).sort_values('soma_joinid')\n    else:\n        print(f\"✅ Downloading all available cells ({len(obs_df):,})...\")\n        obs_df_sub = obs_df.sort_values('soma_joinid')"
            )
            if new_source != source:
                lines = new_source.split('\n')
                cell['source'] = [line + '\n' for line in lines[:-1]] + [lines[-1]] if lines[-1] else [line + '\n' for line in lines[:-1]]
                print("Patched downsampling logic successfully!")
                break

with open(nb_path, 'w') as f:
    json.dump(nb, f, indent=1)
