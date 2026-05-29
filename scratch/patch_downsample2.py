import json
import re

nb_path = 'scripts/cancer_cellxgene_integration.ipynb'
with open(nb_path, 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "obs_df_sorted.head(CAP)" in source:
            # We want to replace the whole block starting from "obs_df_sorted = obs_df.sort_values" 
            # to "obs_df_sub = obs_df_sorted.copy()"
            new_source = re.sub(
                r"obs_df_sorted = obs_df\.sort_values\('soma_joinid'\).*?obs_df_sub = obs_df_sorted\.copy\(\)",
                "if CAP is not None and len(obs_df) > CAP:\n        print(f\"⚠️  Downsampling metadata from {len(obs_df):,} to {CAP:,} cells (stratified random sample)...\")\n        obs_df_sub = obs_df.groupby('tissue_general', group_keys=False).apply(lambda x: x.sample(n=min(len(x), CAP // len(obs_df['tissue_general'].unique())), random_state=42))\n        obs_df_sub = obs_df_sub.sort_values('soma_joinid')\n    else:\n        print(f\"✅ Downloading all available cells ({len(obs_df):,})...\")\n        obs_df_sub = obs_df.sort_values('soma_joinid')",
                source,
                flags=re.DOTALL
            )
            if new_source != source:
                lines = new_source.split('\n')
                cell['source'] = [line + '\n' for line in lines[:-1]] + [lines[-1]] if lines[-1] else [line + '\n' for line in lines[:-1]]
                print("Patched downsampling logic successfully!")
                break

with open(nb_path, 'w') as f:
    json.dump(nb, f, indent=1)
