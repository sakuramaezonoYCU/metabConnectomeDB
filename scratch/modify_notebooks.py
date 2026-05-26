import json
import os

# 1. Modify primary_vs_metastasis_comparison.ipynb
notebook_path1 = 'scripts/primary_vs_metastasis_comparison.ipynb'
with open(notebook_path1, 'r', encoding='utf-8') as f:
    nb1 = json.load(f)

modified1 = False
for idx, cell in enumerate(nb1.get('cells', [])):
    if cell.get('cell_type') == 'code':
        source = "".join(cell.get('source', []))
        if "metab_db_path = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')" in source:
            print(f"Modifying notebook 1 cell {idx} for metab_db_path")
            source = source.replace(
                "metab_db_path = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')",
                "metab_db_path = os.path.join(BASE_DIR, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')"
            )
            cell['source'] = [source]
            modified1 = True

if modified1:
    with open(notebook_path1, 'w', encoding='utf-8') as f:
        json.dump(nb1, f, indent=1)
    print("🎉 Comparison notebook updated successfully!")
else:
    print("⚠️ No cells matched the modification target text in Comparison notebook.")


# 2. Modify orphan_metabolic_immune_evasion.ipynb
notebook_path2 = 'scripts/orphan_metabolic_immune_evasion.ipynb'
with open(notebook_path2, 'r', encoding='utf-8') as f:
    nb2 = json.load(f)

modified2 = False
for idx, cell in enumerate(nb2.get('cells', [])):
    if cell.get('cell_type') == 'code':
        source = "".join(cell.get('source', []))
        if "pairs_df = pd.read_csv(os.path.join(output_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'))" in source:
            print(f"Modifying notebook 2 cell {idx} for pairs_df")
            source = source.replace(
                "pairs_df = pd.read_csv(os.path.join(output_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'))",
                "pairs_df = pd.read_csv(os.path.join(workspace_dir, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'))"
            )
            cell['source'] = [source]
            modified2 = True

if modified2:
    with open(notebook_path2, 'w', encoding='utf-8') as f:
        json.dump(nb2, f, indent=1)
    print("🎉 Evasion notebook updated successfully!")
else:
    print("⚠️ No cells matched the modification target text in Evasion notebook.")
