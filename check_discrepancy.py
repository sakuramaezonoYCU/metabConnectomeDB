import pandas as pd
import os

base_dir = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB"
unique_path = os.path.join(base_dir, "output/human_database_merge_unique_metab_with_HMDB_Info.csv")
target_path = os.path.join(base_dir, "output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv")

unique_metab_df = pd.read_csv(unique_path)
target_pairs_df = pd.read_csv(target_path)

print("--- Unique Metabolites File ---")
print(f"Total rows: {len(unique_metab_df)}")
print(f"Unique Metabolite_Name: {unique_metab_df['Metabolite_Name'].nunique()}")
if 'HMDB_ID' in unique_metab_df.columns:
    print(f"Unique HMDB_ID: {unique_metab_df['HMDB_ID'].nunique()}")
print("Columns:", unique_metab_df.columns.tolist())

print("\n--- Target Pairs File ---")
print(f"Total rows: {len(target_pairs_df)}")
print(f"Unique Metabolite_Name: {target_pairs_df['Metabolite_Name'].nunique()}")
if 'HMDB_ID' in target_pairs_df.columns:
    print(f"Unique HMDB_ID: {target_pairs_df['HMDB_ID'].nunique()}")
print("Columns:", target_pairs_df.columns.tolist())

unique_metabs_in_list = set(unique_metab_df['HMDB_ID'].dropna().unique())
unique_metabs_in_pairs = set(target_pairs_df['HMDB_ID'].dropna().unique())

diff_list_only = unique_metabs_in_list - unique_metabs_in_pairs
diff_pairs_only = unique_metabs_in_pairs - unique_metabs_in_list

print(f"\nHMDB IDs in unique list but NOT in target pairs: {len(diff_list_only)}")
if len(diff_list_only) > 0:
    print("Sample unique-list only:", sorted(list(diff_list_only))[:100])

print(f"HMDB IDs in target pairs but NOT in unique list: {len(diff_pairs_only)}")
if len(diff_pairs_only) > 0:
    print("Sample target-pairs only:", sorted(list(diff_pairs_only))[:100])
