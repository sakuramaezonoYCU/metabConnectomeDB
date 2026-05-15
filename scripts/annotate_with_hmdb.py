"""
annotate_with_hmdb.py
---------------------
Enriches {species}_database_merge_unique_metab.csv with full HMDB annotation.

Steps:
  1. Load HMDB_ID_Metab.csv (name → HMDB_ID lookup) and HMDB_metabolites
     (HMDB_ID → full annotations: NAME, SMILES, INCHIKEY, etc.)
  2. For rows where HMDB_ID is missing, attempt a case-insensitive name match
     of the 'metabolite' column against HMDB_ID_Metab 'Metabolite_Name'.
  3. Expand rows whose HMDB_ID contains multiple comma-separated IDs into
     one row per ID.
  4. Left-join HMDB_metabolites onto the expanded table by HMDB_ID, renaming
     'NAME' → 'HMDB_Name'.
  5. Save as {species}_database_merge_unique_metab_with_HMDB_Info.csv.

Note: This script is standalone and does not modify generate_final_outputs.py,
      but may be merged with it later.
"""

import os
import pandas as pd

# BASE_DIR = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/input/databases'
# os.chdir(BASE_DIR)

# ------------------------------------------------------------------
# 1. Load reference files
# ------------------------------------------------------------------
print("Loading HMDB reference files...")

# Full HMDB annotations (primary source — 248k entries)
df_hmdb = pd.read_csv('HMDB_metabolites', sep=',')
df_hmdb = df_hmdb.rename(columns={'NAME': 'HMDB_Name'})
print(f"  HMDB annotations: {len(df_hmdb):,} entries | cols: {df_hmdb.columns.tolist()}")

# Build a normalised (lowercase, stripped) name → HMDB_ID map from HMDB_metabolites
# This is the correct source per step #1: match 'metabolite' against 'NAME'
lookup_map = dict(
    zip(
        df_hmdb['HMDB_Name'].str.lower().str.strip(),
        df_hmdb['HMDB_ID']
    )
)
print(f"  Name lookup map built: {len(lookup_map):,} entries")


# ------------------------------------------------------------------
# Processing function (one species at a time)
# ------------------------------------------------------------------
def process_species(species: str, in_file):
    # derive from in_file
    out_file = in_file.replace('.csv', '_with_HMDB_Info.csv')

    if not os.path.exists(in_file):
        print(f"  [SKIP] {in_file} not found.")
        return

    df = pd.read_csv(in_file)
        # ------------------------------------------------------------------
    # Step 1 — Fill missing HMDB_ID via name lookup (for within file)
    # ------------------------------------------------------------------
    # --- build mapping from valid rows ---
    mapping_df = df[['Metabolite_Name', 'HMDB_ID']].dropna()

    # ensure unique mapping (drop duplicates if any)
    mapping_df = mapping_df.drop_duplicates(subset='Metabolite_Name')

    name_to_hmdb = dict(zip(mapping_df['Metabolite_Name'],
                            mapping_df['HMDB_ID']))

    # --- fill missing HMDB_ID using Metabolite_Name ---
    mask = df['HMDB_ID'].isna() & df['Metabolite_Name'].notna()

    df.loc[mask, 'HMDB_ID'] = df.loc[mask, 'Metabolite_Name'].map(name_to_hmdb)
   # ------------------------------------------------------------------

    
    print(f"\n--- {species.capitalize()} ---")
    print(f"  Input rows: {len(df):,} | Missing HMDB_ID: {df['HMDB_ID'].isna().sum()}")

    # ------------------------------------------------------------------
    # Step 1 — Fill missing HMDB_ID via name lookup
    # ------------------------------------------------------------------
    def fill_hmdb_id(row):
        if pd.notna(row['HMDB_ID']):
            return row['HMDB_ID']
        key = str(row['Metabolite_Name']).lower().strip()
        return lookup_map.get(key, None)

    df['HMDB_ID'] = df.apply(fill_hmdb_id, axis=1)
    still_missing = df['HMDB_ID'].isna().sum()
    filled = df['HMDB_ID'].notna().sum()
    print(f"  After name-match fill → {filled:,} with HMDB_ID, {still_missing:,} still missing")

    # ------------------------------------------------------------------
    # Step 2 — Explode rows with multiple comma-separated HMDB_IDs
    # ------------------------------------------------------------------
    df['HMDB_ID'] = df['HMDB_ID'].apply(
        lambda x: [h.strip() for h in str(x).split(',') if h.strip() != 'nan'] if pd.notna(x) else [None]
    )
    df = df.explode('HMDB_ID').reset_index(drop=True)
    # Normalise empty strings back to NaN
    df['HMDB_ID'] = df['HMDB_ID'].replace('', None)
    print(f"  After HMDB_ID explode: {len(df):,} rows")

    # ------------------------------------------------------------------
    # Step 3 — Left-join full HMDB annotations by HMDB_ID
    # ------------------------------------------------------------------
    df = df.merge(df_hmdb, on='HMDB_ID', how='left')
    print(f"  After HMDB annotation join: {len(df):,} rows")
    matched = df['HMDB_Name'].notna().sum()
    print(f"  Rows with matched HMDB_Name: {matched:,}")

    # ------------------------------------------------------------------
    # Step 4 — Save
    # ------------------------------------------------------------------
    df.to_csv(out_file, index=False)
    print(f"  -> Saved {out_file}")


# ------------------------------------------------------------------
# Run for both species
# ------------------------------------------------------------------

for sp in ['human', 'mouse']:
    in_file = f'{sp}_database_merge_unique_metab.csv'
    process_species(sp, in_file)
    in_file = f'{sp}_database_merge_unique_metab_target_pairs.csv'
    process_species(sp, in_file)

print("\nDone.")
