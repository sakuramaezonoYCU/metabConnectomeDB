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

# Dynamic path resolution relative to this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
db_dir = os.path.join(PROJECT_ROOT, 'input', 'databases')
out_dir = os.path.join(PROJECT_ROOT, 'output')

# ------------------------------------------------------------------
# 1. Load reference files
# ------------------------------------------------------------------
print("Loading HMDB reference files...")

# Full HMDB annotations (primary source — 248k entries)
hmdb_path = os.path.join(db_dir, 'HMDB_metabolites')
df_hmdb = pd.read_csv(hmdb_path, sep=',')
df_hmdb = df_hmdb.rename(columns={
    'NAME': 'HMDB_Name',
    'SuperClass': 'Super_Class',
    'SubClass': 'Sub_Class'
})
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
    print(f"\n--- {species.capitalize()} ({'Target Pairs' if 'target_pairs' in in_file else 'Unique Metabolites'}) ---")
    print(f"  Input rows: {len(df):,} | Missing HMDB_ID: {df['HMDB_ID'].isna().sum()}")

    # Normalize Metabolite_Name
    if 'Metabolite_Name' in df.columns:
        df['Metabolite_Name'] = df['Metabolite_Name'].astype(str).str.strip().str.lower()

    if "target_pairs" in in_file:
        # Load the completed unique metabolite reference file (which was processed first)
        ref_file = in_file.replace("_unique_metab_target_pairs.csv", "_unique_metab_with_HMDB_Info.csv")
        if os.path.exists(ref_file):
            print(f"  Mapping HMDB_IDs directly from clean reference file: {os.path.basename(ref_file)}...")
            df_ref = pd.read_csv(ref_file)
            df_ref['clean_name'] = df_ref['Metabolite_Name'].astype(str).str.lower().str.strip()
            
            # Map each clean name to its list of HMDB_IDs from the reference
            ref_map = df_ref.dropna(subset=['HMDB_ID']).groupby('clean_name')['HMDB_ID'].apply(lambda x: sorted(list(set(x)))).to_dict()
            
            # Clean up target pair names and map
            df['clean_name'] = df['Metabolite_Name'].astype(str).str.lower().str.strip()
            df['HMDB_ID'] = df['clean_name'].map(ref_map)
            
            # For any missing/unmapped, fall back to lookup_map
            def fallback_id(row):
                val = row['HMDB_ID']
                if isinstance(val, list) and len(val) > 0:
                    return val
                lk = lookup_map.get(row['clean_name'])
                if lk:
                    # Clean/split pipe or commas
                    cleaned_lk = [h.strip() for h in str(lk).replace('|', ',').split(',') if h.strip() and h.strip() != 'nan']
                    return cleaned_lk if cleaned_lk else [None]
                return [None]
                
            df['HMDB_ID'] = df.apply(fallback_id, axis=1)
            df = df.explode('HMDB_ID').reset_index(drop=True)
            df = df.drop(columns=['clean_name'])
            print(f"  After mapping and exploding target pair HMDB_IDs: {len(df):,} rows")
        else:
            print("  [WARNING] Reference file not found. Running fallback name-match explode...")
            # Original name-match fill
            def fill_hmdb_id(row):
                if pd.notna(row['HMDB_ID']):
                    # clean any pipe characters
                    return str(row['HMDB_ID']).replace('|', ',')
                key = str(row['Metabolite_Name']).lower().strip()
                return lookup_map.get(key, None)

            df['HMDB_ID'] = df.apply(fill_hmdb_id, axis=1)
            df['HMDB_ID'] = df['HMDB_ID'].apply(
                lambda x: [h.strip() for h in str(x).replace('|', ',').split(',') if h.strip() != 'nan'] if pd.notna(x) else [None]
            )
            df = df.explode('HMDB_ID').reset_index(drop=True)
    else:
        # 1. Fill missing HMDB_ID using Metabolite_Name from name_to_hmdb within the file itself
        mapping_df = df[['Metabolite_Name', 'HMDB_ID']].dropna()
        mapping_df = mapping_df.drop_duplicates(subset='Metabolite_Name')
        name_to_hmdb = dict(zip(mapping_df['Metabolite_Name'].astype(str).str.lower().str.strip(),
                                mapping_df['HMDB_ID']))

        mask = df['HMDB_ID'].isna() & df['Metabolite_Name'].notna()
        df.loc[mask, 'HMDB_ID'] = df.loc[mask, 'Metabolite_Name'].map(name_to_hmdb)

        # 2. Fill remaining via name lookup (lookup_map)
        def fill_hmdb_id(row):
            if pd.notna(row['HMDB_ID']):
                return str(row['HMDB_ID']).replace('|', ',')
            key = str(row['Metabolite_Name']).lower().strip()
            return lookup_map.get(key, None)

        df['HMDB_ID'] = df.apply(fill_hmdb_id, axis=1)
        still_missing = df['HMDB_ID'].isna().sum()
        filled = df['HMDB_ID'].notna().sum()
        print(f"  After name-match fill → {filled:,} with HMDB_ID, {still_missing:,} still missing")

        # 3. Explode rows with multiple comma-separated HMDB_IDs
        df['HMDB_ID'] = df['HMDB_ID'].apply(
            lambda x: [h.strip() for h in str(x).replace('|', ',').split(',') if h.strip() != 'nan'] if pd.notna(x) else [None]
        )
        df = df.explode('HMDB_ID').reset_index(drop=True)

    # Normalize empty strings back to NaN
    df['HMDB_ID'] = df['HMDB_ID'].replace('', None)
    print(f"  After HMDB_ID clean/explode: {len(df):,} rows")

    # ------------------------------------------------------------------
    # Step 3 — Left-join full HMDB annotations by HMDB_ID
    # ------------------------------------------------------------------
    # Identify overlapping columns (excluding the join key) to prevent duplicates/suffixes
    overlapping_cols = [c for c in df_hmdb.columns if c in df.columns and c != 'HMDB_ID']
    
    df = df.merge(df_hmdb, on='HMDB_ID', how='left')
    
    # Surgical coalescing of any duplicate column suffixes (e.g. SMILES, Super_Class, Class)
    for c in overlapping_cols:
        x_col = f"{c}_x"
        y_col = f"{c}_y"
        if x_col in df.columns and y_col in df.columns:
            df[c] = df[x_col].fillna(df[y_col])
            df = df.drop(columns=[x_col, y_col])
            
    print(f"  After HMDB annotation join and coalescing: {len(df):,} rows")
    matched = df['HMDB_Name'].notna().sum()
    print(f"  Rows with matched HMDB_Name: {matched:,}")

    # Map Super_Class and Sub_Class from target pairs for the unique metabolite file
    if "target_pairs" not in in_file:
        pair_file = in_file.replace("_unique_metab.csv", "_unique_metab_target_pairs.csv")
        if os.path.exists(pair_file):
            print(f"  Mapping Super_Class and Sub_Class from {os.path.basename(pair_file)}...")
            df_pair = pd.read_csv(pair_file)
            if 'Super_Class' in df_pair.columns:
                sc_map = df_pair.dropna(subset=['Super_Class']).drop_duplicates(subset='Metabolite_Name').set_index('Metabolite_Name')['Super_Class'].to_dict()
                df['Super_Class'] = df['Metabolite_Name'].map(sc_map)
            if 'Sub_Class' in df_pair.columns:
                sub_map = df_pair.dropna(subset=['Sub_Class']).drop_duplicates(subset='Metabolite_Name').set_index('Metabolite_Name')['Sub_Class'].to_dict()
                df['Sub_Class'] = df['Metabolite_Name'].map(sub_map)

    # Ensure Super_Class is standardized to ClassyFire standard
    if 'Super_Class' in df.columns:
        from standardize_categories import standardize_superclass
        df['Super_Class'] = df['Super_Class'].apply(standardize_superclass)

    # Ensure Sub_Class is standardized to ClassyFire standard
    if 'Sub_Class' in df.columns:
        from standardize_categories import standardize_subclass
        df['Sub_Class'] = df['Sub_Class'].apply(standardize_subclass)

    # ------------------------------------------------------------------
    # Step 4 — Save
    # ------------------------------------------------------------------
    df.to_csv(out_file, index=False)
    print(f"  -> Saved {out_file}")


# ------------------------------------------------------------------
# Run for both species
# ------------------------------------------------------------------

if __name__ == "__main__":
    for sp in ['human', 'mouse']:
        # 1. Unique metabolites
        in_file = os.path.join(out_dir, f'{sp}_database_merge_unique_metab.csv')
        process_species(sp, in_file)
        
        # 2. Target pairs
        in_file = os.path.join(out_dir, f'{sp}_database_merge_unique_metab_target_pairs.csv')
        process_species(sp, in_file)

    print("\nDone.")
