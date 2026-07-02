import pandas as pd
import numpy as np
import os
import glob


def normalize_empty(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # Select only object/string columns (where messy values exist)
    obj_cols = df.select_dtypes(include=["object", "string"]).columns

    # Normalize strings: strip whitespace + uppercase
    df[obj_cols] = df[obj_cols].apply(lambda s: s.str.strip().str.upper())

    # Replace common empty markers with NaN
    df[obj_cols] = df[obj_cols].replace(
        {"": np.nan, "NA": np.nan, "N/A": np.nan, "NULL": np.nan}
    )

    return df

def compute_raw_counts():
    db_counts = {}
    base_dir = os.path.join('input', 'databases')
    if not os.path.exists(base_dir): return db_counts
    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        if os.path.isdir(folder_path) and folder not in ['.', '..', '.DS_Store']:
            folder_files = glob.glob(f'{folder_path}/*.*')
            count = 0
            for f in folder_files:
                if f.endswith('.csv') or f.endswith('.txt') or f.endswith('.tsv'):
                    try:
                        sep = '\t' if f.endswith('.txt') or f.endswith('.tsv') else ','
                        d = pd.read_csv(f, sep=sep, low_memory=False, on_bad_lines='skip', encoding_errors='ignore')
                        count += len(d)
                    except Exception as e:
                        print(f"Error reading {f}: {e}")
                        raise
            if count > 0:
                db_counts[folder] = count
    return db_counts

def process_species(species, out_dir):
    f = os.path.join(out_dir, f'merged_{species}_metabolites.csv')
    if not os.path.exists(f):
        print(f"File {f} not found.")
        return

    print(f"Processing {species} datasets...")
    df = pd.read_csv(f, low_memory=False)

    total_merged_rows = len(df)
    db_breakdown = df['database'].value_counts().to_dict()
    
    # 0. Create a dictionary of HMDB_ID and Metabolite_Name
    
    # Keep relevant columns & drop missing
    df_copy = df[['HMDB_ID', 'Metabolite_Name']].dropna(subset=['HMDB_ID', 'Metabolite_Name']).copy()
    df_copy = df_copy.drop_duplicates()
    
    # --- Step 2: Identify rows needing splitting ---
    mask = df_copy['Metabolite_Name'].str.contains(';', na=False)
    df_split = df_copy[mask].copy()       # rows to split
    df_nosplit = df_copy[~mask].copy()    # rows to keep as-is
    
    # --- Step 3: Split into lists ---
    df_split['HMDB_ID'] = df_split['HMDB_ID'].astype(str).str.split(',')
    df_split['Metabolite_Name'] = df_split['Metabolite_Name'].astype(str).str.split(';')
    
    # Remove whitespace
    df_split['HMDB_ID'] = df_split['HMDB_ID'].apply(lambda lst: [x.strip() for x in lst])
    df_split['Metabolite_Name'] = df_split['Metabolite_Name'].apply(lambda lst: [x.strip() for x in lst])
    
    # --- Step 4: Check that lengths match ---
    mismatch = df_split[df_split['HMDB_ID'].str.len() != df_split['Metabolite_Name'].str.len()]
    if not mismatch.empty:
        print("⚠️ Mismatched rows (ID count != Metabolite count):")
        print(mismatch)
    
    # --- Step 5: Pair 1-to-1 ---
    df_split['pairs'] = [list(zip(i, j)) for i, j in zip(df_split['HMDB_ID'], df_split['Metabolite_Name'])]
    
    # Explode pairs into separate rows
    df_split = df_split.explode('pairs')
    df_split[['HMDB_ID', 'Metabolite_Name']] = pd.DataFrame(df_split['pairs'].tolist(), index=df_split.index)
    df_split = df_split.drop(columns=['pairs'])
    
    # --- Step 6: Combine split and non-split rows ---
    df_final = pd.concat([df_split, df_nosplit], ignore_index=True)
    
    # --- Step 7: Remove duplicates & reset index ---
    df_final = df_final.drop_duplicates().reset_index(drop=True)

    # Correct deduplication for comma-separated HMDB_IDs
    HMDB_dict = (
        df_final.groupby('Metabolite_Name')['HMDB_ID']
        .apply(lambda ids: ','.join(sorted(set(
            [i.strip() for sublist in ids for i in sublist.split(',')]
        ))))
        .reset_index()
    )


    # 1. Filter out scCellFie_value == 0
    if 'scCellFie_value' in df.columns:
        df['scCellFie_value'] = pd.to_numeric(df['scCellFie_value'], errors='coerce')
        df_filtered = df[~(df['scCellFie_value'] == 0)].copy()
    else:
        df_filtered = df.copy()

    filtered_rows = len(df_filtered)

    out1 = os.path.join(out_dir, f'merged_{species}_metabolites_filtered_scCellfie_value_0.csv')
    df_filtered.to_csv(out1, index=False)
    # We will build Target_original using a prioritization strategy
    df_pairs = df_filtered.copy()
    df_pairs['Target_original'] = np.nan
    
    # Standard priority for all databases:
    standard_cols = ['Gene_Name', 'Receptor_Gene_Symbol', 'Ligand_Gene_Symbol', 'Target_Gene', 'Transporter', 'Enzyme', 'Uniprot', 'uniprot']
    target_cols = ['Uniprot', 'uniprot', 'Gene_Name', 'Receptor_Gene_Symbol', 'Ligand_Gene_Symbol', 'Target_Gene', 'Transporter', 'Enzyme']
    avail_targets = [c for c in target_cols if c in df_filtered.columns]

    for tc in standard_cols:
        if tc in avail_targets:
            df_pairs['Target_original'] = df_pairs['Target_original'].fillna(df_pairs[tc])
            
    # For CellPhoneDB, prefer uniprot over Gene_Name because CellPhoneDB's Gene_Name is often complex components
    if 'database' in df_pairs.columns:
        mask = df_pairs['database'].astype(str).str.contains('CellPhoneDB', case=False, na=False)
        for tc in ['Uniprot', 'uniprot']:
            if tc in avail_targets:
                df_pairs.loc[mask, 'Target_original'] = df_pairs.loc[mask, tc].fillna(df_pairs.loc[mask, 'Target_original'])

    df_pairs = df_pairs.dropna(subset=['Metabolite_Name', 'Target_original'])
    
    # 1. Explode Metabolite_Name (split by | or ;)
    df_pairs['Metabolite_Name'] = df_pairs['Metabolite_Name'].astype(str).str.split(r'\s*[;|]\s*')
    df_pairs = df_pairs.explode('Metabolite_Name')
    
    # 2. Explode Target_original (split by , or / or ; or |)
    df_pairs['Target_original'] = df_pairs['Target_original'].astype(str).str.split(r'\s*[,/;|]+\s*')
    df_pairs = df_pairs.explode('Target_original')
    
    # Clean up empty strings
    df_pairs['Metabolite_Name'] = df_pairs['Metabolite_Name'].str.strip().str.lower()
    df_pairs['Target_original'] = df_pairs['Target_original'].str.strip()
    df_pairs = df_pairs[(df_pairs['Metabolite_Name'] != '') & (df_pairs['Target_original'] != '')]

    import json
    import re
    hgnc_path = os.path.join(out_dir, '..', 'input', 'hgnc_approved_genes.json')
    alias_map = {}
    uniprot_map = {}
    uniprot_to_hgnc = {}
    hgnc_canonical = set()
    if os.path.exists(hgnc_path):
        with open(hgnc_path, "r", encoding="utf-8") as f:
            hgnc_data = json.load(f)
        for doc in hgnc_data.get("response", {}).get("docs", []):
            canon = doc.get("symbol")
            if not canon: continue
            
            hgnc_canonical.add(canon)
            alias_map[str(canon).upper()] = canon
            for alias in doc.get("alias_symbol", []):
                alias_map[str(alias).upper()] = canon
            for prev in doc.get("prev_symbol", []):
                alias_map[str(prev).upper()] = canon
                
            u_ids = doc.get("uniprot_ids", [])
            if u_ids:
                uniprot_map[canon] = str(u_ids[0])
                for uid in u_ids:
                    uniprot_to_hgnc[str(uid).upper()] = canon


    unmapped_targets = set()

    def map_hgnc_target(raw_val):
        if pd.isna(raw_val): return np.nan
        p = str(raw_val).strip()
        if not p: return np.nan
        p_up = p.upper()
        
        # Check uniprot first
        if p_up in uniprot_to_hgnc:
            return uniprot_to_hgnc[p_up]
        
        # Check alias
        if p_up in alias_map:
            return alias_map[p_up]
            
        # Unmapped case
        unmapped_targets.add(p)
        if '_' not in p:
            return p.upper()
        return p

    def map_target_uniprot(canon_val):
        if pd.isna(canon_val): return np.nan
        if canon_val in uniprot_map:
            return uniprot_map[canon_val]
        return np.nan


    print("  Canonicalizing Target symbols against HGNC mappings...")
    df_pairs['Target'] = df_pairs['Target_original'].apply(map_hgnc_target)
    df_pairs['Target_Uniprot'] = df_pairs['Target'].apply(map_target_uniprot)

    def flatten_unique(series):
        st = set(series.dropna().astype(str))
        if len(st) == 0: return np.nan
        
        if series.name == 'database':
            all_db = []
            for item in st:
                for slash_split in item.split('/'):
                    for comma_split in slash_split.split(','):
                        c = comma_split.strip()
                        if c: all_db.append(c)
            return ', '.join(sorted(list(set(all_db))))
            
        return ' | '.join(sorted(list(st)))

    group_keys = ['Metabolite_Name', 'Target_original']
    agg_dict = {c: flatten_unique for c in df_pairs.columns if c not in group_keys}

    df_grouped = df_pairs.groupby(group_keys, as_index=False, dropna=False).agg(agg_dict)

    # Normalize metabolite names
    df_grouped['Metabolite_Name'] = df_grouped['Metabolite_Name'].str.strip()
    df_grouped = df_grouped.reset_index(drop=True)

    

    # Clean up duplicated target columns
    for tc in avail_targets:
        if tc in df_grouped.columns:
            df_grouped[tc] = df_grouped.apply(lambda row: np.nan if row['Target'] == row[tc] or not pd.isna(row['Target']) else row[tc], axis=1)


    # Add databases_count
    def count_dbs(db_str):
        if pd.isna(db_str): return 0
        return len([x for x in str(db_str).split(',') if x.strip()])
    
    df_grouped['databases_count'] = df_grouped['database'].apply(count_dbs)

    # --- Apply the exact same HMDB mapping as unique_metabs to target pairs ---
    # This prevents fragmentation (e.g. fumarate mapping to 131,134 here but 134 in unique_metabs)
    # and fills all missing HMDB_IDs perfectly.
    mapped_hmdb_pairs = df_grouped['Metabolite_Name'].str.lower().map(
        HMDB_dict.assign(Metabolite_Name=HMDB_dict['Metabolite_Name'].str.lower())
                 .drop_duplicates(subset=['Metabolite_Name'])
                 .set_index('Metabolite_Name')['HMDB_ID']
    )
    df_grouped['HMDB_ID'] = mapped_hmdb_pairs

    out2 = os.path.join(out_dir, f'{species}_database_merge_unique_metab_target_pairs.csv')
    df_grouped.to_csv(out2, index=False)
    unique_pairs = len(df_grouped)
    print(f"  -> Saved {out2} ({unique_pairs} rows)")
    
    # Save unmapped targets to a CSV
    unmapped_targets_file = os.path.join(out_dir, f'{species}_unmapped_targets.csv')
    pd.DataFrame({'Unmapped_Target': sorted(list(unmapped_targets))}).to_csv(unmapped_targets_file, index=False)
    print(f"  -> Saved unmapped HGNC targets to {unmapped_targets_file} ({len(unmapped_targets)} targets)")


    # 3. Unique Metabolites dictionary
    df_db = df_filtered.dropna(subset=['Metabolite_Name']).copy()
    
    # Normalize metabolite names
    df_db['Metabolite_Name'] = df_db['Metabolite_Name'].str.strip().str.lower()
    df_db['Metabolite_Name'] = df_db['Metabolite_Name'].str.split(r'\s*[;|]\s*')
    df_db = df_db.explode('Metabolite_Name')
    df_db = df_db.reset_index(drop=True)
    
    def extract_dbs(series):
        st = set(series.dropna().astype(str))
        if not st: return np.nan
        all_db = []
        for item in st:
            for slash_split in item.split('/'):
                for comma_split in slash_split.split(','):
                    c = comma_split.strip()
                    if c: all_db.append(c)
        return ', '.join(sorted(list(set(all_db))))
    
    def extract_hmdbs(series):
        all_ids = []
        for v in series.dropna().astype(str):
            for c in v.replace(';', ',').replace('|', ',').split(','):
                c = c.strip()
                if c and c.upper().startswith('HMDB'):
                    all_ids.append(c.upper())
        return ','.join(sorted(list(set(all_ids)))) if all_ids else np.nan

    agg_dict = {'database': extract_dbs}
    if 'HMDB_ID' in df_db.columns:
        agg_dict['HMDB_ID'] = extract_hmdbs

    df_met_db = df_db.groupby('Metabolite_Name', as_index=False, dropna=False).agg(agg_dict)
    
    df_met_db['databases_count'] = df_met_db['database'].apply(count_dbs)
    
    # One-liner mapping (case-insensitive) - strictly 1:1 mapping!
    df_met_db['HMDB_ID'] = df_met_db['Metabolite_Name'].str.lower().map(
        HMDB_dict.assign(Metabolite_Name=HMDB_dict['Metabolite_Name'].str.lower())
                 .drop_duplicates(subset=['Metabolite_Name'])
                 .set_index('Metabolite_Name')['HMDB_ID']
    )

    
    out3 = os.path.join(out_dir, f'{species}_database_merge_unique_metab.csv')
    df_met_db.to_csv(out3, index=False)
    unique_metabs = len(df_met_db)
    print(f"  -> Saved {out3} ({unique_metabs} rows)")
    
    unmapped_hmdb = df_met_db[df_met_db['HMDB_ID'].isna()]['Metabolite_Name'].dropna().unique()
    unmapped_hmdb_file = os.path.join(out_dir, f'{species}_unmapped_hmdb.csv')
    pd.DataFrame({'Unmapped_Metabolite': sorted(list(unmapped_hmdb))}).to_csv(unmapped_hmdb_file, index=False)
    print(f"  -> Saved unmapped HMDB metabolites to {unmapped_hmdb_file} ({len(unmapped_hmdb)} metabolites)")

    # 4. Generate the dedicated Statistics File
    raw_counts = compute_raw_counts()
    stat_file = os.path.join(out_dir, f'{species}_database_merge_metab_statistics.txt')
    with open(stat_file, 'w') as sf:
        sf.write(f"Database Merge Statistics Report - {species.capitalize()} Metabolites\n")
        sf.write("="*60 + "\n\n")

        sf.write("--- Final Output Dataset Totals ---\n")
        sf.write(f"* Total Rows (Pre-filtering): {total_merged_rows}\n")
        sf.write(f"* Final Merged Rows (scCellFie_value != 0): {filtered_rows}\n")
        sf.write(f"* Unique Metabolites: {unique_metabs}\n")
        sf.write(f"* Unique Metabolite-Target Pairs: {unique_pairs}\n\n")

        sf.write("--- Row Counts by Database (Pre-filtering) ---\n")
        for k, v in db_breakdown.items():
            sf.write(f"* {k}: {v} rows\n")
            
        sf.write("\n--- Original Raw Row Counts (Before Deduplication & Merging) ---\n")
        sf.write("Note: These generic file counts contain proteins and overlapping species.\n")
        for k, v in raw_counts.items():
            sf.write(f"* {k}: {v} rows\n")
            
    print(f"  -> Saved {stat_file}\n")


if __name__ == '__main__':
    # Dynamic path resolution relative to this script
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
    base_dir = os.path.join(PROJECT_ROOT, 'input', 'databases')
    out_dir = os.path.join(PROJECT_ROOT, 'output')
    
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(base_dir)
    
    process_species('human', out_dir)
    process_species('mouse', out_dir)
    
    print("Cleanup: Removing obsolete temporary test files...")
    for old_file in glob.glob(os.path.join(out_dir, 'final_*')) + [os.path.join(out_dir, 'database_merge_statistics.txt')]:
        if os.path.exists(old_file):
            os.remove(old_file)
            print(f"Removed {old_file}")
