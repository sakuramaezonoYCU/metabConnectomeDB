import os
import pandas as pd
import glob
import re
import csv
import numpy as np
from collections import defaultdict

base_dir = '/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/input/databases'
out_dir = os.path.abspath(os.path.join(base_dir, '../../output'))
os.makedirs(out_dir, exist_ok=True)
os.chdir(base_dir)

# Full HMDB annotations (primary source — 248k entries)
df_hmdb = pd.read_csv('HMDB_metabolites', sep=',')
# Create lookup dictionary
hmdb_dict = dict(zip(df_hmdb['HMDB_ID'], df_hmdb['NAME']))


def get_species(filename):
    fname = os.path.basename(filename).lower()
    if 'mouse' in fname or 'mus musculus' in fname:
        return 'mouse'
    elif 'zebrafish' in fname:
        return 'zebrafish'
    else:
        return 'human'

def normalize_empty(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    obj_cols = df.select_dtypes(include=["object", "string"]).columns

    # Safely convert to string before using .str
    df[obj_cols] = df[obj_cols].apply(
        lambda s: s.astype("string").str.strip().str.upper()
    )

    df[obj_cols] = df[obj_cols].replace(
        {"": np.nan, "NA": np.nan, "N/A": np.nan, "NULL": np.nan}
    )

    return df
    return df

# For CellPhoneDB processing
# Skipping complex_input.csv as per new requirement; will use extract_NPL_ETR_from_cellPhoneDB5.csv instead.
cpdb_extract_df = None
if os.path.exists('CellPhoneDBv5/extract_NPL_ETR_from_cellPhoneDB5.csv'):
    cpdb_extract_df = pd.read_csv('CellPhoneDBv5/extract_NPL_ETR_from_cellPhoneDB5.csv')

all_human_metabolites = []
all_human_proteins = []
all_mouse_metabolites = []
all_mouse_proteins = []

for root, dirs, files in os.walk(base_dir):
    folder = os.path.basename(root)
    if folder == 'databases': 
        continue
    # Load Ligand_list once per MetaLigand folder visit — annotation table with one row
    # per metabolite (transporter, class info etc.). Pre-deduplicated on Metabolite_Name
    # to prevent Cartesian explosion when merged against NPLRdb_* interaction tables.
    # One-to-many (metabolite -> many receptors) is preserved because the fan-out
    # lives in NPLRdb_*, not in Ligand_list.
    metaligand_ligand_list = None
    if folder == 'MetaLigand':
        ligand_list_path = os.path.join(root, 'Ligand_list.csv')
        if os.path.exists(ligand_list_path):
            ll = pd.read_csv(ligand_list_path)
            ll = ll.rename(columns={'Compounds': 'Metabolite_Name'})
            # Define merge_unique_genes here so it's in scope for LR_data_* merge below
            def merge_unique_genes(row):
                genes = []
                for col in ['Transporter_genes', 'Transporter_genes_in']:
                    val = row.get(col)
                    if pd.notna(val):
                        genes.extend([g.strip() for g in str(val).split(',') if g.strip()])
                seen = set()
                unique_genes = []
                for g in genes:
                    if g not in seen:
                        seen.add(g)
                        unique_genes.append(g)
                return ','.join(unique_genes) if unique_genes else np.nan
            # Drop transporter and receptor columns — LR_data_* is the authoritative
            # source for these. Ligand_list only contributes Class/SuperClass/Synonyms etc.
            ll = ll.drop(columns=[c for c in [
                'Receptor_genes', 'Receptor genes',
                'Transporter_genes', 'Transporter_genes_in'
            ] if c in ll.columns])
            # Deduplicate on Metabolite_Name: Ligand_list is annotation-only,
            # so one row per metabolite is correct here
            metaligand_ligand_list = (
                ll.drop_duplicates(subset=['Metabolite_Name'])
                if 'Metabolite_Name' in ll.columns
                else ll.drop_duplicates()
            )
    for f in files:
        # Skip previously generated output files to prevent recursion loop.
        # Exempt MEBOCOST base files whose names legitimately start with species prefix.
        is_mebocost_base = (
            os.path.basename(root) == 'MEBOCOST'
            and bool(re.search(r'(human|mouse)_met_sensor_update', f))
        )
        if not is_mebocost_base:
            if f.startswith('human_') or f.startswith('mouse_') or f.startswith('merged_'):
                continue
            
            if not (f.endswith('.csv') or f.endswith('.txt') or f.endswith('.tsv')):
                continue
            
        filepath = os.path.join(root, f)
        species = get_species(filepath)
        if species == 'zebrafish':
            continue  # Skip zebrafish
            
        db_name = folder
        
        sep = '\t' if filepath.endswith('.txt') or filepath.endswith('.tsv') else ','
        # Supplementary_Table_1.csv is comma-separated despite the name — sep stays ','
        
        try:
            if folder == 'NeuronChat':
                with open(filepath, 'r', encoding='utf-8') as f_nc:
                    content = f_nc.read()
                content = content.replace('Synaptic adhesion molecule', 'Synaptic_adhesion_molecule')
                content = content.replace('Gap junction protein', 'Gap_junction_protein')
                content = content.replace('ligand degradation', 'ligand_degradation')
                content = content.replace('ligand uptake', 'ligand_uptake')
                content = content.replace('ligand activation', 'ligand_activation')
                content = content.replace('ligand-receptor', 'ligand_receptor')
                
                import io
                content = content.strip()
                df = pd.read_csv(io.StringIO(content), sep=r'\s+', header=None, skiprows=1)
                df.columns = ['interaction_name', 'ligand_type', 'interaction_type', 'lig_contributor', 'lig_contributor_group', 'lig_contributor_coeff', 'target_subunit', 'target_subunit_group', 'target_subunit_coeff']
                
                if 'lig_contributor' in df.columns:
                    df = df.rename(columns={'lig_contributor': 'Metabolite_Name'})
                if 'interaction_name' in df.columns:
                    # e.g. "VIP_VIPR1" -> VIPR1
                    df['Receptor_Gene_Symbol'] = df['interaction_name'].apply(lambda x: str(x).split('_')[1] if '_' in str(x) else np.nan)
            
            # scCellFie requires special multi-table merging, so load it later as a single unit per species
            elif folder == 'scCellFie':
                continue
            
            else:
                try:
                    df = pd.read_csv(filepath, sep=sep, encoding='utf-8', on_bad_lines='skip', low_memory=False)
                except UnicodeDecodeError:
                    df = pd.read_csv(filepath, sep=sep, encoding='latin1', on_bad_lines='skip', low_memory=False)
        except Exception as e:
            print(f"Failed to read {filepath}: {e}")
            continue

        if df is None or df.empty:
            continue

        df['database'] = db_name
        is_metab = False
        is_prot = False

        # ---- Cellinker2 & MRCLinkDB ----
        # Both folders share the same file structure. Annotation files (enzyme,
        # transporter, Metabolite-cell interaction) are skipped here and merged
        # inside the base branch so that both Cellinker2 AND MRCLinkDB get the
        # same enrichment before deduplicate_cellinker_mrclink combines them.
        if folder in ['Cellinker2', 'MRCLinkDB']:
            if 'protein L-R interaction' in f:
                is_prot = True

            elif 'metabolite L-R interaction' in f:
                # ── Base metabolite file ─────────────────────────────────────
                df = df.rename(columns={
                    'Metabolite name':     'Metabolite_Name',
                    'HMDB ID':             'HMDB_ID',
                    'Receptor_gene ID':    'Receptor_Gene_ID',
                    'Receptor gene ID':    'Receptor_Gene_ID',
                    'Receptor uniprot_ id':'Receptor_Uniprot',
                    'Receptor uniprot id': 'Receptor_Uniprot',
                    'Receptor_symbol':     'Receptor_Gene_Symbol',
                    'Receptor symbol':     'Receptor_Gene_Symbol',
                    'protein name':        'Protein_Name',
                })

                species_prefix_f = f.split(' metabolite')[0]  # e.g. "Homo sapiens"

                # ── Helper: read annotation file robustly ────────────────────
                def _read_anno(path):
                    try:
                        return pd.read_csv(path, sep='	', encoding='utf-8', on_bad_lines='skip')
                    except UnicodeDecodeError:
                        return pd.read_csv(path, sep='	', encoding='latin1', on_bad_lines='skip')

                # ── Enzyme annotation: GENE_NAME (or Mouse_gene symbol) -> Enzyme
                enzyme_path = os.path.join(root, f"{species_prefix_f} enzyme.txt")
                if os.path.exists(enzyme_path):
                    enz = _read_anno(enzyme_path)
                    enz = enz.rename(columns={
                        'METABOLITE_NAME':   'Metabolite_Name',
                        'mETABOLITE_NAME':   'Metabolite_Name',
                        'Mouse_gene symbol': 'GENE_NAME',  # normalise mouse column
                    })
                    if 'type' in enz.columns:
                        gene_col = 'GENE_NAME' if 'GENE_NAME' in enz.columns else None
                        if gene_col:
                            enz['Enzyme'] = enz.apply(
                                lambda r: r[gene_col] if str(r['type']).strip() == 'Enzyme' else np.nan, axis=1
                            )
                            enz = enz.drop(columns=['type', gene_col])
                    cols_to_merge = [c for c in enz.columns
                                     if c in ['HMDB_ID', 'Metabolite_Name'] or c not in df.columns]
                    df = df.merge(
                        enz[cols_to_merge],
                        on=['HMDB_ID', 'Metabolite_Name'], how='left'
                    )

                # ── Transporter annotation: GENE_NAME (or Mouse_gene symbol) -> Transporter
                transporter_path = os.path.join(root, f"{species_prefix_f} transporter protein.txt")
                if os.path.exists(transporter_path):
                    trp = _read_anno(transporter_path)
                    trp = trp.rename(columns={
                        'METABOLITE_NAME':   'Metabolite_Name',
                        'mETABOLITE_NAME':   'Metabolite_Name',
                        'Mouse_gene symbol': 'GENE_NAME',  # normalise mouse column
                    })
                    if 'type' in trp.columns:
                        gene_col = 'GENE_NAME' if 'GENE_NAME' in trp.columns else None
                        if gene_col:
                            trp['Transporter'] = trp.apply(
                                lambda r: r[gene_col] if str(r['type']).strip() == 'Transporter' else np.nan, axis=1
                            )
                            trp = trp.drop(columns=['type', gene_col])
                    cols_to_merge = [c for c in trp.columns
                                     if c in ['HMDB_ID', 'Metabolite_Name'] or c not in df.columns]
                    df = df.merge(
                        trp[cols_to_merge],
                        on=['HMDB_ID', 'Metabolite_Name'], how='left'
                    )

                # ── Metabolite-cell interaction (only exists in MRCLinkDB) ───
                cell_path = os.path.join(root, 'Metabolite-cell interaction.txt')
                if os.path.exists(cell_path):
                    cell = _read_anno(cell_path)
                    cell = cell.rename(columns={
                        'HMDB ID':       'HMDB_ID',
                        'Metabolite name':'Metabolite_Name',
                    })
                    cols_to_merge = [c for c in cell.columns
                                     if c in ['HMDB_ID', 'Metabolite_Name'] or c not in df.columns]
                    df = df.merge(
                        cell[cols_to_merge],
                        on=['HMDB_ID', 'Metabolite_Name'], how='left'
                    )

                df = df.dropna(axis=1, how='all')
                df['database'] = db_name
                is_metab = True

            elif 'enzyme' in f or 'transporter protein' in f or 'Metabolite-cell interaction' in f:
                # Merged inside the base branch above — skip independent processing
                continue
                
        # ---- MetaLigand ----
        # NPLRdb_* is the base interaction table: one row per metabolite-receptor pair
        # (L -> Metabolite_Name, R -> Receptor_Gene_Symbol). One-to-many relationships
        # are preserved here — the same metabolite can appear across multiple rows with
        # different receptors.
        # Ligand_list.csv is annotation only (one row per metabolite) — loaded once
        # before this loop into metaligand_ligand_list and merged into all NPLRdb_*
        # species files. Skip its own loop iteration.
        # ---- MetaLigand ----
        # NPLRdb_* is the base interaction table: one row per metabolite-receptor pair
        # (L -> Metabolite_Name, R -> Receptor_Gene_Symbol).
        # LR_data_* is the primary annotation table for NPLRdb_* (one row per metabolite).
        # Ligand_list.csv is secondary annotation — loaded once before this loop.
        # All three are merged per species. One-to-many (metabolite -> many receptors)
        # is preserved because the fan-out lives in NPLRdb_*, not in the annotation tables.
        elif folder == 'MetaLigand':
            if f == 'Ligand_list.csv':
                continue

            elif f.startswith('NPLRdb_'):
                if species == 'zebrafish':
                    continue
                # Rename both human (L/R) and mouse (Compounds/Receptor_genes) column names
                df = df.rename(columns={
                    'L':              'Metabolite_Name',
                    'Compounds':      'Metabolite_Name',
                    'R':              'Receptor_Gene_Symbol',
                    'Receptor_genes': 'Receptor_Gene_Symbol',
                })
                # Drop index column if present
                df = df.drop(columns=[c for c in df.columns if str(c).startswith('Unnamed')])

                # Merge primary annotation (LR_data_*) — same species
                lr_data_path = os.path.join(root, f.replace('NPLRdb_', 'LR_data_'))
                if os.path.exists(lr_data_path):
                    try:
                        lr = pd.read_csv(lr_data_path, encoding='utf-8', index_col=0, on_bad_lines='skip', low_memory=False)
                    except UnicodeDecodeError:
                        lr = pd.read_csv(lr_data_path, encoding='latin1', index_col=0, on_bad_lines='skip', low_memory=False)
                    lr = lr.rename(columns={'Compounds': 'Metabolite_Name'})
                    lr = lr.drop(columns=[c for c in ['Receptor_genes', 'Receptor genes'] if c in lr.columns])
                    
                    # Combine Transporter_genes + Transporter_genes_in into unified Transporter column
                    if 'Transporter_genes' in lr.columns or 'Transporter_genes_in' in lr.columns:
                        lr['Transporter'] = lr.apply(merge_unique_genes, axis=1)
                        lr = lr.drop(columns=[c for c in ['Transporter_genes', 'Transporter_genes_in'] if c in lr.columns])
                    
                    lr_dedup = lr.drop_duplicates(subset=['Metabolite_Name'])
                    cols_to_merge = ['Metabolite_Name'] + [
                        c for c in lr_dedup.columns
                        if c != 'Metabolite_Name' and c not in df.columns
                    ]
                    df = df.merge(lr_dedup[cols_to_merge], on='Metabolite_Name', how='left')

                # Merge secondary annotation (Ligand_list) for any remaining gaps
                if metaligand_ligand_list is not None:
                    cols_to_merge = ['Metabolite_Name'] + [
                        c for c in metaligand_ligand_list.columns
                        if c != 'Metabolite_Name' and c not in df.columns
                    ]
                    df = df.merge(metaligand_ligand_list[cols_to_merge], on='Metabolite_Name', how='left')

                df = df.dropna(axis=1, how='all')
                df['database'] = db_name
                is_metab = True

            # LR_data_* — merged inside NPLRdb_* branch above, skip standalone
            elif f.startswith('LR_data_'):
                continue

            # PLRdb_* — zebrafish only in current files, skip
            elif f.startswith('PLRdb_'):
                continue

            else:
                print(f"[MetaLigand WARNING] Unrecognised file skipped: {f!r}")
                continue
        
        # ---- MEBOCOST ----
        elif folder == 'MEBOCOST':
            # Skip xlsx entirely — use Supplementary_Table_1.csv instead
            if f.endswith('.xlsx'):
                continue

            # ── Base file: {species}_met_sensor_update_{anything}.tsv ──────
            if re.search(r'(human|mouse)_met_sensor_update', f):
                df = df.rename(columns={
                    'standard_metName': 'Metabolite_Name',
                    'Gene_name':        'Sensor_Gene',
                })
                # 'Annotation' is richer than Sensor_Type — promote it
                if 'Annotation' in df.columns:
                    df = df.rename(columns={'Annotation': 'Sensor_Type'})
                # Drop metName — standard_metName (now Metabolite_Name) takes priority
                df = df.drop(columns=[c for c in ['metName'] if c in df.columns])

                # Merge Supplementary_Table_1.csv on [Metabolite_Name, Sensor_Gene]
                supp_path = os.path.join(root, 'Supplementary_Table_1.csv')
                if os.path.exists(supp_path):
                    supp = pd.read_csv(supp_path, sep=',', encoding='latin1')
                    # Drop Sensor_Type from supp — base Annotation takes priority
                    supp = supp.drop(columns=[c for c in ['Sensor_Type'] if c in supp.columns])
                    cols_to_merge = [c for c in supp.columns
                                     if c in ['Metabolite_Name', 'Sensor_Gene']
                                     or c not in df.columns]
                    supp_dedup = supp[cols_to_merge].drop_duplicates(
                        subset=['Metabolite_Name', 'Sensor_Gene']
                    )
                    df = df.merge(supp_dedup, on=['Metabolite_Name', 'Sensor_Gene'], how='left')

                # Merge HMDB annotation files on Metabolite_Name
                for hmdb_fname in [
                    'hmdb_blood_metabolite_concentration.tsv',
                    'metabolite_annotation_HMDB_summary.tsv',
                ]:
                    hmdb_path = os.path.join(root, hmdb_fname)
                    if os.path.exists(hmdb_path):
                        try:
                            hmdb_df = pd.read_csv(hmdb_path, sep='\t', encoding='utf-8', on_bad_lines='skip')
                        except UnicodeDecodeError:
                            hmdb_df = pd.read_csv(hmdb_path, sep='\t', encoding='latin1', on_bad_lines='skip')
                        # Rename 'metabolite' column to Metabolite_Name for joining
                        hmdb_df = hmdb_df.rename(columns={'metabolite': 'Metabolite_Name'})
                        # Only bring in columns not already in base
                        cols_to_merge = [c for c in hmdb_df.columns
                                         if c == 'Metabolite_Name' or c not in df.columns]
                        hmdb_dedup = hmdb_df[cols_to_merge].drop_duplicates(subset=['Metabolite_Name'])
                        df = df.merge(hmdb_dedup, on='Metabolite_Name', how='left')

                df = df.dropna(axis=1, how='all')
                df['database'] = db_name

                # ── Assign Receptor_Gene_Symbol / Transporter / Enzyme ──────
                # Sensor_Type partial matching:
                #   'Receptor' / 'Nuclear Receptor' → Receptor_Gene_Symbol
                #   'Transporter' / 'Channel'       → Transporter
                #   'Enzyme'                         → Enzyme
                # Multi-type rows (e.g. 'Receptor; Transporter') fill multiple columns
                def assign_sensor_cols(row):
                    stype = str(row.get('Sensor_Type', '')).strip()
                    gene  = row.get('Sensor_Gene', np.nan)
                    receptor    = np.nan
                    transporter = np.nan
                    enzyme      = np.nan
                    if not stype or stype == 'nan':
                        return pd.Series([receptor, transporter, enzyme])
                    if re.search(r'(?i)receptor', stype):
                        receptor = gene
                    if re.search(r'(?i)transporter|channel', stype):
                        transporter = gene
                    if re.search(r'(?i)enzyme', stype):
                        enzyme = gene
                    return pd.Series([receptor, transporter, enzyme])

                df[['Receptor_Gene_Symbol', 'Transporter', 'Enzyme']] = df.apply(
                    assign_sensor_cols, axis=1
                )
                is_metab = True

            # ── All other MEBOCOST files are merged inside the base branch ──
            else:
                continue

        # ---- NeuronChat ----
        elif folder == 'NeuronChat':
            if 'ligand_type' in df.columns:
                metab_types = ['Gas', 'Neurotransmitter'] 
                df_metab = df[df['ligand_type'].isin(metab_types)].copy()
                df_prot = df[~df['ligand_type'].isin(metab_types)].copy()
                
                if species == 'human':
                    if not df_metab.empty: all_human_metabolites.append(df_metab)
                    if not df_prot.empty: all_human_proteins.append(df_prot)
                else:
                    if not df_metab.empty: all_mouse_metabolites.append(df_metab)
                    if not df_prot.empty: all_mouse_proteins.append(df_prot)
            continue

        # ---- CellPhoneDBv5 ----
        elif folder == 'CellPhoneDBv5':
            # Pivot Uniprot columns for any CellPhoneDB file containing them
            uniprot_cols = [c for c in df.columns if re.match(r'^uniprot_\d+$', c)]
            if uniprot_cols:
                id_vars = [c for c in df.columns if c not in uniprot_cols]
                df = df.melt(id_vars=id_vars, value_vars=uniprot_cols, value_name='uniprot')
                df = df.dropna(subset=['uniprot'])
                df = df[df['uniprot'].astype(str).str.strip() != '']
                if 'variable' in df.columns:
                     df = df.drop(columns=['variable'])
                # Rename complex_name -> Task to align with the shared column name across databases
                if 'complex_name' in df.columns:
                    df = df.rename(columns={'complex_name': 'Task'})
                # Merge extra cols from extract_NPL_ETR on Task (was complex_name)
                if cpdb_extract_df is not None and 'Task' in df.columns:
                    extract_renamed = cpdb_extract_df.rename(columns={'complex_name': 'Task'})
                    cols_to_merge = [c for c in extract_renamed.columns
                                         if c == 'Task' or c not in df.columns]
                    df = df.merge(extract_renamed[cols_to_merge], on='Task', how='left')
                # Prune columns that are entirely empty after the merge
                df = df.dropna(axis=1, how='all')
                is_metab = True
                
            elif f == 'interaction_input.csv':
                if 'directionality' in df.columns:
                    df = df[df['directionality'] == 'Ligand_Receptor'].copy()
                    df = df.drop(columns=['directionality'])
                if 'annotation_strategy' in df.columns:
                    df = df.drop(columns=['annotation_strategy'])
                
                rows = []
                for _, row in df.iterrows():
                    interactors = str(row.get('interactors', ''))
                    if '-' in interactors:
                        lig, recs = interactors.split('-', 1)
                        for rec in recs.split('+'):
                            r = row.copy()
                            r['Ligand_Gene_Symbol'] = lig
                            r['Receptor_Gene_Symbol'] = rec
                            rows.append(r)
                    else:
                        rows.append(row)
                        
                if rows:
                    df = pd.DataFrame(rows)
                if 'interactors' in df.columns:
                    df = df.drop(columns=['interactors'])
                    
                is_prot = True
            else:
                continue
                
            if is_metab:
                # Split on Metabolite_Name presence: matched rows -> metabolites,
                # unmatched (no Metabolite_Name after extract_NPL_ETR merge) -> proteins
                if 'Metabolite_Name' in df.columns:
                    mask = df['Metabolite_Name'].notna() & df['Metabolite_Name'].astype(str).str.strip().ne('')
                    df_metab = df[mask].copy()
                    def add_hyphen(name):
                        # Add hyphen after numbers not followed by "-", end of string, or "(s)"
                        name = re.sub(r'(\d+)(?!(?:-|$|\(s\)))', r'\1-', name)
                        # Add hyphen after "oxo" with same exceptions
                        name = re.sub(r'(oxo)(?!(?:-|$|\(s\)))', r'\1-', name)
                        return name
                    # Apply function to the column
                    df_metab['Metabolite_Name'] = df_metab['Metabolite_Name'].apply(add_hyphen)
                    df_metab['Metabolite_Name'] = df_metab['Metabolite_Name'].str.replace(
                        r'(?<=[a-z])([A-Z])',
                        r' \1',
                        regex=True
                    )
                    df_prot  = df[~mask].copy()
                else:
                    df_metab = df.copy()
                    df_prot  = pd.DataFrame(columns=df.columns)
                if not df_metab.empty:
                    # 1) Sensor_type from secreted_highlight
                    if 'secreted_highlight' in df_metab.columns:
                        df_metab["Sensor_type"] = df_metab["secreted_highlight"].apply(
                            lambda x: "Cytokines, hormones, growth factors and other immune-related proteins" if x is True else np.nan
                        )
                        df_metab = df_metab.drop(columns=["secreted_highlight"])
                
                    # 3) Rename "other"
                    if 'other' in df_metab.columns:
                        df_metab = df_metab.rename(columns={
                            "other": "cellphoneDB5_CCC_excluded"
                        })
                
                    # 4) Receptor_location
                    if 'secreted' in df_metab.columns or 'transmembrane' in df_metab.columns:
                        df_metab["Receptor_location"] = np.select(
                            [
                                df_metab.get("secreted", False) == True,
                                df_metab.get("transmembrane", False) == True
                            ],
                            [
                                "secreted",
                                "transmembrane"
                            ],
                            default=None
                        )
                        df_metab = df_metab.drop(columns=[
                            c for c in ["secreted", "transmembrane"]
                            if c in df_metab.columns
                        ])
                    if species == 'human': all_human_metabolites.append(df_metab)
                    else: all_mouse_metabolites.append(df_metab)
                if not df_prot.empty:
                    if species == 'human': all_human_proteins.append(df_prot)
                    else: all_mouse_proteins.append(df_prot)
                
            elif is_prot:
                if species == 'human': all_human_proteins.append(df)
                else: all_mouse_proteins.append(df)
            continue

        if is_metab:
            # Filter rows to ensure Metabolite_Name is present
            if 'Metabolite_Name' in df.columns:
                df = df[df['Metabolite_Name'].notna() & df['Metabolite_Name'].astype(str).str.strip().ne('')]
            if species == 'human':
                all_human_metabolites.append(df)
            else:
                all_mouse_metabolites.append(df)
        elif is_prot:
            if species == 'human':
                all_human_proteins.append(df)
            else:
                all_mouse_proteins.append(df)

# ---- Handle scCellFie independently ----
def process_sccellfie(species_str):
    task_by_gene_path = f'scCellFie/Task_by_Gene_{species_str}.csv'
    task_info_path = f'scCellFie/Task-Info_{species_str}.csv'
    thresholds_path = f'scCellFie/Thresholds_{species_str}.csv'
    
    if os.path.exists(task_by_gene_path) and os.path.exists(task_info_path) and os.path.exists(thresholds_path):
        df_tasks = pd.read_csv(task_by_gene_path, on_bad_lines='skip')
        df_info = pd.read_csv(task_info_path, on_bad_lines='skip')
        df_thresh = pd.read_csv(thresholds_path, on_bad_lines='skip')
        
        gene_cols = [c for c in df_tasks.columns if c != 'Task']
        df_melt = df_tasks.melt(id_vars=['Task'], value_vars=gene_cols, var_name='Gene_Name', value_name='scCellFie_value')
        
        # Merge info
        df_melt = df_melt.merge(df_info, on='Task', how='left')
        
        # Merge thresholds
        df_thresh = df_thresh.rename(columns={'symbol': 'Gene_Name'})
        df_melt = df_melt.merge(df_thresh, on='Gene_Name', how='left')
        
        def parse_sccellfie_task(row):
            task = str(row.get('Task', ''))
            subs = str(row.get('Subsystem', ''))
            
            meta = None
            enzyme = None
            if task == 'Co-translational translocation':
                meta = 'GTP'
            elif ' of ' in task and ' to ' in task:
                # XXX of YYY to ZZZ -> Metabolite=YYY;ZZZ
                m = re.search(r' of (.*?) to (.*)', task)
                if m:
                    meta = m.group(1).strip() + ";" + m.group(2).strip()
            elif '(to ' in task:
                # "X ... (to Y)"
                m = re.search(r'^(.*?)\s+(?:[\w\s]+?)\s*\(to\s+(.*?)\)$', task)
                if m:
                    meta = m.group(1).strip() + ";" + m.group(2).strip()
            elif ' of ' in task and ' from ' in task:
                # XXX of YYY from ZZZ -> Metabolite=ZZZ;YYY
                m = re.search(r' of (.*?) from (.*)', task)
                if m:
                    meta = m.group(2).strip() + ";" + m.group(1).strip()
            elif ' of ' in task and '(' in task and ')' in task and ('ynthesis' in task):
                # [Bio]synthesis of ZZZ(YYY) -> Metabolite=ZZZ, enzyme=YYY
                m = re.search(r' of (.*?)\((.*?)\)', task)
                if m:
                    meta = m.group(1).strip()
                    enzyme = m.group(2).strip()
            elif re.search(r'(?i)(bio)?synthesis$', task):
                # XXX synthesis / XXX biosynthesis (plain, no parentheses) -> meta=XXX
                m = re.search(r'^(.*?)\s+(?:bio)?synthesis$', task, re.IGNORECASE)
                if m:
                    meta = m.group(1).strip()
            elif 'Branching' in task:
                #  Branching (XXX) -> enzyme=XXX
                m = re.search(r'Branching[\s\(]+(.*?)[\)\s]*$', task)
                if m:
                    enzyme = m.group(1).strip().rstrip(')')
            if pd.isna(meta) or not str(meta).strip():
                if pd.notna(subs) and str(subs).strip():
                    s = str(subs).strip()
                    # remove last word if multiple words
                    if ' ' in s:
                        s = ' '.join(s.split()[:-1])
                    # normalize case: lowercase then capitalize first letter
                    #meta = s.lower().capitalize()
                else:
                    meta = None

                
            return pd.Series([meta, enzyme])  
        df_melt[['Metabolite_Name', 'enzyme']] = df_melt.apply(parse_sccellfie_task, axis=1)
        # Using a single replace call with a regex pattern that handles both replacements
        df_melt['Metabolite_Name'] = df_melt['Metabolite_Name'].str.replace(r"(processing|conversion|degradation)\s*", "", regex=True)
        # to make estrone (e1) & estrone sulfate (e1s) & estradiol-17beta (e2) same with others
        df_melt['Metabolite_Name'] = df_melt['Metabolite_Name'].str.replace(
            r"(estrone|estrone sulfate|estradiol-17beta)\s*\([^)]*\)",
            r"\1",
            regex=True
        )
        # make it all the lowercase
        df_melt['Metabolite_Name'] = df_melt['Metabolite_Name'].str.lower()
        df_melt['enzyme'] = df_melt['enzyme'].str.replace(r"(link with|metabolism)\s*", "", regex=True)
        df_melt['database'] = 'scCellFie'
        return df_melt
    return None

scc_human = process_sccellfie('human')
if scc_human is not None:
    all_human_metabolites.append(scc_human)
scc_mouse = process_sccellfie('mouse')
if scc_mouse is not None:
    all_mouse_metabolites.append(scc_mouse)

def remove_index_cols(df):
    drop_cols = []
    for c in df.columns:
        cl = str(c).lower()
        if cl.startswith('unnamed') or cl == 'index' or cl == 'id' or c == '':
            drop_cols.append(c)
        elif 'calcitriol' in cl or ('.' in str(c) and 'e-05' in str(c)): 
            # Hard-catch any corrupted headers like "Calcitriol,7.99e-5" or float strings that bled into headers
            try:
                float(str(c))
                drop_cols.append(c)
            except ValueError:
                if 'calcitriol' in cl: drop_cols.append(c)
                
    return df.drop(columns=drop_cols)


# ── The explicit rename_map for genuinely different column names ────────────
rename_map = {
    'Metabolite name': 'Metabolite_Name',
    'mETABOLITE_NAME': 'Metabolite_Name',
    'METABOLITE_NAME': 'Metabolite_Name',
    'Compounds':       'Metabolite_Name',
    'standard_metName': 'Metabolite_Name',
    'metabolite':      'Metabolite_Name',
    'L':               'Metabolite_Name',

    'HMDB ID':  'HMDB_ID',
    'HMDB_id':  'HMDB_ID',

    'Receptor gene ID':  'Receptor_Gene_ID',
    'Receptor_gene ID':  'Receptor_Gene_ID',
    'Receptor_geneid':   'Receptor_Gene_ID',

    'Receptor symbol':   'Receptor_Gene_Symbol',
    'Receptor_symbol':   'Receptor_Gene_Symbol',
    'Receptor_genes':    'Receptor_Gene_Symbol',
    'R':                 'Receptor_Gene_Symbol',

    'Receptor uniprot id':  'Receptor_Uniprot',
    'Receptor uniprot_ id': 'Receptor_Uniprot',
    'Receptor_Uniprot':     'Receptor_Uniprot',

    'ligand_symbol': 'Ligand_Gene_Symbol',
    'Ligand':        'Ligand_Gene_Symbol',
    'Ligand_geneid': 'Ligand_Gene_ID',
    'Ligand_Uniprot':'Ligand_Uniprot',

    'gene':      'Gene_Name',
    'Gene':      'Gene_Name',
    'Gene_name': 'Gene_Name',

    'protein name': 'Protein_Name',
    'Protein name': 'Protein_Name',
    'Protein_name': 'Protein_Name',

    'uniprot': 'Uniprot',
}

# ── Populated later by audit_and_standardize_headers ───────────────────────
case_rename_map = {}

# Columns used internally for pipeline tracking — never rename these
PROTECTED_COLS = {'database'}


def audit_and_standardize_headers(df_lists):
    """
    1. Collects every column header across all dataframes.
    2. Groups them by their lowercased/stripped form.
    3. Prints any groups that differ only in casing so you can review them.
    4. Builds and returns a case_rename_map  {bad_variant -> canonical}
       that is applied on top of rename_map inside standardize_cols.

    Canonical priority:
      a) already a value in rename_map  (e.g. Metabolite_Name)
      b) mixed-case and starts uppercase (e.g. Gene_Name, Class)
      c) first alphabetical variant as last resort

    Protected columns (e.g. 'database') are never renamed.
    """
    all_cols = set()
    for df_list in df_lists:
        for df in df_list:
            all_cols.update(df.columns.tolist())

    groups = defaultdict(list)
    for col in sorted(all_cols):
        groups[col.strip().lower()].append(col)

    result_map = {}
    print("\n=== Case-insensitive column conflicts ===")
    found = False
    for norm, variants in sorted(groups.items()):
        # Never touch protected internal columns
        if norm in PROTECTED_COLS:
            continue
        if len(variants) > 1:
            found = True
            canonical = (
                next((v for v in variants if v in rename_map.values()), None)
                or next((v for v in variants if v != v.upper() and v != v.lower() and v[0].isupper()), None)
                or sorted(variants)[0]
            )
            print(f"  {variants}  →  '{canonical}'")
            for v in variants:
                if v != canonical:
                    result_map[v] = canonical

    if not found:
        print("  No conflicts found.")
    print(f"\n{len(result_map)} columns will be renamed for case consistency.\n")
    return result_map


def standardize_cols(df):
    df = remove_index_cols(df)

    # Step 1: apply the explicit rename_map for genuinely different column names
    df = df.rename(columns=rename_map)

    # Step 2: apply the case-consistency map built by audit_and_standardize_headers
    df = df.rename(columns=case_rename_map)

    # Step 3: coalesce any remaining duplicates (safety net) — skip protected columns
    groups = defaultdict(list)
    for col in df.columns:
        if col.strip().lower() in PROTECTED_COLS:
            groups[col].append(col)  # keep as-is under its own key
        else:
            groups[col.strip().lower()].append(col)

    new_df = pd.DataFrame(index=df.index)
    for norm_name, variants in groups.items():
        if len(variants) == 1:
            new_df[variants[0]] = df[variants[0]]
        else:
            canonical = (
                next((v for v in variants if v in rename_map.values()), None)
                or next((v for v in variants if v != v.upper() and v != v.lower() and v[0].isupper()), None)
                or variants[0]
            )
            combined = df[variants].bfill(axis=1).iloc[:, 0]
            new_df[canonical] = combined

    return new_df


def enrich_pmids(df):
    pmid_list = []
    cols_to_check = ['PMID', 'source', 'comments_complex', 'comments', 'Database', 'Other.DB', 'Evidence', 'Annotation']
    avail_cols = [c for c in cols_to_check if c in df.columns]
    
    if not avail_cols:
        if 'PMID' not in df.columns:
            df['PMID'] = np.nan
        return df

    # We iterate over rows explicitly because we want to merge pmids from multiple columns
    for idx, row in df.iterrows():
        row_pmids = set()
        for c in avail_cols:
            val = row[c]
            if pd.notna(val):
                v_str = str(val).strip()
                if not v_str: continue
                
                matches = re.findall(r'(?i)(?:pmid|pubmed)\s*:?\s*(\d+)', v_str)
                row_pmids.update(matches)
                
                if c == 'PMID':
                    if v_str.endswith('.0'): v_str = v_str[:-2]
                    raw_nums = re.findall(r'\b\d+\b', v_str)
                    for n in raw_nums:
                        if len(n) >= 5: # basic check to avoid catching small random numbers
                            row_pmids.add(n)
        
        if row_pmids:
            pmid_list.append(';'.join(sorted(list(row_pmids))))
        else:
            pmid_list.append(np.nan)
            
    df['PMID'] = pmid_list
    return df

def deduplicate_cellinker_mrclink(df):
    if df.empty: return df
    mask = df['database'].isin(['Cellinker2', 'MRCLinkDB'])
    if not mask.any(): return df
    
    df_other = df[~mask].copy()
    df_cm = df[mask].copy()
    
    # Fuzzy Deduplication (Round 6)
    core_keys = ['Metabolite_Name', 'HMDB_ID', 'Gene_Name', 'Receptor_Gene_Symbol', 'Ligand_Gene_Symbol', 'Uniprot', 'Sensor_Gene', 'Enzyme', 'Transporter']
    group_keys = [k for k in core_keys if k in df_cm.columns]
    
    if not group_keys:
        return df
        
    def combine_dbs(x):
        s = set(x.dropna())
        if 'Cellinker2' in s and 'MRCLinkDB' in s:
            return 'MRCLinkDB/Cellinker2'
        return '/'.join(sorted(s))

    agg_dict = {'database': combine_dbs}
    for c in df_cm.columns:
        if c not in group_keys and c != 'database':
            agg_dict[c] = lambda val: val.dropna().iloc[0] if not val.dropna().empty else np.nan

    df_cm_agg = df_cm.groupby(group_keys, dropna=False, as_index=False).agg(agg_dict)
    
    # We must enforce the columns match the exact layout they were before aggregating
    # so pd.concat doesn't misalign types or create duplicate headers. Groupby sometimes reorders them.
    df_cm_agg = df_cm_agg[df_cm.columns]
    
    return pd.concat([df_other, df_cm_agg], ignore_index=True)

def deduplicate_keep_most_info(df):
    """
    Removes partial-duplicate rows within each database.

    A row B is a partial duplicate of row A if every non-null key value in B
    also exists in A with the same value (B's key set is a subset of A's).
    The sparser row is dropped; any non-null values it had that are missing
    in A are first coalesced into A so no information is lost.

    Complexity: O(n * k) via an inverted index over (col, val) pairs,
    avoiding the O(n^2) pairwise loop.
    """
    if df.empty:
        return df

    key_cols = [c for c in [
        'Metabolite_Name', 'HMDB_ID', 'Receptor_Gene_Symbol',
        'Ligand_Gene_Symbol', 'Gene_Name', 'Uniprot', 'Sensor_Gene',
        'Enzyme', 'Transporter', 'database'
    ] if c in df.columns]

    if not key_cols:
        return df

    def _is_empty(v):
        try:
            return pd.isna(v) or str(v).strip() == ''
        except Exception:
            return False

    df = df.copy().reset_index(drop=True)

    # Sort most-complete rows first so earlier rows always dominate later ones
    df['_nc'] = df.isnull().sum(axis=1)
    df = df.sort_values('_nc').reset_index(drop=True)
    df = df.drop(columns=['_nc'])

    # Normalise key columns into a parallel frame (lowercase + stripped)
    def _norm(x):
        if _is_empty(x):
            return np.nan
        return str(x).strip().lower()

    norm_keys = df[key_cols].apply(lambda col: col.map(_norm))

    # Build key signature (frozenset of (col, val) pairs) for each row
    sigs = [
        frozenset(
            (c, norm_keys.at[i, c])
            for c in key_cols
            if pd.notna(norm_keys.at[i, c])
        )
        for i in range(len(df))
    ]

    # Inverted index: (col, val) -> set of *kept* row indices that have that pair
    pair_to_rows = defaultdict(set)

    records = df.to_dict('records')
    drop_set = set()

    for j, sig_j in enumerate(sigs):
        if not sig_j:
            # No key info at all: cannot be identified as a partial duplicate
            continue

        # Candidate kept rows = intersection of all sets for each pair in sig_j.
        # Any row in the intersection has every pair of sig_j, i.e. sig_j is a subset.
        pair_iter = iter(sig_j)
        candidates = pair_to_rows[next(pair_iter)].copy()
        for pair in pair_iter:
            candidates &= pair_to_rows[pair]
            if not candidates:
                break

        if candidates:
            # Absorb j into the first matching kept row; coalesce missing values
            absorber = next(iter(candidates))
            for col, val in records[j].items():
                if not _is_empty(val) and _is_empty(records[absorber].get(col)):
                    records[absorber][col] = val
            drop_set.add(j)
        else:
            # j is unique — keep it and register its pairs in the index
            for pair in sig_j:
                pair_to_rows[pair].add(j)

    result = pd.DataFrame([rec for i, rec in enumerate(records) if i not in drop_set])
    return result.reset_index(drop=True)


def save_per_database(df_list, mol_type, species_prefix):
    if not df_list:
        return None
    
    std_dfs = [standardize_cols(d) for d in df_list]
    
    # Base vs Annotation Relational Segregation
    db_groups = defaultdict(list)
    for d in std_dfs:
        if not d.empty and 'database' in d.columns:
            db_groups[d['database'].iloc[0]].append(d)
            
    merged_groups = []
    outer_merge_keys = ['Metabolite_Name', 'HMDB_ID', 'Gene_Name', 'Receptor_Gene_Symbol', 'Ligand_Gene_Symbol', 'Uniprot']
    # Interaction target keys indicating a 'Base' (Fact) table
    fact_indicators = ['Gene_Name', 'Receptor_Gene_Symbol', 'Ligand_Gene_Symbol', 'Uniprot', 'Sensor_Gene']
    
    for db_name, dfs in db_groups.items():
        if len(dfs) == 1:
            merged_groups.append(dfs[0])
            continue
            
        # Segregate Base (Fact) vs Annotation (Dimension) tables
        base_tables = []
        anno_tables = []
        for d in dfs:
            if any(hc in d.columns for hc in fact_indicators):
                base_tables.append(d)
            else:
                anno_tables.append(d)
                
        # 1. Vertically append all Base/Fact tables so interactions don't explode horizontally
        if base_tables:
            res = pd.concat(base_tables, ignore_index=True)
        else:
            # If no base tables exist at all, just fallback to concatenating everything
            res = pd.concat(anno_tables, ignore_index=True)
            anno_tables = []
            
        # 2. Extract uniquely deduplicated Annotation/Dimension tables and Left-Join them
        for nxt in anno_tables:
            common_cols = list(set(res.columns) & set(nxt.columns))
            merge_on = [k for k in outer_merge_keys if k in common_cols]
            
            if not merge_on:
                res = pd.concat([res, nxt], ignore_index=True)
            else:
                # Deduplicate the dimension table on the join keys so it doesn't accidentally Cartesian product
                nxt_dedup = nxt.drop_duplicates(subset=merge_on, ignore_index=True)
                
                # Left join to broadcast properties without adding orphan parameter rows!
                res = pd.merge(res, nxt_dedup, on=merge_on, how='left')
                
                # Coalesce identical `_x` and `_y` duplicate columns created by pandas automatically
                for c in list(res.columns):
                    if str(c).endswith('_x'):
                        orig = str(c)[:-2]
                        y_col = orig + '_y'
                        if y_col in res.columns:
                            res[orig] = res[c].fillna(res[y_col])
                            res = res.drop(columns=[c, y_col])
                            
        res['database'] = db_name
        merged_groups.append(res)
        
    full_df = pd.concat(merged_groups, ignore_index=True)
    full_df = enrich_pmids(full_df)
    full_df = deduplicate_cellinker_mrclink(full_df)
    
    # Drop entirely empty columns to keep schema clean
    full_df = full_df.dropna(axis=1, how='all')

    # keep columns with >1 unique value OR named 'database'
    full_df = full_df.loc[:, (full_df.nunique() > 1) | (full_df.columns == 'database')]
    
    # Drop completely duplicated identical rows
    full_df = full_df.drop_duplicates(ignore_index=True)

    # Drop partial duplicates: rows that are a subset of another row
    # (same key values, but fewer columns filled in). Keep the richer row,
    # coalescing any unique non-null values from the sparser one into it.
    full_df = deduplicate_keep_most_info(full_df)
    if 'Metabolite_Name' in full_df.columns:
        full_df['Metabolite_Name'] = full_df['Metabolite_Name'].str.replace(r"(?i)SEROTONIN\s+DOPAMIN", "SEROTONIN;DOPAMINE", regex=True)
        full_df['Metabolite_Name'] = full_df['Metabolite_Name'].str.replace(r"(?i)ALDOSTERONE\s+CORTICOSTERONE", "ALDOSTERONE;CORTICOSTERONE", regex=True)
    # to make separation of metabolite by ";" consistently 
        # for metaligand, the "," has to be ";"
        mask = full_df['database'] == "MetaLigand"
        
        full_df.loc[mask, 'Metabolite_Name'] = (
            full_df.loc[mask, 'Metabolite_Name']
            .str.replace(r'(?<=[a-z\)]),\s*', ';', regex=True)
        )
        
        
        full_df['Metabolite_Name'] = full_df['Metabolite_Name'].str.replace(r"(\salpha\b|alpha\b|-a\b)", "a", regex=True)
        full_df['Metabolite_Name'] = full_df['Metabolite_Name'].str.replace(r"(\sbeta\b|beta\b|-b\b)", "b", regex=True)
        full_df['Metabolite_Name'] = full_df['Metabolite_Name'].str.replace(r"(?i)(5-SHp ETE|5SHp ETE|5SHpETE|5s-hpete|5S-HPETE)\s*", "5(S)-hydroperoxyeicosatetraenoic acid", regex=True)
        metab_dict = pd.read_csv("metab_dict.csv")
        lookup = {k.lower(): v.lower() for k, v in zip(metab_dict["metabolite_acronym"], metab_dict["Metabolite_Name"])}

    # --- Strict exact-match pass (metab_dict_strict.csv) ---
    # These entries require a full case-insensitive exact match on the entire
    # Metabolite_Name value (or sub-value if separated by ';') — no partial/substring substitution.
    if os.path.exists("metab_dict_strict.csv"):
        metab_dict_strict = pd.read_csv("metab_dict_strict.csv")
        strict_lookup = {
            k.lower().strip(): v.lower().strip()
            for k, v in zip(metab_dict_strict["metabolite_acronym"], metab_dict_strict["Metabolite_Name"])
        }
        if 'Metabolite_Name' in full_df.columns:
            def apply_strict(name):
                if not isinstance(name, str):
                    return name
                parts = str(name).split(';')
                mapped = []
                for p in parts:
                    clean_p = p.strip()
                    mapped.append(strict_lookup.get(clean_p.lower(), clean_p))
                return ';'.join(mapped)
            full_df['Metabolite_Name'] = full_df['Metabolite_Name'].apply(apply_strict)

    if 'Metabolite_Name' in full_df.columns:
        def apply_metab_lookup(name):
            if not isinstance(name, str):
                return name
            # Apply longest-key-first to avoid partial matches on shorter substrings.
            for key in sorted(lookup.keys(), key=len, reverse=True):
                canonical = lookup[key]
                # For entries that add an 'l-' prefix (bare amino acid → l-amino acid),
                # use a stricter lookbehind: (?<!l-)(?<!\w) — this prevents re-matching
                # 'histidine' inside an already-correct 'l-histidine' (l-l- problem),
                # but still blocks digit-hyphen prefixes too (e.g. '5-histidine' is fine).
                # For all other entries (acronyms, acid→anion, etc.) only (?<!\w) is used,
                # so '5-hete', '12-hpete', etc. still match after their digit-hyphen prefix.
                if re.match(r'^[a-zA-Z]-', canonical) and not re.match(r'^[a-zA-Z]-', key):
                    # Block matching when the key appears after ANY letter-hyphen prefix
                    # (d-, l-, n-, etc.) to prevent e.g. d-alanine → d-l-alanine.
                    pat = r'(?<![a-zA-Z]-)(?<!\w)' + re.escape(key) + r'(?!\w)'
                else:
                    pat = r'(?<!\w)' + re.escape(key) + r'(?!\w)'
                replaced = re.sub(pat, canonical, name, flags=re.IGNORECASE)
                if replaced != name:
                    name = replaced
            return name

        full_df['Metabolite_Name'] = full_df['Metabolite_Name'].apply(apply_metab_lookup)

    # ── Column normalisation ─────────────────────────────────────────────────
    # 1. Coalesce the three Super Class columns (from MEBOCOST, MRCLinkDB/Cellinker2,
    #    MetaLigand) into a single 'Super_Class'. ClassyFire vocabulary
    #    (super_class / Super Class) is preferred; MetaLigand's coarser 'SuperClass'
    #    fills gaps.
    _sc_cols = [c for c in ['super_class', 'Super Class', 'SuperClass'] if c in full_df.columns]
    if _sc_cols:
        full_df['Super_Class'] = full_df[_sc_cols[0]].copy()
        for _c in _sc_cols[1:]:
            full_df['Super_Class'] = full_df['Super_Class'].fillna(full_df[_c])
        full_df = full_df.drop(columns=_sc_cols)

    # 2. Coalesce synonym columns (MEBOCOST → synonyms_name, MetaLigand → Synonyms)
    if 'synonyms_name' in full_df.columns and 'Synonyms' in full_df.columns:
        full_df['Synonyms'] = full_df['Synonyms'].fillna(full_df['synonyms_name'])
        full_df = full_df.drop(columns=['synonyms_name'])
    elif 'synonyms_name' in full_df.columns:
        full_df = full_df.rename(columns={'synonyms_name': 'Synonyms'})

    # 3. Coalesce Evidence (MEBOCOST PMIDs) into the unified PMID column
    if 'Evidence' in full_df.columns and 'PMID' in full_df.columns:
        full_df['PMID'] = full_df['PMID'].fillna(full_df['Evidence'])
        full_df = full_df.drop(columns=['Evidence'])
    elif 'Evidence' in full_df.columns:
        full_df = full_df.rename(columns={'Evidence': 'PMID'})

    # 4. Simple renames for clarity
    _col_renames = {
        'ENZYME_NAME':          'Enzyme_Full_Name',           # full protein name vs gene symbol (Enzyme col)
        'HAMDBP_ID':            'HMDB_Protein_ID',            # HMDBP namespace (protein), not metabolite
        'UNIPROT_ID':           'Protein_Uniprot',            # enzyme/transporter Uniprot (distinct from Receptor_Uniprot)
        'Gene_Name':            'Target_Gene',                # scCellFie target gene, distinct from receptor columns
        'Transporter_genes':    'Metaligand_transporter_genes_out',      # MetaLigand GEM-derived transporters
        'Transporter_genes_in': 'Metaligand_transporter_genes_in',  # MetaLigand inward transporters (data-only, not used in R pipeline)
    }
    full_df = full_df.rename(columns={k: v for k, v in _col_renames.items() if k in full_df.columns})
    full_df = normalize_empty(full_df)
    # 5. Drop Sensor_Gene/Receptor — receptor rows already propagated to Receptor_Gene_Symbol
    #    and transporter/enzyme rows to Transporter / Enzyme by assign_sensor_cols().
    #    Sensor_Type is kept for context.
    if 'Sensor_Gene' in full_df.columns:
        full_df = full_df.drop(columns=['Sensor_Gene'])
    
    if 'receptor' in full_df.columns:
        full_df = full_df.drop(columns=['receptor'])
    
    saved_groups = []  # ← add this before the loop
    
    for db_name, group in full_df.groupby('database'):
        folder_to_save = db_name.split('/')[0] if '/' in db_name else db_name

        # ── MetaLigand: expand multi-HMDB_ID rows into compound Metabolite_Name ──
        if db_name == 'METALIGAND' and 'HMDB_ID' in group.columns and 'Metabolite_Name' in group.columns:
            group = group.copy()

            def update_name(row):
                hmdb_val = str(row['HMDB_ID']).strip() if pd.notna(row['HMDB_ID']) else ''
                metab_val = str(row['Metabolite_Name']).strip() if pd.notna(row['Metabolite_Name']) else ''

                if not hmdb_val or not metab_val:
                    return row['Metabolite_Name']

                if metab_val.upper() in ('NAN', '<NA>', 'NONE', 'NULL', 'N/A', 'NA'):
                    return row['Metabolite_Name']

                ids = [x.strip() for x in hmdb_val.split(',')]

                if len(ids) > 1 and ';' not in metab_val:
                    second_id = ids[1]
                    second_name = hmdb_dict.get(second_id)
                    if second_name:
                        return f"{metab_val};{second_name.upper()}"

                return metab_val

            group['Metabolite_Name'] = group.apply(update_name, axis=1)

        saved_groups.append(group)  # ← collect every group (modified or not)

        out_filename = f"{species_prefix}_{db_name.replace('/', '_')}_{mol_type}.csv"
        db_out_dir = os.path.join(out_dir, folder_to_save)
        os.makedirs(db_out_dir, exist_ok=True)
        out_path = os.path.join(db_out_dir, out_filename)

        group.dropna(axis=1, how='all').to_csv(out_path, index=False)
        print(f"Saved {out_path}")

    # Rebuild full_df from the (potentially modified) groups so the returned
    # dataframe reflects the MetaLigand update_name fix
    full_df = pd.concat(saved_groups, ignore_index=True)  # ← add this after the loop

    return full_df


# ── Audit headers across all collected dataframes BEFORE standardizing ──────
print("Auditing column headers for case-insensitive conflicts...")
all_lists = [all_human_metabolites, all_human_proteins, all_mouse_metabolites, all_mouse_proteins]
case_rename_map = audit_and_standardize_headers(all_lists)

print("Standardizing, Concatenating, Extracting PMIDs, Deduplicating, and Saving Dataframes...")

hm_df = save_per_database(all_human_metabolites, 'metabolites', 'human')
if hm_df is not None:
    hm_df.to_csv(os.path.join(out_dir, 'merged_human_metabolites.csv'), index=False)
    print(f"Saved global merged_human_metabolites.csv to {out_dir}")

hp_df = save_per_database(all_human_proteins, 'proteins', 'human')
if hp_df is not None:
    hp_df.to_csv(os.path.join(out_dir, 'merged_human_proteins.csv'), index=False)
    print(f"Saved global merged_human_proteins.csv to {out_dir}")

mm_df = save_per_database(all_mouse_metabolites, 'metabolites', 'mouse')
if mm_df is not None:
    mm_df.to_csv(os.path.join(out_dir, 'merged_mouse_metabolites.csv'), index=False)
    print(f"Saved global merged_mouse_metabolites.csv to {out_dir}")

mp_df = save_per_database(all_mouse_proteins, 'proteins', 'mouse')
if mp_df is not None:
    mp_df.to_csv(os.path.join(out_dir, 'merged_mouse_proteins.csv'), index=False)
    print(f"Saved global merged_mouse_proteins.csv to {out_dir}")

print("Merge processing complete.")