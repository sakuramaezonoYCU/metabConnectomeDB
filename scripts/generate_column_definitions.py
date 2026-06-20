import pandas as pd
import os
import json
import glob

def build_provenance_map(project_root):
    """
    Dynamically scan the raw database output subdirectories to trace exactly 
    which database(s) each column originates from.
    """
    provenance = {}
    output_dir = os.path.join(project_root, 'output')
    
    # Iterate through immediate subdirectories inside output/ (e.g., MEBOCOST, CellPhoneDBv5)
    if os.path.exists(output_dir):
        for item in os.listdir(output_dir):
            item_path = os.path.join(output_dir, item)
            if os.path.isdir(item_path):
                db_name = item
                # Scan all CSVs in this database folder
                for csv_file in glob.glob(os.path.join(item_path, '*.csv')):
                    try:
                        df = pd.read_csv(csv_file, nrows=0)
                        for col in df.columns:
                            if col not in provenance:
                                provenance[col] = set()
                            provenance[col].add(db_name)
                    except Exception:
                        pass
    return provenance

def generate_sub_definitions(cols, definitions_map, provenance_map, output_csv, output_md, title, desc):
    """
    Modular utility to map columns to their definitions, categories, and provenance,
    writing out clean CSV and formatted Markdown summaries.
    """
    data = []
    if not cols:
        print(f"Skipping {output_csv} as no columns were found (files probably don't exist yet).")
        return
        
    for col in sorted(cols):
        # Determine provenance dynamically from the source files
        if col in provenance_map and provenance_map[col]:
            db = ", ".join(sorted(list(provenance_map[col])))
        else:
            db = 'Derived/Computed'

        if col in definitions_map:
            category, definition = definitions_map[col]
        else:
            category = 'Interaction Evidence & Scores'
            definition = 'Consolidated metadata field.'
            
        data.append({
            'header': col,
            'database': db,
            'category': category,
            'definition': definition
        })
        
    df_out = pd.DataFrame(data)
    df_out.to_csv(output_csv, index=False)
    print(f"🎉 Successfully created column definitions at {output_csv}")
    
    # Save beautiful Markdown grouped by Category
    with open(output_md, 'w', encoding='utf-8') as f:
        f.write(f"# {title}\n\n")
        f.write(f"{desc}\n\n")
        
        categories = df_out['category'].unique()
        for cat in sorted(categories):
            f.write(f"## 📁 {cat}\n\n")
            f.write("| Header | Database | Definition |\n")
            f.write("| --- | --- | --- |\n")
            df_cat = df_out[df_out['category'] == cat]
            for _, r in df_cat.iterrows():
                f.write(f"| `{r['header']}` | **{r['database']}** | {r['definition']} |\n")
            f.write("\n")
            
    print(f"🎉 Successfully created markdown documentation at {output_md}")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    # Load configuration from centralized JSON file
    config_path = os.path.join(script_dir, 'column_definitions.json')
    with open(config_path, 'r', encoding='utf-8') as f:
        definitions_map = json.load(f)

    # Core Database Inputs
    pair_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
    metab_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_with_HMDB_Info.csv')
    
    # Unique Metabolite Exploration Notebook Output Inputs
    tier1_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_tier1.csv')
    tier2_csv = os.path.join(project_root, 'output', 'human_database_merge_unique_metab_tier2.csv')
    
    # Metabolite-Target Pair Analysis Notebook Output Inputs
    cancer_csv = os.path.join(project_root, 'output', 'human_metab_target_pairs_cancer_annotated.csv')

    # Main Database Outputs
    main_csv = os.path.join(project_root, 'output', 'merged_database_col_definition.csv')
    main_md = os.path.join(project_root, 'output', 'merged_database_col_definition.md')
    
    # Notebook Output Definitions Outputs
    unique_csv = os.path.join(project_root, 'output', 'unique_metab_exploration_col_definition.csv')
    unique_md = os.path.join(project_root, 'output', 'unique_metab_exploration_col_definition.md')
    
    pair_analysis_csv = os.path.join(project_root, 'output', 'metab_target_pair_analysis_col_definition.csv')
    pair_analysis_md = os.path.join(project_root, 'output', 'metab_target_pair_analysis_col_definition.md')

    # Build dynamic provenance map by scanning the raw output databases
    provenance_map = build_provenance_map(project_root)

    # 1. Main merged database column definitions
    main_cols = set()
    if os.path.exists(pair_csv):
        df_pair = pd.read_csv(pair_csv, nrows=1)
        main_cols.update(df_pair.columns)
    if os.path.exists(metab_csv):
        df_metab = pd.read_csv(metab_csv, nrows=1)
        main_cols.update(df_metab.columns)
        
    generate_sub_definitions(
        main_cols, definitions_map, provenance_map, main_csv, main_md,
        "MetabConnectomeDB Merged Database Column Definitions",
        "This file lists all unified database column headers in the consolidated database, their functional categories, their original source database (or `multiple` if shared/common), and their exact definitions."
    )
    
    # 2. Unique metabolite exploration notebook outputs definitions
    unique_cols = set()
    if os.path.exists(tier1_csv):
        df_t1 = pd.read_csv(tier1_csv, nrows=1)
        unique_cols.update(df_t1.columns)
    if os.path.exists(tier2_csv):
        df_t2 = pd.read_csv(tier2_csv, nrows=1)
        unique_cols.update(df_t2.columns)
        
    generate_sub_definitions(
        unique_cols, definitions_map, provenance_map, unique_csv, unique_md,
        "Unique Metabolite Exploration Column Definitions (Notebook Outputs)",
        "This file lists the column headers, categories, and definitions for the outputs of the unique metabolite exploration Jupyter notebook (`unique_metab_data_exploration.ipynb`), specifically `human_database_merge_unique_metab_tier1.csv` and `human_database_merge_unique_metab_tier2.csv`."
    )
    
    # 3. Metabolite-Target pair analysis notebook outputs definitions
    pair_analysis_cols = set()
    if os.path.exists(cancer_csv):
        df_cancer = pd.read_csv(cancer_csv, nrows=1)
        pair_analysis_cols.update(df_cancer.columns)
        
    generate_sub_definitions(
        pair_analysis_cols, definitions_map, provenance_map, pair_analysis_csv, pair_analysis_md,
        "Metabolite Target Pair Analysis Column Definitions (Notebook Outputs)",
        "This file lists the column headers, categories, and definitions for the outputs of the metabolite-target pair analysis Jupyter notebook (`metab_targetPair_analysis.ipynb`), specifically the clinical and cancer-annotated file `human_metab_target_pairs_cancer_annotated.csv`."
    )

if __name__ == '__main__':
    main()
