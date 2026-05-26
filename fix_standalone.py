import os
import nbformat as nbf

notebook_path1 = 'scripts/primary_vs_metastasis_comparison.ipynb'
with open(notebook_path1, 'r', encoding='utf-8') as f:
    nb1 = nbf.read(f, as_version=4)

for cell in nb1.cells:
    if cell.cell_type == 'code':
        source = cell.source
        
        # 1. Add GLOBAL_OUTPUT_DIR
        if "OUTPUT_DIR = globals().get('OUTPUT_DIR'" in source:
            source = source.replace(
                "OUTPUT_DIR = globals().get('OUTPUT_DIR', os.path.join(BASE_DIR, 'output'))",
                "GLOBAL_OUTPUT_DIR = os.path.join(BASE_DIR, 'output')\nOUTPUT_DIR = globals().get('OUTPUT_DIR', os.path.join(GLOBAL_OUTPUT_DIR, 'breast_cancer'))"
            )
            
        # 2. Fix metab_db_path to use GLOBAL_OUTPUT_DIR
        if "metab_db_path = os.path.join(OUTPUT_DIR, 'human_database_merge" in source:
            source = source.replace(
                "metab_db_path = os.path.join(OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')",
                "metab_db_path = os.path.join(GLOBAL_OUTPUT_DIR, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')"
            )
            
        # 3. Update h5ad_path default
        if "h5ad_path = os.path.join(OUTPUT_DIR, 'cancer_breast-cancer" in source:
            source = source.replace(
                "h5ad_path = os.path.join(OUTPUT_DIR, 'cancer_breast-cancer_breast_100k_whole_transcriptome_2025-11-08.h5ad')",
                "h5ad_path = globals().get('h5ad_path', os.path.join(OUTPUT_DIR, 'breast_mammary-gland_liver_axilla_chest-wall_100k_whole_transcriptome_2025-11-08.h5ad'))"
            )
            
        cell.source = source

with open(notebook_path1, 'w', encoding='utf-8') as f:
    nbf.write(nb1, f)
print("Updated primary_vs_metastasis_comparison.ipynb")


notebook_path2 = 'scripts/orphan_metabolic_immune_evasion.ipynb'
with open(notebook_path2, 'r', encoding='utf-8') as f:
    nb2 = nbf.read(f, as_version=4)

for cell in nb2.cells:
    if cell.cell_type == 'code':
        source = cell.source
        
        # 1. Add GLOBAL_OUTPUT_DIR
        if "output_dir = globals().get('OUTPUT_DIR'" in source:
            source = source.replace(
                "output_dir = globals().get('OUTPUT_DIR', os.path.join(workspace_dir, 'output'))",
                "global_output_dir = os.path.join(workspace_dir, 'output')\noutput_dir = globals().get('OUTPUT_DIR', os.path.join(global_output_dir, 'breast_cancer'))"
            )
            
        # 2. Fix pairs_df to use global_output_dir
        if "pairs_df = pd.read_csv(os.path.join(output_dir, 'human_database_merge" in source:
            source = source.replace(
                "pairs_df = pd.read_csv(os.path.join(output_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'))",
                "pairs_df = pd.read_csv(os.path.join(global_output_dir, 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv'))"
            )
            
        cell.source = source

with open(notebook_path2, 'w', encoding='utf-8') as f:
    nbf.write(nb2, f)
print("Updated orphan_metabolic_immune_evasion.ipynb")
