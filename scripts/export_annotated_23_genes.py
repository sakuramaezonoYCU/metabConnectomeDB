import sys
if '..' not in sys.path: sys.path.append('..')
from pan_cancer_config import ANALYSIS_SUFFIX
import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(BASE_DIR), 'output')
META_RESULTS_DIR = os.path.join(OUTPUT_DIR, 'pan_cancer_meta_results')

CANCERS = ['breast', 'colorectal', 'lung', 'melanoma', 'ovarian']

# 1. Read 23 genes
genes_file = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv')
if not os.path.exists(genes_file):
    print(f"pan_cancer_conserved_genes{ANALYSIS_SUFFIX}.csv not found!")
    exit(1)
    
df_genes = pd.read_csv(genes_file)
gene_list = df_genes['Strictly_Conserved_Gene'].tolist()

# 2. Build the result rows
merged_data = {}
for g in gene_list:
    merged_data[g] = {'Gene': g, 'Metabolites': '', 'HMDB_ID': ''}

# 3. Read each cancer's DE file and extract LFC and Score
for cancer in CANCERS:
    res_file = os.path.join(OUTPUT_DIR, f"{cancer}_results", f"primary_vs_metastasis_{cancer}_results_DE_metabolic_targets.csv")
    if os.path.exists(res_file):
        df_cancer = pd.read_csv(res_file)
        # Filter for our 23 genes
        df_cancer = df_cancer[df_cancer['names'].isin(gene_list)]
        
        for _, row in df_cancer.iterrows():
            g = row['names']
            # Populate LFC and Score
            merged_data[g][f'{cancer.capitalize()}_LFC'] = row['logfoldchanges']
            merged_data[g][f'{cancer.capitalize()}_Score'] = row['scores']
            
            # They all should have the same Metabolites/HMDB_ID, we can grab it from any cancer
            if merged_data[g]['Metabolites'] == '':
                merged_data[g]['Metabolites'] = row['Metabolite_Name']
                merged_data[g]['HMDB_ID'] = row['HMDB_ID']

# 4. Convert to DataFrame and save
df_annotated = pd.DataFrame(list(merged_data.values()))

# Reorder columns to be nice
cols = ['Gene', 'Metabolites', 'HMDB_ID']
for cancer in CANCERS:
    c = cancer.capitalize()
    if f'{c}_LFC' in df_annotated.columns:
        cols.append(f'{c}_LFC')
        cols.append(f'{c}_Score')

df_annotated = df_annotated[cols]

out_file = os.path.join(META_RESULTS_DIR, f'pan_cancer_conserved_genes_with_annotation{ANALYSIS_SUFFIX}.csv')
df_annotated.to_csv(out_file, index=False)
print(f"Created annotated 23 genes file at {out_file}")
