import nbformat as nbf
import os
import argparse
import sys
import subprocess

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_DIR = os.path.join(BASE_DIR, 'output', 'camp_integration')

def create_notebook(signature_name, genes):
    nb = nbf.v4.new_notebook()
    
    # 1. Title and Goal Markdown
    nb.cells.append(nbf.v4.new_markdown_cell(f"""# Pan-Cancer CAMP Metabolomics Integration

### Goal
Integrate the {signature_name} metabolic signature with the pan-cancer mass spectrometry metabolomics dataset (CAMP, PMID: 37337120).

### Purpose
To dynamically analyze the abundance of metabolites associated with a specific gene signature across multiple cancer cohorts. It directly maps genes to metabolites using the human database, extracts matching metabolomic/transcriptomic data, and evaluates differential abundance, gene-metabolite covariation, and immune microenvironment correlation.

### Interpretation
- **Tumor vs Normal Differential Abundance**: Identifies if the signature metabolites are consistently enriched or depleted in tumor tissue compared to normal adjacent tissue.
- **Gene-Metabolite Covariation**: A strong Spearman correlation between a gene's expression and its associated metabolite's abundance in tumors provides direct multi-omic validation of the metabolic axis.
- **Immune Microenvironment (TME) Correlation**: Significant associations between metabolite levels and specific immune cell populations (e.g., Macrophages, Mast cells) highlight potential mechanisms of metabolic immune evasion or interaction.

### Inputs
- **CAMP Data Directory:** `input/pancancer_metabolomics_2023_PMID37337120/data/` (contains `metabolomics_processed`, `transcriptomics_processed`, `TME_deconvolution_processed`)
- **Metabolite Database:** `input/databases/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`
- **Gene Signature Name:** {signature_name}
- **Genes:** {', '.join(genes)}

### Outputs
- **Plots & Output CSV:** `output/camp_integration/`

### Example Usage
To run the analysis for this {signature_name} signature:
1. Ensure the parameters cell below is configured properly (it is already pre-filled with the `{signature_name}` parameters).
2. Click **Run All** in Jupyter.
3. The integrated plots will be displayed inline and the master dataset will be saved as `master_integrated_camp_data.csv`.
"""))

    # 2. Config Code
    nb.cells.append(nbf.v4.new_code_cell(f"""import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ==========================================
# PARAMETERS
# ==========================================
# Gene signature name and the genes themselves
SIGNATURE_NAME = "{signature_name}"
TARGET_GENES = {genes}

# Which cohorts to load. 
# Available TCGA-like in CAMP: 'BRCA1', 'BRCA2', 'COAD', 'DLBCL', 'GBM', 'HCC', 'HurthleCC', 'ICC', 'OV', 'PDAC', 'PRAD', 'ccRCC1', 'ccRCC2', 'ccRCC3', 'ccRCC4'
# Leave as None to load all available inside the directory.
TARGET_COHORTS = ['BRCA1', 'BRCA2', 'COAD', 'OV', 'PRAD', 'GBM', 'PDAC']

# Paths
try:
    BASE_DIR = os.path.abspath(os.path.dirname(__file__))
    if os.path.basename(BASE_DIR) == 'scripts':
        BASE_DIR = os.path.dirname(BASE_DIR)
except NameError:
    BASE_DIR = os.path.abspath('.')
    if os.path.basename(BASE_DIR) == 'scripts':
        BASE_DIR = os.path.dirname(BASE_DIR)

INPUT_DIR = os.path.join(BASE_DIR, 'input')
CAMP_DIR = os.path.join(INPUT_DIR, 'pancancer_metabolomics_2023_PMID37337120', 'data')
DB_PATH = os.path.join(INPUT_DIR, 'databases', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
OUTPUT_DIR = os.path.join(BASE_DIR, 'output', 'camp_integration')
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Signature: {{SIGNATURE_NAME}}")
print(f"Target Genes: {{TARGET_GENES}}")
"""))

    nb.cells.append(nbf.v4.new_code_cell("""# ==========================================
# 1. LOAD DATABASE AND EXTRACT TARGET METABOLITES
# ==========================================
print("Loading database...")
db = pd.read_csv(DB_PATH, low_memory=False)

# Filter for the target genes
sig_db = db[db['Target'].isin(TARGET_GENES)].copy()
signature_metabolites = sig_db['Metabolite_Name'].dropna().unique().tolist()

print(f"Found {len(signature_metabolites)} metabolites associated with the target genes:")
for g in TARGET_GENES:
    assoc_mets = sig_db[sig_db['Target'] == g]['Metabolite_Name'].dropna().unique()
    print(f"  - {g}: {len(assoc_mets)} metabolites -> {', '.join(assoc_mets[:5])}...")
"""))

    nb.cells.append(nbf.v4.new_code_cell("""# ==========================================
# 2. DATA LOADING LOOP (Transcriptomics, Metabolomics, TME)
# ==========================================

metab_dir = os.path.join(CAMP_DIR, 'metabolomics_processed')
trans_dir = os.path.join(CAMP_DIR, 'transcriptomics_processed')
tme_dir = os.path.join(CAMP_DIR, 'TME_deconvolution_processed')
master_mapping_file = os.path.join(CAMP_DIR, 'MasterMapping_MetImmune_03_16_2022_release.csv')

# Load the master mapping to find the correct file names for each cohort
mapping_df = pd.read_csv(master_mapping_file, low_memory=False)
cohort_to_rna = dict(zip(mapping_df['Dataset'], mapping_df['RNAFile']))
cohort_to_tme = dict(zip(mapping_df['Dataset'], mapping_df['ITHFile']))
cohort_to_metab = dict(zip(mapping_df['Dataset'], mapping_df['MetabFile']))

all_cohorts_data = []

if TARGET_COHORTS is None:
    TARGET_COHORTS = list(cohort_to_metab.keys())

print(f"Processing Cohorts: {TARGET_COHORTS}")

for cohort in TARGET_COHORTS:
    cohort_mapping = mapping_df[mapping_df['Dataset'] == cohort]
    metab_to_common = dict(zip(cohort_mapping['MetabID'], cohort_mapping['CommonID']))
    rna_to_common = dict(zip(cohort_mapping['RNAID'], cohort_mapping['CommonID']))
    ith_to_common = dict(zip(cohort_mapping['ITHID'], cohort_mapping['CommonID']))
    common_to_tn = dict(zip(cohort_mapping['CommonID'], cohort_mapping['TN']))

    # Load Metabolomics
    metab_filename = cohort_to_metab.get(cohort, f"PreprocessedData_{cohort}.xlsx")
    metab_file = os.path.join(metab_dir, metab_filename)
    if not os.path.exists(metab_file):
        print(f"Skipping {cohort} - no metabolomics file found.")
        continue
    
    df_metab = pd.read_excel(metab_file, index_col=0)
    if df_metab.shape[0] >= df_metab.shape[1]:
        df_metab = df_metab.T

    keep_metabs = [m for m in df_metab.columns if m in signature_metabolites]
    df_metab = df_metab[keep_metabs]
    df_metab.columns = [f"METAB_{c}" for c in df_metab.columns]
    df_metab.index = df_metab.index.map(lambda x: metab_to_common.get(x, x))
    
    # Load Transcriptomics
    trans_filename = cohort_to_rna.get(cohort, None)
    if trans_filename and os.path.exists(os.path.join(trans_dir, trans_filename)):
        df_trans = pd.read_csv(os.path.join(trans_dir, trans_filename), index_col=0)
        if df_trans.shape[0] >= df_trans.shape[1]:
            df_trans = df_trans.T
        keep_genes = [g for g in df_trans.columns if g in TARGET_GENES]
        df_trans = df_trans[keep_genes]
        df_trans.columns = [f"GENE_{c}" for c in df_trans.columns]
        df_trans.index = df_trans.index.map(lambda x: rna_to_common.get(x, x))
    else:
        df_trans = pd.DataFrame(index=df_metab.index)
        print(f"  Warning: No TPM transcriptomics found for {cohort}")

    # Load TME Deconvolution
    tme_filename = cohort_to_tme.get(cohort, None)
    if tme_filename and os.path.exists(os.path.join(tme_dir, tme_filename)):
        df_tme = pd.read_csv(os.path.join(tme_dir, tme_filename), index_col=0)
        df_tme.columns = [f"TME_{c}" for c in df_tme.columns]
        df_tme.index = df_tme.index.map(lambda x: ith_to_common.get(x, x))
    else:
        df_tme = pd.DataFrame(index=df_metab.index)
        
    df_cohort = df_metab.join(df_trans, how='outer').join(df_tme, how='outer')
    df_cohort['Cohort'] = cohort
    
    # Determine Tissue_Type using the master mapping TN column.
    # If not found, check if it was originally marked as Normal, else Tumor
    def get_tissue_type(x):
        tn = common_to_tn.get(x)
        if pd.isna(tn): tn = None
        if tn:
            return str(tn).capitalize()
        return 'Normal' if str(x).endswith('N') or '_N' in str(x) else 'Tumor'

    df_cohort['Tissue_Type'] = df_cohort.index.map(get_tissue_type)
    all_cohorts_data.append(df_cohort)

master_df = pd.concat(all_cohorts_data)
print(f"\\nSuccessfully assembled multi-omics master table: {master_df.shape[0]} samples, {master_df.shape[1]} features.")
"""))

    nb.cells.append(nbf.v4.new_code_cell("""# ==========================================
# 4. TUMOR vs NORMAL ANALYSIS
# ==========================================
tumor_normal_df = master_df.dropna(subset=[c for c in master_df.columns if c.startswith('METAB_') or c.startswith('GENE_')], how='all')
features_to_plot = [c for c in master_df.columns if c.startswith('METAB_') or c.startswith('GENE_')]

if len(features_to_plot) > 0 and len(tumor_normal_df['Tissue_Type'].unique()) > 1:
    fig, axes = plt.subplots(len(features_to_plot), 1, figsize=(10, 4*len(features_to_plot)))
    if len(features_to_plot) == 1: axes = [axes]
    
    for i, feature in enumerate(features_to_plot):
        sns.boxplot(data=tumor_normal_df, x='Cohort', y=feature, hue='Tissue_Type', ax=axes[i])
        axes[i].set_title(f"Tumor vs Normal: {feature}")
        axes[i].tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'tumor_vs_normal_boxplots_{SIGNATURE_NAME}.png'))
    plt.show()
else:
    print("Not enough overlapping Tumor/Normal labels or features to plot.")
"""))

    nb.cells.append(nbf.v4.new_code_cell("""# ==========================================
# 4. GENE-METABOLITE COVARIATION (CONCORDANCE)
# ==========================================
tumor_df = master_df[master_df['Tissue_Type'] == 'Tumor']
gene_cols = [c for c in tumor_df.columns if c.startswith('GENE_')]
metab_cols = [c for c in tumor_df.columns if c.startswith('METAB_')]

if gene_cols and metab_cols:
    corr_matrix = pd.DataFrame(index=gene_cols, columns=metab_cols)
    for g in gene_cols:
        for m in metab_cols:
            valid = tumor_df[[g, m]].dropna()
            if len(valid) > 10:
                rho, p = stats.spearmanr(valid[g], valid[m])
                corr_matrix.loc[g, m] = rho

    corr_matrix = corr_matrix.astype(float)
    plt.figure(figsize=(10, max(4, len(gene_cols)*0.5)))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, vmin=-1, vmax=1)
    plt.title("Gene-Metabolite Spearman Correlation (Tumors only)")
    plt.savefig(os.path.join(OUTPUT_DIR, f'gene_metabolite_covariation_{SIGNATURE_NAME}.png'))
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""# ==========================================
# 5. METABOLITE - TME (IMMUNE CELL) ASSOCIATIONS
# ==========================================
tme_cols = [c for c in tumor_df.columns if c.startswith('TME_')]
immune_cells = [c for c in tme_cols if any(x in c.lower() for x in ['mast', 'macrophage', 'dc', 'dendritic', 't cell', 'b cell', 'nk'])]

if metab_cols and immune_cells:
    tme_corr = pd.DataFrame(index=metab_cols, columns=immune_cells)
    for m in metab_cols:
        for tme in immune_cells:
            valid = tumor_df[[m, tme]].dropna()
            if len(valid) > 10:
                rho, p = stats.spearmanr(valid[m], valid[tme])
                tme_corr.loc[m, tme] = rho

    tme_corr = tme_corr.astype(float)
    plt.figure(figsize=(12, max(5, len(metab_cols)*0.8)))
    sns.heatmap(tme_corr, annot=False, cmap='PRGn', center=0)
    plt.title("Metabolite vs Immune Cell Population Correlation")
    plt.savefig(os.path.join(OUTPUT_DIR, f'metabolite_immune_covariation_{SIGNATURE_NAME}.png'))
    plt.show()
"""))

    nb.cells.append(nbf.v4.new_code_cell("""# Save Master Data to Output
output_csv = os.path.join(OUTPUT_DIR, f'master_integrated_camp_data_{SIGNATURE_NAME}.csv')
master_df.to_csv(output_csv)
print(f"Analysis complete. Full merged data saved to {output_csv}")
"""))

    nb.cells.append(nbf.v4.new_code_cell(f"""# ==========================================
# 6. HTML EXPORT
# ==========================================
import subprocess
import sys
import os

notebook_filename = 'camp_pancancer_integration_{signature_name}.ipynb'
output_base = 'camp_pancancer_integration_{signature_name}'
output_dir = OUTPUT_DIR
os.makedirs(output_dir, exist_ok=True)

jupyter_bin = os.path.join(os.path.dirname(sys.executable), 'jupyter')
if not os.path.exists(jupyter_bin): jupyter_bin = 'jupyter'

cmd_html = [jupyter_bin, "nbconvert", "--to", "html", notebook_filename, "--output-dir", output_dir, "--output", output_base]
res_html = subprocess.run(cmd_html, capture_output=True, text=True)

if res_html.returncode == 0:
    print(f"🎉 SUCCESS: Notebook successfully exported to '{{os.path.join(output_dir, output_base)}}.html'")
else:
    print("❌ HTML export failed.")
    print(res_html.stderr)
"""))

    notebook_filename = os.path.join(BASE_DIR, 'scripts', f'camp_pancancer_integration_{signature_name}.ipynb')
    with open(notebook_filename, 'w') as f:
        nbf.write(nb, f)
    print(f"Notebook generated successfully at {notebook_filename}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate Mass Spec Notebook")
    parser.add_argument('--signature-name', default="Directed Metastatic Signature", help="Name of the signature")
    parser.add_argument('--genes', nargs='+', default="GLS SGMS1 SPTLC1 GBE1 SLC16A7 AUH FZD6 NR1D2 CD46 MTMR1 ESRRG ITGA4 SLC11A2 ERAP1 C1GALT1 ADAM10 TRPM8 SLC22A1 AMDHD1 EPOR PDE3B".split(), help="List of genes in the signature")
    args = parser.parse_args()
    
    create_notebook(args.signature_name, args.genes)
