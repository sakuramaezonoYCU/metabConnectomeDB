import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ==========================================
# PARAMETERS
# ==========================================
# Gene signature name and the genes themselves
SIGNATURE_NAME = "3-gene"
TARGET_GENES = ['MYC', 'E2F1', 'TP53']

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
print(f"Signature: {SIGNATURE_NAME}")
print(f"Target Genes: {TARGET_GENES}")

# ==========================================
# 1. LOAD DATABASE AND EXTRACT TARGET METABOLITES
# ==========================================
print("Loading database...")
db = pd.read_csv(DB_PATH, low_memory=False)

# Filter for the target genes
sig_db = db[db['Target_Gene'].isin(TARGET_GENES)].copy()
signature_metabolites = sig_db['Metabolite_Name'].dropna().unique().tolist()

print(f"Found {len(signature_metabolites)} metabolites associated with the target genes:")
for g in TARGET_GENES:
    assoc_mets = sig_db[sig_db['Target_Gene'] == g]['Metabolite_Name'].dropna().unique()
    print(f"  - {g}: {len(assoc_mets)} metabolites -> {', '.join(assoc_mets[:5])}...")

# ==========================================
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
    
    # Load Transcriptomics
    trans_filename = cohort_to_rna.get(cohort, None)
    if trans_filename and os.path.exists(os.path.join(trans_dir, trans_filename)):
        df_trans = pd.read_csv(os.path.join(trans_dir, trans_filename), index_col=0)
        if df_trans.shape[0] >= df_trans.shape[1]:
            df_trans = df_trans.T
        keep_genes = [g for g in df_trans.columns if g in TARGET_GENES]
        df_trans = df_trans[keep_genes]
        df_trans.columns = [f"GENE_{c}" for c in df_trans.columns]
    else:
        df_trans = pd.DataFrame(index=df_metab.index)
        print(f"  Warning: No TPM transcriptomics found for {cohort}")

    # Load TME Deconvolution
    tme_filename = cohort_to_tme.get(cohort, None)
    if tme_filename and os.path.exists(os.path.join(tme_dir, tme_filename)):
        df_tme = pd.read_csv(os.path.join(tme_dir, tme_filename), index_col=0)
        df_tme.columns = [f"TME_{c}" for c in df_tme.columns]
    else:
        df_tme = pd.DataFrame(index=df_metab.index)
        
    df_cohort = df_metab.join(df_trans, how='outer').join(df_tme, how='outer')
    df_cohort['Cohort'] = cohort
    df_cohort['Tissue_Type'] = df_cohort.index.map(lambda x: 'Normal' if str(x).endswith('N') or '_N' in str(x) else 'Tumor')
    all_cohorts_data.append(df_cohort)

master_df = pd.concat(all_cohorts_data)
print(f"\nSuccessfully assembled multi-omics master table: {master_df.shape[0]} samples, {master_df.shape[1]} features.")

# ==========================================
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

# ==========================================
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

# ==========================================
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

# Save Master Data to Output
output_csv = os.path.join(OUTPUT_DIR, f'master_integrated_camp_data_{SIGNATURE_NAME}.csv')
master_df.to_csv(output_csv)
print(f"Analysis complete. Full merged data saved to {output_csv}")
