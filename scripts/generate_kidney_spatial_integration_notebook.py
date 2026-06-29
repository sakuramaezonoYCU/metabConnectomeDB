import os
import sys
import nbformat as nbf

def generate_notebook():
    nb = nbf.v4.new_notebook()
    
    cells = []
    
    cells.append(nbf.v4.new_markdown_cell("""# Kidney Cancer (RCC): Mass Spec & Spatial Transcriptomics Integration

This notebook dynamically identifies top altered metabolites from bulk Mass Spec data for Kidney Renal Clear Cell Carcinoma (KIRC/ccRCC), maps them to their cognate enzymes/receptors using our integrated database, and then projects the expression of these genes onto RCC spatial transcriptomics coordinates (Kalogirou et al. 2025). This serves as a proof-of-concept for spatially resolving metabolic adaptations in kidney cancer."""))

    cells.append(nbf.v4.new_code_cell("""import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
import glob
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname('__file__')))
if os.path.basename(BASE_DIR) == 'scripts':
    BASE_DIR = os.path.dirname(BASE_DIR)

MASS_SPEC_PATH = os.path.join(BASE_DIR, 'input', 'massSpecDataMetabolicData_7panCancer_PMID29396322.csv')
DB_PATH = os.path.join(BASE_DIR, 'output', 'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
SPATIAL_DIR = os.path.join(BASE_DIR, 'input', 'spatial', 'Zenodo_16833780_RCC')
"""))

    cells.append(nbf.v4.new_markdown_cell("""## 1. Load and Process Mass Spec Data"""))
    cells.append(nbf.v4.new_code_cell("""# Load mass spec data
ms_df = pd.read_csv(MASS_SPEC_PATH)
ms_df = ms_df.rename(columns={'Unnamed: 0': 'Index', 'X': 'Metabolite'})
print(f"Loaded {ms_df.shape[0]} metabolites and {ms_df.shape[1]} columns.")

# Filter columns for Kidney (KIRC)
kirc_cols = [c for c in ms_df.columns if c.startswith('KIRC.')]
tumor_cols = [c for c in kirc_cols if 'Tumor' in c]
normal_cols = [c for c in kirc_cols if 'Normal' in c]

print(f"Found {len(tumor_cols)} Tumor samples and {len(normal_cols)} Normal samples for KIRC.")

# Calculate fold change
tumor_mean = ms_df[tumor_cols].mean(axis=1)
normal_mean = ms_df[normal_cols].mean(axis=1)

# Add small epsilon to avoid div by zero
epsilon = 1e-6
ms_df['Log2FC_Tumor_vs_Normal'] = np.log2((tumor_mean + epsilon) / (normal_mean + epsilon))

# Filter out rows with NaNs or inf
ms_df = ms_df.replace([np.inf, -np.inf], np.nan).dropna(subset=['Log2FC_Tumor_vs_Normal'])

# Select top altered metabolites (both up and down)
top_n = 10
top_up = ms_df.nlargest(top_n, 'Log2FC_Tumor_vs_Normal')[['Metabolite', 'Log2FC_Tumor_vs_Normal']]
top_down = ms_df.nsmallest(top_n, 'Log2FC_Tumor_vs_Normal')[['Metabolite', 'Log2FC_Tumor_vs_Normal']]

print("\\nTop Upregulated Metabolites:")
display(top_up)
print("\\nTop Downregulated Metabolites:")
display(top_down)

top_metabolites = pd.concat([top_up, top_down])['Metabolite'].tolist()

# Save for downstream scraping
out_dir = os.path.join(BASE_DIR, 'output')
os.makedirs(out_dir, exist_ok=True)
top_up.to_csv(os.path.join(out_dir, 'kidney_top_up_metabolites.csv'), index=False)
top_down.to_csv(os.path.join(out_dir, 'kidney_top_down_metabolites.csv'), index=False)
"""))

    cells.append(nbf.v4.new_markdown_cell("""## 2. Gene-Metabolite Mapping"""))
    cells.append(nbf.v4.new_code_cell("""# Load the integrated database
db_df = pd.read_csv(DB_PATH)

# Basic fuzzy matching since names might slightly differ
def match_metabolite(m_name, db_names):
    m_name_lower = str(m_name).lower().strip()
    # Direct match
    if m_name_lower in db_names:
        return [m_name_lower]
    # Substring match
    matches = [name for name in db_names if type(name) == str and (m_name_lower in name.lower() or name.lower() in m_name_lower)]
    return matches

db_metabolites = set(db_df['Metabolite_Name'].astype(str).str.lower().dropna())
matched_genes = set()

print("Mapping top metabolites to genes...")
for m in top_metabolites:
    matches = match_metabolite(m, db_metabolites)
    if matches:
        genes = db_df[db_df['Metabolite_Name'].astype(str).str.lower().isin(matches)]['Target'].dropna().unique()
        matched_genes.update(genes)
        print(f"  {m} -> {len(genes)} genes found")
    else:
        print(f"  {m} -> No matches found in DB")

matched_genes = list(matched_genes)
print(f"\\nTotal unique genes mapped to top altered metabolites: {len(matched_genes)}")
"""))

    cells.append(nbf.v4.new_markdown_cell("""## 3. Load Spatial Transcriptomics Data"""))
    cells.append(nbf.v4.new_code_cell("""# Process Spatial Transcriptomics Data
import zipfile
import subprocess
import glob

# 1. Unzip COSMx_Seurat_objects.zip or Visium_objects.zip if they exist
for zip_file in glob.glob(os.path.join(SPATIAL_DIR, "*.zip")):
    extract_dir = zip_file.replace(".zip", "")
    if not os.path.exists(extract_dir):
        print(f"Extracting {os.path.basename(zip_file)}...")
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)

# 2. Look for .mtx files
mtx_files = glob.glob(os.path.join(SPATIAL_DIR, "**", "*.mtx"), recursive=True)

# 3. If no .mtx, but .rds exists, attempt conversion
if len(mtx_files) == 0:
    rds_files = glob.glob(os.path.join(SPATIAL_DIR, "**", "*.rds"), recursive=True)
    if len(rds_files) > 0:
        print(f"Found {len(rds_files)} .rds files. Will attempt to extract matrices using R...")
        for target_rds in rds_files:
            print(f"\\nAttempting to extract {os.path.basename(target_rds)}...")
            base_path = target_rds.replace('.rds', '')
            target_mtx = base_path + ".mtx"
            
            r_code = f\"\"\"
            library(Seurat)
            library(Matrix)
            seu <- readRDS("{target_rds}")
            
            assay_obj <- seu@assays[[DefaultAssay(seu)]]
            counts_mat <- tryCatch({{ LayerData(seu, layer="counts") }}, error = function(e) {{
                if (.hasSlot(assay_obj, "counts")) assay_obj@counts else assay_obj$counts
            }})
            
            if (inherits(counts_mat, "BPCells")) {{
                counts_mat <- as(counts_mat, "dgCMatrix")
            }}
            
            writeMM(counts_mat, "{base_path}.mtx")
            write.csv(colnames(counts_mat), "{base_path}.barcodes.csv", row.names=FALSE)
            write.csv(rownames(counts_mat), "{base_path}.features.csv", row.names=FALSE)
            write.csv(seu@meta.data, "{base_path}.meta.csv", row.names=TRUE)
            
            for (red in names(seu@reductions)) {{
                write.csv(seu@reductions[[red]]@cell.embeddings, paste0("{base_path}.obsm.", red, ".csv"), row.names=FALSE)
            }}
            \"\"\"
            
            script_path = os.path.join(SPATIAL_DIR, "convert_rds.R")
            with open(script_path, "w") as f:
                f.write(r_code)
                
            try:
                print("Running Rscript convert_rds.R ...")
                subprocess.run(["Rscript", script_path], check=True)
                if os.path.exists(target_mtx):
                    mtx_files = [target_mtx]
                    print(f"Successfully extracted {os.path.basename(target_rds)} to .mtx!")
                    break
                else:
                    print(f"R script finished for {os.path.basename(target_rds)}, but .mtx was not created. Trying next...")
            except Exception as e:
                print(f"R extraction failed for {os.path.basename(target_rds)}. Error: {e}")
                continue
                
        if len(mtx_files) == 0:
            import sys
            sys.exit("CRITICAL ERROR - All R extraction attempts failed.")
"""))

# 4. Load all available spatial datasets (Visium + MTX)
    cells.append(nbf.v4.new_code_cell("""
import scipy.sparse as sp
import pandas as pd
import anndata as ad
import os
import glob

adatas = {}

# 4A. Find and load Visium objects
print("Looking for Visium directories...")
for root, dirs, files in os.walk(SPATIAL_DIR):
    if '__MACOSX' in root:
        continue
    if 'spatial' in dirs and 'filtered_feature_bc_matrix' in dirs:
        name = os.path.basename(root)
        print(f"\\nLoading Visium object: {name}")
        try:
            ad_obj = sc.read_visium(root)
            ad_obj.var_names_make_unique()
            adatas[name] = ad_obj
        except Exception as e:
            print(f"Failed to load Visium {name}: {e}")
            raise RuntimeError(f"CRITICAL FAILURE: Failed to load Visium {name}. Halting execution as per invariant.") from e

# 4B. Find and load MTX objects (from CosMx/Nanostring RDS extraction)
mtx_files = glob.glob(os.path.join(SPATIAL_DIR, "**", "*.mtx"), recursive=True)
for mtx in mtx_files:
    name = os.path.basename(mtx).replace('.mtx', '')
    print(f"\\nLoading extracted MTX object: {name}")
    try:
        base_path = mtx.replace('.mtx', '')
        counts = sc.read_mtx(mtx).T
        counts.obs_names = pd.read_csv(f"{base_path}.barcodes.csv", header=0).iloc[:, 0].values
        counts.var_names = pd.read_csv(f"{base_path}.features.csv", header=0).iloc[:, 0].values
        counts.obs = pd.read_csv(f"{base_path}.meta.csv", index_col=0)
        
        for obsm_file in glob.glob(f"{base_path}.obsm.*.csv"):
            red_name = obsm_file.split('.obsm.')[1].replace('.csv', '')
            coords = pd.read_csv(obsm_file).values
            counts.obsm[f'X_{red_name}'] = coords
            if red_name.lower() in ['spatial', 'fov']:
                counts.obsm['spatial'] = coords
        
        if 'spatial' not in counts.obsm:
            if 'CenterX_global_px' in counts.obs and 'CenterY_global_px' in counts.obs:
                counts.obsm['spatial'] = counts.obs[['CenterX_global_px', 'CenterY_global_px']].values
            elif 'CenterX_local_px' in counts.obs and 'CenterY_local_px' in counts.obs:
                counts.obsm['spatial'] = counts.obs[['CenterX_local_px', 'CenterY_local_px']].values
            elif 'x' in counts.obs and 'y' in counts.obs:
                counts.obsm['spatial'] = counts.obs[['x', 'y']].values
        
        counts.var_names_make_unique()
        adatas[name] = counts
    except Exception as e:
        print(f"Failed to load MTX {name}: {e}")
        raise RuntimeError(f"CRITICAL FAILURE: Failed to load MTX {name}. Halting execution as per invariant.") from e

# Basic QC and normalization for all raw counts
for name, adata in adatas.items():
    if 'n_genes_by_counts' not in adata.obs:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    if adata.X.max() > 1000: # Heuristic for raw counts
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    print(f"[{name}] Ready: {adata.shape[0]} spots/cells, {adata.shape[1]} genes")

if len(adatas) == 0:
    print(f"No usable spatial datasets found in {SPATIAL_DIR}. Please check the downloaded ZIPs and extraction.")
"""))

    cells.append(nbf.v4.new_markdown_cell("""## 4. Spatial Inferred Metabolomics (Module Score)"""))
    cells.append(nbf.v4.new_code_cell("""
if len(adatas) > 0 and len(matched_genes) > 0:
    for name, adata in adatas.items():
        print(f"\\n{'='*50}\\nProcessing {name}\\n{'='*50}")
        # Filter for genes actually present in this spatial dataset
        present_genes = [g for g in matched_genes if g in adata.var_names]
        print(f"Genes present in spatial data: {len(present_genes)} / {len(matched_genes)}")
        
        if len(present_genes) > 0:
            # Score the spots based on the metabolic gene signature
            sc.tl.score_genes(adata, gene_list=present_genes, score_name='Metabolic_Activity_Score')
            
            # Determine plotting method based on available spatial coordinates
            if 'spatial' in adata.obsm:
                sc.pl.spatial(adata, color='Metabolic_Activity_Score', cmap='magma', title=f'Metabolic Activity ({name})', show=False)
            elif 'X_spatial' in adata.obsm:
                sc.pl.embedding(adata, basis='spatial', color='Metabolic_Activity_Score', cmap='magma', title=f'Metabolic Activity ({name})', show=False)
            elif 'X_umap' in adata.obsm:
                sc.pl.umap(adata, color='Metabolic_Activity_Score', cmap='magma', title=f'UMAP Inferred Metabolic Activity ({name})', show=False)
            else:
                print("No spatial or UMAP coordinates found in adata.obsm.")
                
            SAVE_AS_HTML = True # Placeholder for execute_pancancer_notebooks to strip
            plt.show()
        else:
            print("None of the mapped genes are present in this spatial dataset.")
else:
    print("Cannot proceed with spatial plotting due to missing data or genes.")

"""))

    nb['cells'] = cells
    
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kidney_spatial_massspec_integration.ipynb')
    with open(output_path, 'w') as f:
        nbf.write(nb, f)
        
    print(f"Generated notebook at {output_path}")

if __name__ == '__main__':
    generate_notebook()
