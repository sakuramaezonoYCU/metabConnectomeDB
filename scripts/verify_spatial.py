"""
verify_spatial.py
=================
Validates the HTR7+ TAM -> IP4 -> HR repair spatial axis using true 
Visium spatial transcriptomics data from HGSOC (GSE211956).
"""
import os
import sys
import glob
import shutil
import zipfile
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

from serotonin_config import HTR7_TAM_SIGNATURE, HR_REPAIR_GENES

def format_visium_dir(base_raw_dir, sample_prefix, out_dir):
    """
    sc.read_visium requires a strict folder structure:
    dir/
      matrix.mtx.gz
      features.tsv.gz
      barcodes.tsv.gz
      spatial/
        tissue_positions_list.csv (or tissue_positions.csv)
        scalefactors_json.json
        tissue_hires_image.png
        tissue_lowres_image.png
        
    GEO data often prepends the accession (e.g., GSM6506110_SP1_matrix.mtx.gz).
    This function creates a temporary symlinked directory matching the 10x format.
    """
    tmp_dir = os.path.join(out_dir, "tmp_visium_loaders", f"tmp_10x_{sample_prefix}")
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Symlink matrix files
    files_to_link = {
        "matrix.mtx.gz": "matrix.mtx.gz",
        "features.tsv.gz": "features.tsv.gz",
        "barcodes.tsv.gz": "barcodes.tsv.gz"
    }
    
    for f_target, f_dest in files_to_link.items():
        # Find the matching file in the raw dir
        pattern = os.path.join(base_raw_dir, f"*{sample_prefix}*{f_target}")
        matches = glob.glob(pattern)
        if not matches:
            print(f"  -> Missing {f_target} for {sample_prefix}")
            return None
            
        src = matches[0]
        dst = os.path.join(tmp_dir, f_dest)
        if not os.path.exists(dst):
            os.symlink(src, dst)
            
    # Handle the spatial zip
    spatial_dir = os.path.join(tmp_dir, "spatial")
    if not os.path.exists(spatial_dir):
        zip_pattern = os.path.join(base_raw_dir, f"*{sample_prefix}*spatial.zip")
        zip_matches = glob.glob(zip_pattern)
        if not zip_matches:
            print(f"  -> Missing spatial.zip for {sample_prefix}")
            return None
            
        with zipfile.ZipFile(zip_matches[0], 'r') as zip_ref:
            zip_ref.extractall(tmp_dir) # Often zips contain the 'spatial' folder already
            
        # Check if it extracted a 'spatial' folder
        if not os.path.exists(spatial_dir):
            print(f"  -> Warning: Extracted zip for {sample_prefix} did not create 'spatial' folder.")
            return None
            
    return tmp_dir

def calculate_module_score(adata, gene_list, score_name):
    valid_genes = [g for g in gene_list if g in adata.var_names]
    if len(valid_genes) == 0:
        return np.zeros(adata.n_obs)
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=score_name)
    score = adata.obs[score_name].values
    del adata.obs[score_name]
    return score

def load_visium_from_mtx(sample_dir):
    """
    Custom loader for Visium data from matrix.mtx.gz format, 
    as sc.read_visium relies on the .h5 file format.
    """
    import json
    from matplotlib.image import imread
    
    # 1. Read counts
    adata = sc.read_10x_mtx(sample_dir)
    
    # 2. Read spatial coordinates
    spatial_dir = os.path.join(sample_dir, "spatial")
    tissue_pos_file = os.path.join(spatial_dir, "tissue_positions.csv")
    has_header = True
    if not os.path.exists(tissue_pos_file):
        tissue_pos_file = os.path.join(spatial_dir, "tissue_positions_list.csv")
        has_header = False
        
    positions = pd.read_csv(
        tissue_pos_file,
        header=0 if has_header else None,
        index_col=0
    )
    positions.columns = [
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]
    
    # Join positions to obs (barcodes are indices)
    adata.obs = adata.obs.join(positions, how="left")
    
    # Set obsm spatial coordinates
    adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
    
    # 3. Read images and scalefactors
    library_id = "sample"
    adata.uns["spatial"] = {library_id: {"images": {}, "scalefactors": {}, "metadata": {}}}
    
    hires_path = os.path.join(spatial_dir, "tissue_hires_image.png")
    if os.path.exists(hires_path):
        adata.uns["spatial"][library_id]["images"]["hires"] = imread(hires_path)
        
    lowres_path = os.path.join(spatial_dir, "tissue_lowres_image.png")
    if os.path.exists(lowres_path):
        adata.uns["spatial"][library_id]["images"]["lowres"] = imread(lowres_path)
        
    scale_path = os.path.join(spatial_dir, "scalefactors_json.json")
    if os.path.exists(scale_path):
        with open(scale_path, 'r') as f:
            adata.uns["spatial"][library_id]["scalefactors"] = json.load(f)
            
    return adata

def analyze_sample(sample_dir, sample_name, out_dir):
    print(f"\n[{sample_name}] Loading Visium data...")
    try:
        adata = load_visium_from_mtx(sample_dir)
    except Exception as e:
        import traceback
        print(f"  -> Failed to load {sample_name}:")
        traceback.print_exc()
        return
        
    adata.var_names_make_unique()
    
    # Compute signatures
    print(f"[{sample_name}] Computing spatial signatures...")
    adata.obs['HTR7_TAM_Score'] = calculate_module_score(adata, HTR7_TAM_SIGNATURE, 'temp_tam')
    adata.obs['HR_Repair_Score'] = calculate_module_score(adata, HR_REPAIR_GENES, 'temp_hr')
    
    # Calculate spatial neighbors
    print(f"[{sample_name}] Building spatial graph...")
    sq.gr.spatial_neighbors(adata)
    
    # Generate plots
    print(f"[{sample_name}] Generating spatial plots...")
    sc.pl.spatial(adata, color=['HTR7_TAM_Score', 'HR_Repair_Score'], cmap='viridis', 
                  title=[f'{sample_name}: HTR7+ TAM Score', f'{sample_name}: HR Repair Score'], 
                  show=False)
                  
    plot_path = os.path.join(out_dir, f"visium_{sample_name}_signatures.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate spot-to-spot correlation
    r, p = stats.pearsonr(adata.obs['HTR7_TAM_Score'], adata.obs['HR_Repair_Score'])
    print(f"[{sample_name}] Spot-level correlation: r={r:.3f}, p={p:.3e}")
    
    # Binarize scores (top 20% spots) for co-occurrence analysis
    thresh_tam = np.percentile(adata.obs['HTR7_TAM_Score'], 80)
    thresh_hr = np.percentile(adata.obs['HR_Repair_Score'], 80)
    
    adata.obs['High_TAM'] = adata.obs['HTR7_TAM_Score'] > thresh_tam
    adata.obs['High_HR'] = adata.obs['HR_Repair_Score'] > thresh_hr
    
    # Create categorical for co-occurrence
    def assign_category(row):
        if row['High_TAM'] and row['High_HR']: return 'Co-localized'
        elif row['High_TAM']: return 'TAM_Only'
        elif row['High_HR']: return 'HR_Only'
        return 'Neither'
        
    adata.obs['Spot_Category'] = adata.obs.apply(assign_category, axis=1).astype('category')
    
    sc.pl.spatial(adata, color='Spot_Category', 
                  title=f'{sample_name}: High TAM / High HR Repair Co-localization',
                  show=False)
    coloc_path = os.path.join(out_dir, f"visium_{sample_name}_colocalization.png")
    plt.savefig(coloc_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        'Sample': sample_name,
        'Spots': adata.n_obs,
        'Pearson_R': r,
        'Pearson_P': p,
        'Co_localized_Spots': (adata.obs['Spot_Category'] == 'Co-localized').sum()
    }

def main():
    import scipy.stats as stats
    # Make stats available globally for the analyze function
    global stats 
    
    raw_dir = os.path.join(BASE_DIR, "input", "spatial", "GSE211956_RAW_Forrest")
    out_dir = os.path.join(BASE_DIR, "output", "serotonin_axis_spatial_mapping")
    os.makedirs(out_dir, exist_ok=True)
    
    if not os.path.exists(raw_dir):
        print(f"Error: Directory {raw_dir} not found.")
        print("Please ensure the GSE211956 spatial data is downloaded.")
        return
        
    samples = [f"SP{i}" for i in range(1, 9)]
    results = []
    
    for sp in samples:
        print(f"\nProcessing {sp}...")
        tmp_dir = format_visium_dir(raw_dir, sp, out_dir)
        if tmp_dir:
            res = analyze_sample(tmp_dir, sp, out_dir)
            if res:
                results.append(res)
                
    if results:
        df_res = pd.DataFrame(results)
        csv_path = os.path.join(out_dir, "visium_colocalization_summary.csv")
        df_res.to_csv(csv_path, index=False)
        print(f"\nSuccess! Analyzed {len(results)} samples. Summary saved to {csv_path}.")
        print("Generated spatial plots are in the output directory.")
    else:
        print("No samples were successfully analyzed.")

if __name__ == "__main__":
    main()
