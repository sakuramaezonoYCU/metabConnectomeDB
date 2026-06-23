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
import scipy.stats as stats
import argparse

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

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

def analyze_sample(sample_dir, sample_name, out_dir, sig_genes, sig_name):
    print(f"\n[{sample_name}] Loading Visium data...")
    adata = load_visium_from_mtx(sample_dir)
        
    adata.var_names_make_unique()
    
    # Compute signatures
    print(f"[{sample_name}] Computing spatial signatures...")
    if len(sig_genes) == 0:
        print(f"  -> Signature has 0 genes. Skipping analysis for {sample_name}.")
        return {
            'Sample': sample_name,
            'Spots': adata.n_obs,
            'Morans_I': np.nan,
            'Morans_Pval': np.nan
        }
    adata.obs['Signature_Score'] = calculate_module_score(adata, sig_genes, 'temp_sig')
    
    # Calculate spatial neighbors
    print(f"[{sample_name}] Building spatial graph...")
    sq.gr.spatial_neighbors(adata)
    
    # Generate plots
    print(f"[{sample_name}] Generating spatial plots...")
    sc.pl.spatial(adata, color=['Signature_Score'], cmap='viridis', 
                  title=[f'{sample_name}: {sig_name} Score'], 
                  show=False)
                  
    plot_path = os.path.join(out_dir, f"visium_{sample_name}_{sig_name}.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Calculate spot-to-spot spatial autocorrelation (Moran's I) to measure spatial clustering
    # Squidpy spatial_autocorr expects genes in adata.X, so we create a temp AnnData
    tmp_adata = sc.AnnData(X=adata.obs[['Signature_Score']].values)
    tmp_adata.obsp['spatial_connectivities'] = adata.obsp['spatial_connectivities']
    tmp_adata.obsp['spatial_distances'] = adata.obsp['spatial_distances']
    tmp_adata.var_names = ['Signature_Score']
    
    sq.gr.spatial_autocorr(tmp_adata, mode='moran', genes=['Signature_Score'])
    moran_df = tmp_adata.uns['moranI']
    moran_i = moran_df.loc['Signature_Score', 'I']
    
    # Squidpy versions differ in their p-value column names (e.g., pval_norm, pval_sim)
    pval_cols = [c for c in moran_df.columns if 'pval' in c.lower()]
    moran_pval = moran_df.loc['Signature_Score', pval_cols[0]] if pval_cols else np.nan
    
    print(f"[{sample_name}] Moran's I (Spatial Clustering): I={moran_i:.3f}, p={moran_pval:.3e}")
    
    return {
        'Sample': sample_name,
        'Spots': adata.n_obs,
        'Morans_I': moran_i,
        'Morans_Pval': moran_pval
    }

def main():
    parser = argparse.ArgumentParser(description="Visium Spatial Transcriptomics Signature Verification")
    parser.add_argument('--signature_csv', required=True, help="Path to the signature CSV file")
    args = parser.parse_args()
    
    if not os.path.exists(args.signature_csv):
        raise FileNotFoundError(f"CRITICAL ERROR: Signature CSV {args.signature_csv} does not exist.")
        
    sig_name = os.path.basename(args.signature_csv).replace('.csv', '')
        
    df_sig = pd.read_csv(args.signature_csv)
    if 'Strictly_Conserved_Gene' in df_sig.columns:
        sig_genes = df_sig['Strictly_Conserved_Gene'].dropna().unique().tolist()
    elif 'Gene' in df_sig.columns:
        sig_genes = df_sig['Gene'].dropna().unique().tolist()
    elif 'Target' in df_sig.columns:
        sig_genes = df_sig['Target'].dropna().unique().tolist()
    else:
        raise ValueError(f"CRITICAL ERROR: Could not find gene column in {args.signature_csv}")

    raw_dir = os.path.join(BASE_DIR, "input", "spatial", "GSE211956_RAW_Forrest")
    out_dir = os.path.join(BASE_DIR, "output", "spatial_verification", sig_name)
    os.makedirs(out_dir, exist_ok=True)
    
    if not os.path.exists(raw_dir):
        print(f"Warning: Directory {raw_dir} not found. Skipping spatial verification.")
        return
        
    samples = [f"SP{i}" for i in range(1, 9)]
    results = []
    
    for sp in samples:
        print(f"\nProcessing {sp} for {sig_name}...")
        tmp_dir = format_visium_dir(raw_dir, sp, out_dir)
        if tmp_dir:
            res = analyze_sample(tmp_dir, sp, out_dir, sig_genes, sig_name)
            if res:
                results.append(res)
                
    if results:
        df_res = pd.DataFrame(results)
        csv_path = os.path.join(out_dir, f"visium_spatial_clustering_summary.csv")
        df_res.to_csv(csv_path, index=False)
        print(f"\nSuccess! Analyzed {len(results)} samples. Summary saved to {csv_path}.")
        print("Generated spatial plots are in the output directory.")
    else:
        print("No samples were successfully analyzed.")

if __name__ == "__main__":
    main()
