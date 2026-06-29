import os
import sys
import glob
import zipfile
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.stats as stats

# Config
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if BASE_DIR not in sys.path:
    sys.path.append(BASE_DIR)

import serotonin_config as cfg

out_dir = os.path.join(BASE_DIR, "output", "serotonin_axis_spatial_mapping")
os.makedirs(out_dir, exist_ok=True)

def format_visium_dir(base_raw_dir, sample_prefix):
    """Formats Visium symlinks from raw Geo data."""
    tmp_dir = os.path.join(out_dir, "tmp_visium_loaders", f"tmp_10x_{sample_prefix}")
    os.makedirs(tmp_dir, exist_ok=True)
    
    files_to_link = {
        "matrix.mtx.gz": "matrix.mtx.gz",
        "features.tsv.gz": "features.tsv.gz",
        "barcodes.tsv.gz": "barcodes.tsv.gz"
    }
    
    for f_target, f_dest in files_to_link.items():
        pattern = os.path.join(base_raw_dir, f"*{sample_prefix}*{f_target}")
        matches = glob.glob(pattern)
        if not matches:
            return None
        src = matches[0]
        dst = os.path.join(tmp_dir, f_dest)
        if not os.path.exists(dst):
            os.symlink(src, dst)
            
    spatial_dir = os.path.join(tmp_dir, "spatial")
    if not os.path.exists(spatial_dir):
        zip_pattern = os.path.join(base_raw_dir, f"*{sample_prefix}*spatial.zip")
        zip_matches = glob.glob(zip_pattern)
        if not zip_matches:
            return None
        with zipfile.ZipFile(zip_matches[0], 'r') as zip_ref:
            zip_ref.extractall(tmp_dir)
            
    return tmp_dir

def load_visium_from_mtx(sample_dir):
    """Loads Visium from MTX and extracts coordinates."""
    import json
    from matplotlib.image import imread
    
    adata = sc.read_10x_mtx(sample_dir)
    spatial_dir = os.path.join(sample_dir, "spatial")
    
    tissue_pos_file = os.path.join(spatial_dir, "tissue_positions.csv")
    has_header = True
    if not os.path.exists(tissue_pos_file):
        tissue_pos_file = os.path.join(spatial_dir, "tissue_positions_list.csv")
        has_header = False
        
    positions = pd.read_csv(tissue_pos_file, header=0 if has_header else None, index_col=0)
    positions.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
    adata.obs = adata.obs.join(positions, how="left")
    adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
    
    library_id = "sample"
    adata.uns["spatial"] = {library_id: {"images": {}, "scalefactors": {}, "metadata": {}}}
    
    for key, filename in [("hires", "tissue_hires_image.png"), ("lowres", "tissue_lowres_image.png")]:
        path = os.path.join(spatial_dir, filename)
        if os.path.exists(path):
            adata.uns["spatial"][library_id]["images"][key] = imread(path)
            
    scale_path = os.path.join(spatial_dir, "scalefactors_json.json")
    if os.path.exists(scale_path):
        with open(scale_path, 'r') as f:
            adata.uns["spatial"][library_id]["scalefactors"] = json.load(f)
            
    return adata

def safely_score_genes(adata, gene_list, score_name):
    """Scores genes safely."""
    valid_genes = [g for g in gene_list if g in adata.var_names]
    if not valid_genes:
        adata.obs[score_name] = 0.0
        return
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=score_name)

def main():
    raw_dir = os.path.join(BASE_DIR, "input", "spatial", "GSE211956_RAW_Forrest")
    if not os.path.exists(raw_dir):
        print(f"Error: Directory {raw_dir} not found.")
        return
        
    samples = [f"SP{i}" for i in range(1, 9)]
    results = []
    
    for sp in samples:
        print(f"\\nProcessing {sp}...")
        tmp_dir = format_visium_dir(raw_dir, sp)
        if not tmp_dir:
            print(f"  -> Format failed for {sp}")
            continue
            
        try:
            adata = load_visium_from_mtx(tmp_dir)
            adata.var_names_make_unique()
            
            # Combine the markers into two logical Niche signatures
            suppressive_niche_genes = cfg.HTR7_TAM_SIGNATURE + cfg.SUPPRESSIVE_LIGAND_TARGETS
            exhaustion_niche_genes = cfg.EXHAUSTION_TARGETS
            
            safely_score_genes(adata, suppressive_niche_genes, 'Suppressive_Niche_Score')
            safely_score_genes(adata, exhaustion_niche_genes, 'Exhaustion_Niche_Score')
            
            # Filter tissue spots only
            adata = adata[adata.obs['in_tissue'] == 1].copy()
            
            # Calculate correlation
            if len(adata.obs) > 0:
                r, p = stats.pearsonr(adata.obs['Suppressive_Niche_Score'], adata.obs['Exhaustion_Niche_Score'])
                print(f"  -> Pearson R: {r:.3f}, p-value: {p:.3e}")
                
                # Plot
                sc.pl.spatial(adata, color=['Suppressive_Niche_Score', 'Exhaustion_Niche_Score'], 
                              cmap='magma', show=False,
                              title=[f'{sp}: HTR7+ Suppressive Niche', f'{sp}: Exhaustion Niche'])
                plot_path = os.path.join(out_dir, f"{sp}_immune_evasion_niches.png")
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                results.append({
                    'Sample': sp,
                    'Spots': adata.n_obs,
                    'Pearson_R': r,
                    'Pearson_P': p
                })
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"Failed {sp}: {e}")
            raise
            
    if results:
        df_res = pd.DataFrame(results)
        df_res.to_csv(os.path.join(out_dir, 'visium_immune_evasion_summary.csv'), index=False)
        print("\\nDone! Spatial summaries saved.")

if __name__ == '__main__':
    main()
