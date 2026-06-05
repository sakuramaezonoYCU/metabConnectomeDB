import json
import os

NOTEBOOK_PATH = "/Users/sakuramaezono/Library/CloudStorage/OneDrive-YokohamaCityUniversity/Personal/05_Python_repositories/metabConnectomeDB/scripts/visium_spatial_validation.ipynb"

def create_markdown_cell(source):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [line + "\n" for line in source]
    }

def create_code_cell(source):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [line + "\n" for line in source]
    }

cells = []

# Cell 0: Markdown - Title
cells.append(create_markdown_cell([
    "# True Spatial Validation (Visium GSE211956)",
    "",
    "This notebook validates the Serotonin-HTR7-TAM -> IP4 -> HR repair spatial axis using true High-Grade Serous Ovarian Cancer (HGSOC) Visium spatial transcriptomics data.",
    "",
    "**⚠️ IMPORTANT ONEDRIVE WARNING:** If you get `[Errno 60] Operation timed out` when running these cells, it means your Visium `.gz` files are still 'Cloud-only'. You MUST open Finder, go to `input/spatial/GSE211956_RAW_Forrest`, right-click the folder, and select **'Always keep on this device'** to force OneDrive to download them to your physical disk."
]))

# Cell 1: Code - Configuration and Setup
cells.append(create_code_cell([
    "import os",
    "import sys",
    "import glob",
    "import zipfile",
    "import json",
    "import numpy as np",
    "import pandas as pd",
    "import scanpy as sc",
    "import squidpy as sq",
    "import matplotlib.pyplot as plt",
    "from matplotlib.image import imread",
    "import scipy.stats as stats",
    "",
    "base_dir = '../'",
    "if base_dir not in sys.path:",
    "    sys.path.append(base_dir)",
    "",
    "from serotonin_config import HTR7_TAM_SIGNATURE, HR_REPAIR_GENES",
    "",
    "raw_dir = os.path.join(base_dir, 'input', 'spatial', 'GSE211956_RAW_Forrest')",
    "out_dir = os.path.join(base_dir, 'output', 'serotonin_axis_spatial_mapping')",
    "os.makedirs(out_dir, exist_ok=True)"
]))

# Cell 2: Code - Custom Loaders
cells.append(create_code_cell([
    "def format_visium_dir(base_raw_dir, sample_prefix, out_dir):",
    "    tmp_dir = os.path.join(out_dir, 'tmp_visium_loaders', f'tmp_10x_{sample_prefix}')",
    "    os.makedirs(tmp_dir, exist_ok=True)",
    "    ",
    "    files_to_link = {",
    "        'matrix.mtx.gz': 'matrix.mtx.gz',",
    "        'features.tsv.gz': 'features.tsv.gz',",
    "        'barcodes.tsv.gz': 'barcodes.tsv.gz'",
    "    }",
    "    ",
    "    for f_target, f_dest in files_to_link.items():",
    "        pattern = os.path.join(base_raw_dir, f'*{sample_prefix}*{f_target}')",
    "        matches = glob.glob(pattern)",
    "        if not matches: return None",
    "        src = matches[0]",
    "        dst = os.path.join(tmp_dir, f_dest)",
    "        if not os.path.exists(dst):",
    "            os.symlink(src, dst)",
    "            ",
    "    spatial_dir = os.path.join(tmp_dir, 'spatial')",
    "    if not os.path.exists(spatial_dir):",
    "        zip_pattern = os.path.join(base_raw_dir, f'*{sample_prefix}*spatial.zip')",
    "        zip_matches = glob.glob(zip_pattern)",
    "        if not zip_matches: return None",
    "        with zipfile.ZipFile(zip_matches[0], 'r') as zip_ref:",
    "            zip_ref.extractall(tmp_dir)",
    "            ",
    "    return tmp_dir",
    "",
    "def load_visium_from_mtx(sample_dir):",
    "    # 1. Read counts",
    "    adata = sc.read_10x_mtx(sample_dir)",
    "    ",
    "    # 2. Read spatial coordinates",
    "    spatial_dir = os.path.join(sample_dir, 'spatial')",
    "    tissue_pos_file = os.path.join(spatial_dir, 'tissue_positions.csv')",
    "    has_header = True",
    "    if not os.path.exists(tissue_pos_file):",
    "        tissue_pos_file = os.path.join(spatial_dir, 'tissue_positions_list.csv')",
    "        has_header = False",
    "        ",
    "    positions = pd.read_csv(tissue_pos_file, header=0 if has_header else None, index_col=0)",
    "    positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']",
    "    ",
    "    adata.obs = adata.obs.join(positions, how='left')",
    "    adata.obsm['spatial'] = adata.obs[['pxl_row_in_fullres', 'pxl_col_in_fullres']].to_numpy()",
    "    ",
    "    # 3. Read images",
    "    library_id = 'sample'",
    "    adata.uns['spatial'] = {library_id: {'images': {}, 'scalefactors': {}, 'metadata': {}}}",
    "    ",
    "    hires_path = os.path.join(spatial_dir, 'tissue_hires_image.png')",
    "    if os.path.exists(hires_path): adata.uns['spatial'][library_id]['images']['hires'] = imread(hires_path)",
    "        ",
    "    lowres_path = os.path.join(spatial_dir, 'tissue_lowres_image.png')",
    "    if os.path.exists(lowres_path): adata.uns['spatial'][library_id]['images']['lowres'] = imread(lowres_path)",
    "        ",
    "    scale_path = os.path.join(spatial_dir, 'scalefactors_json.json')",
    "    if os.path.exists(scale_path):",
    "        with open(scale_path, 'r') as f:",
    "            adata.uns['spatial'][library_id]['scalefactors'] = json.load(f)",
    "            ",
    "    return adata",
    "",
    "def calculate_module_score(adata, gene_list, score_name):",
    "    valid_genes = [g for g in gene_list if g in adata.var_names]",
    "    if len(valid_genes) == 0:",
    "        return np.zeros(adata.n_obs)",
    "    sc.tl.score_genes(adata, gene_list=valid_genes, score_name=score_name)",
    "    score = adata.obs[score_name].values",
    "    del adata.obs[score_name]",
    "    return score"
]))

# Cell 3: Markdown - Sample Selection
cells.append(create_markdown_cell([
    "## Process All Samples",
    "This loop will process all 8 spatial samples (`SP1` through `SP8`), compute the signatures, and render their respective colocalization plots."
]))

# Cell 4: Code - Loop over all samples
cells.append(create_code_cell([
    "samples = [f'SP{i}' for i in range(1, 9)]",
    "",
    "for sample_name in samples:",
    "    print(f\"\\n{'='*50}\")",
    "    print(f\"Processing {sample_name}...\")",
    "    tmp_dir = format_visium_dir(raw_dir, sample_name, out_dir)",
    "    if not tmp_dir:",
    "        print(f\"Could not format directory for {sample_name}. Skipping.\")",
    "        continue",
    "        ",
    "    print(f\"Loading matrix and spatial metadata...\")",
    "    try:",
    "        adata = load_visium_from_mtx(tmp_dir)",
    "    except Exception as e:",
    "        print(f\"Failed to load {sample_name}: {e}\")",
    "        continue",
    "        ",
    "    adata.var_names_make_unique()",
    "    print(f\"Loaded {adata.n_obs} spots and {adata.n_vars} genes.\")",
    "    ",
    "    print(\"Scoring HTR7+ TAMs and HR Repair signatures...\")",
    "    adata.obs['HTR7_TAM_Score'] = calculate_module_score(adata, HTR7_TAM_SIGNATURE, 'temp_tam')",
    "    adata.obs['HR_Repair_Score'] = calculate_module_score(adata, HR_REPAIR_GENES, 'temp_hr')",
    "    ",
    "    print(\"Building spatial kNN graph...\")",
    "    sq.gr.spatial_neighbors(adata)",
    "    ",
    "    # Plot signatures",
    "    sc.pl.spatial(adata, color=['HTR7_TAM_Score', 'HR_Repair_Score'], cmap='viridis', ",
    "                  title=[f'{sample_name}: HTR7+ TAM Score', f'{sample_name}: HR Repair Score'], ",
    "                  size=1.5)",
    "                  ",
    "    # Co-localization analysis",
    "    r, p = stats.pearsonr(adata.obs['HTR7_TAM_Score'], adata.obs['HR_Repair_Score'])",
    "    print(f\"Spot-level Pearson correlation: r={r:.3f}, p={p:.3e}\")",
    "    ",
    "    thresh_tam = np.percentile(adata.obs['HTR7_TAM_Score'], 80)",
    "    thresh_hr = np.percentile(adata.obs['HR_Repair_Score'], 80)",
    "    ",
    "    adata.obs['High_TAM'] = adata.obs['HTR7_TAM_Score'] > thresh_tam",
    "    adata.obs['High_HR'] = adata.obs['HR_Repair_Score'] > thresh_hr",
    "    ",
    "    def assign_category(row):",
    "        if row['High_TAM'] and row['High_HR']: return 'Co-localized'",
    "        elif row['High_TAM']: return 'TAM_Only'",
    "        elif row['High_HR']: return 'HR_Only'",
    "        return 'Neither'",
    "        ",
    "    adata.obs['Spot_Category'] = adata.obs.apply(assign_category, axis=1).astype('category')",
    "    ",
    "    palette = {",
    "        'Neither': '#e0e0e0', ",
    "        'TAM_Only': '#1f77b4', ",
    "        'HR_Only': '#ff7f0e', ",
    "        'Co-localized': '#d62728'",
    "    }",
    "    ",
    "    sc.pl.spatial(adata, color='Spot_Category', palette=palette,",
    "                  title=f'{sample_name}: High TAM / High HR Repair Co-localization',",
    "                  size=1.5)",
    "    plt.show()"
]))


notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3 (ipykernel)",
            "language": "python",
            "name": "python3"
        },
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "nbconvert_exporter": "python",
            "pygments_lexer": "ipython3",
            "version": "3.10.0"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 5
}

with open(NOTEBOOK_PATH, "w", encoding="utf-8") as f:
    json.dump(notebook, f, indent=1)

print(f"Successfully generated notebook at {NOTEBOOK_PATH}")
