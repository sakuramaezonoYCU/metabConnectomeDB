#!/usr/bin/env python3
"""Generate cancer_cellxgene_integration.ipynb notebook programmatically."""
import json, os

def md(source): return {"cell_type":"markdown","metadata":{},"source":source.splitlines(keepends=True)}
def code(source): return {"cell_type":"code","execution_count":None,"metadata":{},"outputs":[],"source":source.splitlines(keepends=True)}

cells = []

# ── Section 1: Setup ──
cells.append(md("""# Cancer Microenvironment × MetabConnectomeDB Integration

This notebook bridges **MetabConnectomeDB metabolite-target gene pairs** with cancer-specific single-cell RNA-seq datasets from [CellxGene Census](https://cellxgene.cziscience.com/datasets).

## Goals
1. Query CellxGene Census for cancer scRNA-seq data filtered by disease, tissue, and organism
2. Map MetabConnectomeDB target genes onto cell-type-resolved expression profiles
3. Infer intercellular metabolic signaling using **LIANA+** cell-cell communication analysis
4. Contextualize findings within known cancer metabolic pathways (IDO1/Kynurenine, xCT/Glutamate, CD73/Adenosine, etc.)

> **Design:** All parameters are set in Section 1. Changing `DISEASE_FILTER` or `TISSUE_FILTER` re-runs the entire analysis without code changes."""))

cells.append(md("## 1. Setup & User Parameters"))

cells.append(code("""# ============================================================
# USER PARAMETERS — Edit these to control the analysis scope
# ============================================================

import os

# Path to MetabConnectomeDB target pairs (relative to this notebook)
TARGET_PAIRS_CSV = os.path.join('..', 'output',
    'human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv')
output_dir     = os.path.join('..', 'output')
os.makedirs(output_dir, exist_ok=True)

# CellxGene Census filters
ORGANISM   = "Homo sapiens"           # Match organism in filename
DISEASE_FILTER = [                     # Cancer types to query
    "lung adenocarcinoma",
    "breast cancer",
    "colorectal cancer",
]
TISSUE_FILTER  = None                  # None = all tissues; or e.g. ["lung"]

# ANALYSIS MODE CONFIGURATION:
# - "whole_transcriptome": Downloads 10,000 cells for ALL ~60,000 genes (genome-wide background).
#   Highly recommended! Allows proper library size normalization and genome-wide Differential Expression (DE) to rank metabolic targets against all other genes.
# - "target_only": Downloads 50,000 cells for ONLY the ~500 target genes (faster query but restricted to target genes).
DOWNLOAD_MODE = "whole_transcriptome"  # "whole_transcriptome" or "target_only"

# Set CAP to an integer to cap cells (e.g. 10000), or None to query ALL available cells
# Can be overridden dynamically via CELLXGENE_CAP environment variable
CAP = os.environ.get("CELLXGENE_CAP", None)
if CAP is not None:
    try:
        CAP = int(CAP)
    except ValueError:
        pass
CENSUS_VERSION = "stable"              # "stable" or specific date string

# Helper to construct a dynamic, parameter-dependent filename for caching
import re, hashlib
def get_cache_filename(diseases, tissue, cap, download_mode, census_version):
    def slugify(text):
        if not text:
            return "all"
        if isinstance(text, list):
            return "_".join(slugify(t) for t in text)
        return re.sub(r'[^a-z0-9]+', '-', str(text).lower()).strip('-')
    
    disease_slug = slugify(diseases)
    if len(disease_slug) > 100:
        h = hashlib.md5("_".join(diseases).encode('utf-8')).hexdigest()[:6]
        disease_slug = f"{disease_slug[:90]}-{h}"
    tissue_slug = slugify(tissue) if tissue else "all"
    cells_str = f"{cap//1000}k" if (cap is not None and cap >= 1000) else (str(cap) if cap is not None else "all")
    return f"cancer_{disease_slug}_{tissue_slug}_{cells_str}_{download_mode}_{census_version}.h5ad"

h5ad_filename = get_cache_filename(DISEASE_FILTER, TISSUE_FILTER, CAP, DOWNLOAD_MODE, CENSUS_VERSION)
h5ad_path = os.path.join(output_dir, h5ad_filename)

# Plot style
%matplotlib inline
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style='whitegrid', context='notebook', font_scale=1.1)
plt.rcParams['figure.dpi'] = 120

print("✅ Parameters configured.")
print(f"   Organism:  {ORGANISM}")
print(f"   Diseases:  {DISEASE_FILTER}")
print(f"   Tissues:   {TISSUE_FILTER or 'all'}")
print(f"   Download Mode: {DOWNLOAD_MODE}")
print(f"   Cap:       {CAP if CAP is not None else 'None (All cells)'}")
print(f"   Cache path: {h5ad_path}")"""))

cells.append(code("""# Dependency check
import importlib, sys
required = ['pandas','numpy','scanpy','anndata','cellxgene_census']
missing = [m for m in required if importlib.util.find_spec(m) is None]
if missing:
    print(f"❌ Missing packages: {missing}")
    print("   Install with: pip install " + " ".join(
        m.replace('_','-') for m in missing))
    sys.exit(1)

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad

# Try importing liana (optional for CCC)
try:
    import liana
    LIANA_AVAILABLE = True
    print("✅ LIANA+ available for cell-cell communication analysis.")
except ImportError:
    LIANA_AVAILABLE = False
    print("⚠️  liana not installed. CCC section will be skipped.")
    print("   Install with: pip install liana")

print("✅ All core dependencies loaded.")"""))

# ── Section 2: Load target pairs ──
cells.append(md("""---
## 2. Load MetabConnectomeDB Target Gene Pairs

### Purpose
Load the curated metabolite-target gene pair dataset produced by our upstream pipeline. These pairs define the "metabolic communication vocabulary" — the set of genes through which metabolites exert their effects in the tumor microenvironment (TME).

### How to Interpret
- **Target genes** include receptors, enzymes, and transporters that physically interact with metabolites
- A high number of unique targets indicates broad metabolic communication capacity
- The database provenance breakdown shows which source databases contribute each interaction"""))

cells.append(code("""# Load target pairs and split composite targets (e.g. multi-subunit complexes)
target_df = pd.read_csv(TARGET_PAIRS_CSV)
print(f"Loaded {len(target_df):,} raw metabolite-target pairs")

# Explode comma- and semicolon-separated target gene symbols to get individual clean symbols
target_df['Target'] = target_df['Target'].astype(str).str.split(r'[,;]')
target_df = target_df.explode('Target')
target_df['Target'] = target_df['Target'].str.strip()
target_df = target_df[target_df['Target'] != 'nan'].dropna(subset=['Target'])

print(f"Exploded into {len(target_df):,} individual target-gene rows")
print(f"  Unique HMDB metabolites: {target_df['HMDB_ID'].nunique():,} ({target_df['Metabolite_Name'].nunique():,} names)")
print(f"  Unique target genes: {target_df['Target'].nunique():,}")

# Extract unique gene symbols for CellxGene query
target_genes = sorted(target_df['Target'].dropna().unique().tolist())
print(f"\\n→ {len(target_genes)} unique target gene symbols ready for Census query")

# Database provenance
if 'database' in target_df.columns:
    print("\\nDatabase provenance:")
    for db, cnt in target_df['database'].value_counts().items():
        print(f"  {db}: {cnt:,}")"""))

# ── Section 3: Census metadata exploration ──
cells.append(md("""---
## 3. CellxGene Census Metadata Exploration

### Purpose
Before downloading expression data (which can be large), we first query only the **cell metadata** from CellxGene Census. This lightweight query lets us preview dataset sizes, available cancer types, and cell type distributions so we can make informed filtering decisions.

### What is CellxGene Census?
[CellxGene Census](https://chanzuckerberg.github.io/cellxgene-census/) provides programmatic access to a standardized, curated collection of single-cell RNA-seq datasets. All datasets use consistent cell ontology annotations, enabling cross-study comparisons.

### How to Interpret
- The **disease** table shows which cancer types are available and how many cells each has
- The **cell type** breakdown reveals the cellular composition of the TME
- Large datasets (>100k cells) may require the `CAP` parameter to cap cell counts for speed/memory"""))

cells.append(code("""import cellxgene_census
import os
import anndata as ad
output_dir = os.path.join('..', 'output')
os.makedirs(output_dir, exist_ok=True)
LOAD_FROM_CACHE = os.path.exists(h5ad_path)

if LOAD_FROM_CACHE:
    print(f"✅ Local H5AD cache found at: {h5ad_path}")
    print("   Loading cell metadata directly from cached AnnData...")
    cache_adata = ad.read_h5ad(h5ad_path)
    obs_df = cache_adata.obs.copy()
    # Ensure soma_joinid is present in the DataFrame if it's referenced downstream
    if "soma_joinid" not in obs_df.columns:
        obs_df["soma_joinid"] = obs_df.index
    print(f"\\n✅ Retrieved metadata for {len(obs_df):,} cells from cache")
else:
    # Build obs filter string
    disease_str = " or ".join(f'disease == \"{d}\"' for d in DISEASE_FILTER)
    obs_filter = f"is_primary_data == True and ({disease_str})"
    if TISSUE_FILTER:
        tissue_str = " or ".join(f'tissue_general == \"{t}\"' for t in TISSUE_FILTER)
        obs_filter += f" and ({tissue_str})"

    print(f"Census version: {CENSUS_VERSION}")
    print(f"Obs filter: {obs_filter}")

    # Query metadata only (lightweight)
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        obs_df = cellxgene_census.get_obs(
            census, ORGANISM,
            value_filter=obs_filter,
            column_names=["soma_joinid","cell_type","tissue_general",
                           "disease","dataset_id","assay"]
        )

    obs_df = obs_df.to_pandas() if hasattr(obs_df, 'to_pandas') else pd.DataFrame(obs_df)
    print(f"\\n✅ Retrieved metadata for {len(obs_df):,} cells")"""))

cells.append(code("""# Disease distribution
print("=== Disease Distribution ===")
disease_counts = obs_df['disease'].value_counts()
display(disease_counts.to_frame('cell_count'))

# Cell type distribution
print("\\n=== Top 20 Cell Types ===")
ct_counts = obs_df['cell_type'].value_counts().head(20)
fig, ax = plt.subplots(figsize=(10, 6))
ct_counts.plot.barh(ax=ax, color=sns.color_palette('viridis', len(ct_counts)))
ax.set_xlabel('Number of Cells')
ax.set_title(f'Top 20 Cell Types in Cancer Datasets ({len(obs_df):,} total cells)')
ax.invert_yaxis()
plt.tight_layout()
plt.show()

# Tissue breakdown
if 'tissue_general' in obs_df.columns:
    print("\\n=== Tissue Distribution ===")
    display(obs_df['tissue_general'].value_counts().to_frame('cell_count'))"""))

# ── Section 4: Fetch expression data ──
cells.append(md("""---
## 4. Fetch Cancer scRNA-seq Expression Data

### Purpose
Now we download the actual gene expression matrix, but **only for our target genes**. This "gene-restricted" query is critical for memory efficiency — instead of fetching ~30,000 genes per cell, we fetch only the ~200-500 genes from our MetabConnectomeDB target list.

### Methodology
- `cellxgene_census.get_anndata()` returns a standard AnnData object compatible with scanpy
- We filter to `is_primary_data == True` to avoid duplicate cells across datasets
- The `CAP` parameter provides a safety cap to downsample the dataset if specified

### How to Interpret
- The AnnData `.obs` contains cell metadata (cell type, disease, tissue)
- The `.var` contains gene metadata
- The `.X` matrix contains raw expression counts"""))

cells.append(code("""try:
    from tqdm.notebook import tqdm
except ImportError:
    from tqdm import tqdm
import anndata as ad

# 1. Subsample cell IDs FIRST to select a CONTIGUOUS range of cells (prevents remote seek lag)
obs_df_sorted = obs_df.sort_values('soma_joinid')
if CAP is not None and len(obs_df_sorted) > CAP:
    print(f"⚠️  Downsampling metadata from {len(obs_df_sorted):,} to {CAP:,} cells (contiguous slice)...")
    obs_df_sub = obs_df_sorted.head(CAP)
else:
    print(f"✅ Downloading all available cells ({len(obs_df_sorted):,})...")
    obs_df_sub = obs_df_sorted.copy()

target_cell_ids = obs_df_sub['soma_joinid'].tolist()
print(f"Ready to download {len(target_cell_ids):,} contiguous cells...")

# 2. Map target genes to SOMA feature IDs and download expression data in batches
BATCH_SIZE = 2500
adata_batches = []
var_ids = []

if LOAD_FROM_CACHE:
    print("✅ Local cache active. Extracting gene symbols and SOMA coordinates from cached AnnData...")
    # Extract var coordinates/ids from the cached AnnData directly if present
    if 'soma_joinid' in cache_adata.var.columns:
        var_ids = cache_adata.var['soma_joinid'].dropna().tolist()
    else:
        var_ids = list(range(cache_adata.n_vars))
    print(f"✅ Extracted {len(var_ids):,} genes from cached AnnData.")
else:
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        # Get the var dataframe to map target gene symbols to their SOMA coordinates
        print("Mapping target genes to CellxGene Census feature IDs...")
        var_df = cellxgene_census.get_var(
            census,
            organism=ORGANISM,
            column_names=["soma_joinid", "feature_name"]
        )
        
        # Identify which genes are unmapped in this initial pass
        census_genes_set = set(var_df["feature_name"].tolist())
        unmapped_genes = [g for g in target_genes[:500] if g not in census_genes_set]
        
        print(f"Initial mapping: Mapped {500 - len(unmapped_genes)} / 500 target genes.")
        if unmapped_genes:
            print(f"\\n⚠️ Unmapped target genes ({len(unmapped_genes)}):")
            print(", ".join(unmapped_genes))
            
            print("\\n🔍 Querying MyGene.info to find if they need the latest gene names to map...")
            try:
                import requests
                url = "https://mygene.info/v3/query"
                payload = {
                    "q": ",".join(unmapped_genes),
                    "scopes": "symbol,alias",
                    "fields": "symbol",
                    "species": "human",
                    "size": 1
                }
                response = requests.post(url, data=payload, timeout=10)
                if response.status_code == 200:
                    results = response.json()
                    symbol_mapping = {}
                    for item in results:
                        q = item.get("query")
                        approved = item.get("symbol")
                        if q and approved and approved != q:
                            symbol_mapping[q] = approved
                    
                    # Manual translations for known mitochondrial or legacy symbols
                    manual_mapping = {
                        "CYTB": "MT-CYB"
                    }
                    for old, new in manual_mapping.items():
                        if old in unmapped_genes:
                            symbol_mapping[old] = new
                    
                    if symbol_mapping:
                        print(f"Found {len(symbol_mapping)} updated / alternative approved gene symbols:")
                        resolved_count = 0
                        actual_replacements = {}
                        for old, new in symbol_mapping.items():
                            if new in census_genes_set:
                                print(f"  * '{old}' can be successfully mapped using latest symbol: '{new}'")
                                resolved_count += 1
                                actual_replacements[old] = new
                            else:
                                print(f"  * '{old}' has approved symbol '{new}' but still missing in Census")
                        
                        if actual_replacements:
                            # Remap unmapped genes in-place in global target_df and rebuild target_genes list
                            target_df['Target'] = target_df['Target'].replace(actual_replacements)
                            target_genes = sorted(target_df['Target'].dropna().unique().tolist())
                            print(f"\\nSummary: Successfully updated {resolved_count} unmapped genes to their latest HGNC approved symbols in memory!")
                    else:
                        print("No gene symbol updates found on MyGene.info.")
                else:
                    print(f"MyGene.info API returned status code {response.status_code}")
            except Exception as e:
                print(f"Could not connect to MyGene.info API: {e}")
                
        # Filter using pandas to find our target genes
        var_df_filtered = var_df[var_df["feature_name"].isin(target_genes[:500])]
        var_ids = var_df_filtered["soma_joinid"].tolist()
        print(f"\\nFinal SOMA Coord Mapping: Mapped {len(var_ids):,} / {min(500, len(target_genes)):,} target genes to SOMA coords.")
        
        if len(var_ids) == 0:
            raise ValueError("❌ None of the target genes were found in the CellxGene Census var table!")

# 3. Load from local cache if present, otherwise download cell expression data using high-speed SOMA coordinate streaming
import time
adata = None

if os.path.exists(h5ad_path):
    print(f"✅ Local AnnData cache found at: {h5ad_path}")
    print("   Loading pre-downloaded dataset from disk...")
    adata = cache_adata
    print("✅ Loaded successfully!")
else:
    max_retries = 3
    for attempt in range(1, max_retries + 1):
        try:
            print(f"Downloading cell expression data (attempt {attempt}/{max_retries})...")
            with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
                if DOWNLOAD_MODE == "target_only":
                    print(f"Downloading target-restricted matrix ({len(var_ids)} genes)...")
                    adata = cellxgene_census.get_anndata(
                        census, ORGANISM,
                        obs_coords=target_cell_ids,
                        var_coords=var_ids,
                        obs_column_names=["cell_type","tissue_general","disease",
                                          "dataset_id","assay","donor_id"],
                    )
                else:
                    print("Downloading whole transcriptome matrix (~60,000 genes)...")
                    adata = cellxgene_census.get_anndata(
                        census, ORGANISM,
                        obs_coords=target_cell_ids,
                        # var_coords=None downloads ALL genes!
                        obs_column_names=["cell_type","tissue_general","disease",
                                          "dataset_id","assay","donor_id"],
                    )
            break
        except Exception as e:
            if attempt == max_retries:
                print(f"❌ Failed to download expression data after {max_retries} attempts.")
                raise e
            print(f"⚠️ S3 connection timed out or interrupted. Retrying in 5 seconds... Error: {e}")
            time.sleep(5)

if adata is not None:
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name'].astype(str)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.var.index.name = None  # Prevent DataFrame index name conflicts during write_h5ad
    
    # Remove unused categories from categorical metadata columns to prevent downstream KeyErrors (e.g. in LIANA or plotting)
    for col in adata.obs.columns:
        if hasattr(adata.obs[col], 'cat'):
            adata.obs[col] = adata.obs[col].cat.remove_unused_categories()

print(f"\\n✅ AnnData object created: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"   Cell types: {adata.obs['cell_type'].nunique()}")
print(f"   Diseases:   {adata.obs['disease'].nunique()}")
print(adata)

# Save the raw downloaded AnnData object to disk if it was not loaded locally
if not os.path.exists(h5ad_path):
    print(f"\\nSaving downloaded AnnData to: {h5ad_path} ...")
    adata.write_h5ad(h5ad_path, compression='gzip')
    print("✅ AnnData successfully saved to disk!")"""))

# ── Section 5: QC ──
cells.append(md("""---
## 5. Dataset Overview & Quality Control

### Purpose
Assess the quality and composition of the downloaded data before downstream analysis. Identifying biases (e.g., one dataset dominating) or low-quality cells is essential for reliable interpretation.

### How to Interpret
- **Uneven cell type proportions** may reflect biological reality (e.g., tumors have many epithelial cells) or technical sampling bias
- **Low gene counts per cell** may indicate low-quality cells or ambient RNA
- **Dataset dominance** — if one dataset contributes >80% of cells, results may not generalize"""))

cells.append(code("""# Basic QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(16, 4))

# Cell type distribution
ct = adata.obs['cell_type'].value_counts().head(15)
ct.plot.barh(ax=axes[0], color=sns.color_palette('Set2', len(ct)))
axes[0].set_title('Top 15 Cell Types')
axes[0].set_xlabel('Cells')
axes[0].invert_yaxis()

# Total counts distribution
axes[1].hist(adata.obs['total_counts'], bins=50, color='steelblue', edgecolor='white')
axes[1].set_title('Total Counts per Cell')
axes[1].set_xlabel('Total Counts')

# Genes detected per cell
axes[2].hist(adata.obs['n_genes_by_counts'], bins=50, color='coral', edgecolor='white')
axes[2].set_title('Genes Detected per Cell')
axes[2].set_xlabel('N Genes')

plt.tight_layout()
plt.show()

# Disease breakdown
print("\\n=== Cells per Disease ===")
display(adata.obs['disease'].value_counts().to_frame('cells'))"""))

# ── Section 6: Expression landscape ──
cells.append(md("""---
## 6. Target Gene Expression Landscape

### Purpose
Visualize how MetabConnectomeDB target genes are expressed across cell types in the cancer TME. This reveals which cell types are "metabolically active" for specific signaling axes.

### How to Interpret
- **High expression in immune cells** (e.g., macrophages, T cells) suggests immune-metabolic crosstalk
- **Tumor-specific expression** identifies metabolite pathways the tumor uses for its own benefit
- **Ubiquitous expression** indicates housekeeping metabolic functions vs. **cell-type-specific** expression indicating specialized signaling"""))

cells.append(code("""# Normalize for visualization
adata_norm = adata.copy()
sc.pp.normalize_total(adata_norm, target_sum=1e4)
sc.pp.log1p(adata_norm)

# Mean expression per cell type (top variable genes)
sc.pp.highly_variable_genes(adata_norm, n_top_genes=min(50, adata_norm.n_vars), 
                             flavor='seurat_v3', span=1.0, subset=False)
hvg = adata_norm.var[adata_norm.var['highly_variable']].index.tolist()

if len(hvg) > 5:
    sc.pl.dotplot(adata_norm, var_names=hvg[:30], groupby='cell_type',
                  standard_scale='var', show=True)
    print(f"Dotplot: top {min(30, len(hvg))} variable target genes × cell types")
else:
    print("⚠️  Fewer than 5 variable genes found; skipping dotplot.")
    print("   This may happen with very few target genes in the dataset.")"""))

# Dimensionality Reduction & Proportions Visualizations
cells.append(md("""### Dimensionality Reduction & Tissue Microenvironment Composition
Map the high-dimensional metabolic target gene space to a 2D UMAP projection and analyze how different cell types infiltrate different cancer cohorts."""))

cells.append(code("""# 1. Run PCA & UMAP Dimensional Reduction
print("Running dimensionality reduction (PCA & UMAP)...")
sc.tl.pca(adata_norm, svd_solver='arpack')
sc.pp.neighbors(adata_norm, n_neighbors=15, n_pcs=20)
sc.tl.umap(adata_norm)

# Set premium figure styling
sc.set_figure_params(dpi=150, facecolor='white', frameon=False)
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(adata_norm, color='cell_type', ax=axes[0], show=False, 
           title='UMAP: Cell Type Composition', palette='tab20', size=25)

sc.pl.umap(adata_norm, color='disease', ax=axes[1], show=False, 
           title='UMAP: Disease Cohorts', palette='Set1', size=25)

plt.tight_layout()
plt.show()

# 2. Compute and Plot Cell-Type Proportions Across Disease Cohorts
print("\\nComputing cell type composition percentages per cohort...")
props = pd.crosstab(adata_norm.obs['disease'], adata_norm.obs['cell_type'], normalize='index') * 100

fig, ax = plt.subplots(figsize=(12, 7), dpi=150)
props.plot(kind='barh', stacked=True, ax=ax, colormap='tab20')

ax.set_title('Cell Type Microenvironment Proportions Across Cancer Cohorts', fontsize=14, fontweight='bold', pad=15)
ax.set_xlabel('Percentage of Cells (%)', fontsize=12)
ax.set_ylabel('Disease Group', fontsize=12)
ax.legend(title='Cell Type', bbox_to_anchor=(1.02, 1), loc='upper left', frameon=True, shadow=True)
plt.grid(axis='x', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()"""))

# Key Immune-Metabolic Target Expression Panel
cells.append(md("""### Expression Profile of Key Immune-Metabolic Targets
Visualize the detailed expression levels and detection rates of high-priority metabolic target genes (e.g. prostaglandin synthases, adenosine receptors, adrenergic receptors, and housekeeping/electron transport chain targets) across cell types in the cancer microenvironment."""))

cells.append(code("""# Select key metabolic target genes relevant to MetabConnectomeDB
key_metabolic_genes = [g for g in ['PTGS1', 'PTGS2', 'GRK2', 'MT-CYB', 'ADORA2A', 'ADORA2B'] if g in adata_norm.var_names]

if key_metabolic_genes:
    print(f"Plotting expression distribution of key targets: {key_metabolic_genes}")
    
    # 1. Violin Plot Panel
    sc.pl.violin(adata_norm, keys=key_metabolic_genes, groupby='cell_type', rotation=90, show=True)
    
    # 2. Dotplot mapping expression and detection rate
    sc.pl.dotplot(adata_norm, var_names=key_metabolic_genes, groupby='cell_type',
                  standard_scale='var', cmap='plasma', 
                  title='Expression & Detection of Key Metabolic Targets', show=True)
else:
    print("⚠️  No key target matches found in target genes; skipping panel.")"""))

# Genome-wide Differential Expression and Target Ranking
cells.append(md("""### Genome-wide Differential Expression & Target Ranking
If the notebook was run in **whole_transcriptome** mode, we can compute genome-wide differential expression (DE) across all cell types in the cancer microenvironment. This allows us to compare our curated MetabConnectomeDB target genes against the rest of the entire transcriptome to see where they rank and whether they serve as significant cell-type marker genes."""))

cells.append(code("""if DOWNLOAD_MODE == "whole_transcriptome":
    print("Running genome-wide differential expression (Wilcoxon rank-sum test)...")
    # Compute DE on all genes
    sc.tl.rank_genes_groups(adata_norm, groupby='cell_type', method='wilcoxon', key_added='de_results')
    
    # Extract results into a clean dataframe
    de_df = sc.get.rank_genes_groups_df(adata_norm, group=None, key='de_results')
    
    # Filter to show how MetabConnectomeDB targets rank in the transcriptome
    curated_targets = target_df['Target'].unique().tolist()
    target_de = de_df[de_df['names'].isin(curated_targets)].copy()
    
    # Display the top 15 most significant metabolic target markers
    print("\\n✅ Top 15 Highly Enriched MetabConnectomeDB Target Genes (against the whole transcriptome):")
    display(target_de.sort_values(by=['group', 'pvals_adj']).head(15))
else:
    print("⚠️  Differential expression skipped because DOWNLOAD_MODE is set to 'target_only'.")
    print("   To perform genome-wide background ranking, set DOWNLOAD_MODE = 'whole_transcriptome' in parameters.")"""))

# ── Section 7: Communication potential ──
cells.append(md("""---
## 7. Metabolite-Target Pair Cell-Type Mapping

### Purpose
Cross-reference each metabolite-target pair from MetabConnectomeDB with the cell types that actually express the target gene. This defines the **communication potential** — which cell types can "receive" each metabolite signal.

### How to Interpret
- A metabolite whose target is expressed in **macrophages** suggests tumor-to-macrophage metabolic signaling
- **Cell-type-specific targets** are high-value candidates for targeted therapy
- **Ubiquitous targets** may serve housekeeping roles and are less specific"""))

cells.append(code("""# Get gene names from adata
census_genes = set(adata_norm.var_names.tolist())

# Cross-reference with MetabConnectomeDB pairs
matched = target_df[target_df['Target'].isin(census_genes)].copy()
print(f"Matched {matched['Target'].nunique()} / {len(target_genes)} target genes in scRNA-seq data")
print(f"Covering {len(matched):,} metabolite-target pairs")

# For each matched gene, find expressing cell types (>10% detection)
# Optimization: Convert sparse matrix to CSR format and compute detection rates efficiently using tqdm
try:
    from tqdm.notebook import tqdm
except ImportError:
    from tqdm import tqdm
import scipy.sparse as sp

gene_detection = {}
unique_matched_genes = [g for g in matched['Target'].unique() if g in adata_norm.var_names]

print("Calculating gene detection rates across cell types...")
# Pre-extract the dense matrix or slice efficiently once
X_data = adata_norm.X
if sp.issparse(X_data):
    # Convert to CSR format for fast row/column slicing
    X_data = X_data.tocsr()

cell_types = adata_norm.obs['cell_type'].unique()
cell_type_masks = {ct: (adata_norm.obs['cell_type'] == ct).values for ct in cell_types}

for gene in tqdm(unique_matched_genes, desc="Mapping target genes to cell types"):
    idx = adata_norm.var_names.get_loc(gene)
    # Extract column slice efficiently
    col_slice = X_data[:, idx]
    if sp.issparse(col_slice):
        expr = col_slice.toarray().flatten()
    else:
        expr = col_slice.flatten()
        
    for ct in cell_types:
        mask = cell_type_masks[ct]
        ct_expr = expr[mask]
        det_rate = (ct_expr > 0).mean()
        if det_rate > 0.1:
            gene_detection.setdefault(gene, []).append(ct)

print(f"\\n✅ Mapping complete! Genes with >10% detection in at least one cell type: {len(gene_detection)}")

# Summary table
rows = []
for gene, cts in sorted(gene_detection.items()):
    matched_metabs = matched[matched['Target']==gene]
    metabs = []
    for _, row in matched_metabs.drop_duplicates(subset=['HMDB_ID']).iterrows():
        name = row['Metabolite_Name']
        hmdb = row['HMDB_ID']
        metabs.append(f"{name} ({hmdb})")
    rows.append({'Target_Gene': gene, 'Expressing_Cell_Types': '; '.join(sorted(cts)),
                 'N_Cell_Types': len(cts), 'Metabolites': '; '.join(sorted(metabs)[:5])})
comm_df = pd.DataFrame(rows).sort_values('N_Cell_Types', ascending=False)
display(comm_df.head(20))"""))

# ── Section 8: LIANA+ CCC ──
cells.append(md("""---
## 8. Intercellular Signaling Network Inference (LIANA+)

### Purpose
Use the **LIANA+ framework** to systematically infer cell-cell communication (CCC) events across the tumor microenvironment. LIANA+ aggregates multiple CCC methods (CellPhoneDB, NATMI, etc.) to provide robust, consensus-based interaction predictions.

### Why LIANA+ for Metabolite Signaling?
Traditional CCC methods focus on protein ligand-receptor pairs. LIANA+ extends this to metabolite-mediated communication by leveraging prior knowledge resources. This is directly relevant to our MetabConnectomeDB data, which curates metabolite-target interactions.

### How to Interpret
- **Source → Target cell type** pairs with high LIANA scores indicate active signaling axes
- Focus on interactions involving known cancer TME metabolites (e.g., Kynurenine, PGE₂, Adenosine)
- Cross-reference with Section 7 to validate that detected interactions use our curated gene targets"""))

cells.append(code("""if LIANA_AVAILABLE:
    print("Running LIANA+ cell-cell communication analysis...")
    
    # Prepare data
    adata_liana = adata_norm.copy()
    
    # Run LIANA with default methods
    try:
        liana.mt.rank_aggregate(
            adata_liana,
            groupby='cell_type',
            resource_name='consensus',
            use_raw=False,
            verbose=True,
        )
        
        # Extract results
        liana_res = adata_liana.uns['liana_res']
        print(f"\\n✅ LIANA+ found {len(liana_res):,} interactions")
        
        # Top interactions
        top = liana_res.sort_values('magnitude_rank').head(20)
        display(top[['source','target','ligand_complex','receptor_complex',
                      'magnitude_rank','specificity_rank']].reset_index(drop=True))
        
        # Filter for our target genes
        our_genes_set = set(target_genes)
        relevant = liana_res[
            liana_res['ligand_complex'].isin(our_genes_set) | 
            liana_res['receptor_complex'].isin(our_genes_set)
        ]
        print(f"\\n→ {len(relevant):,} interactions involve MetabConnectomeDB target genes")
        if len(relevant) > 0:
            display(relevant.sort_values('magnitude_rank').head(15))
            
            # ── Beautiful LIANA+ Cell-Cell Communication Heatmap ──
            print("\\n📊 Generating LIANA+ Cell-Cell Communication Heatmap...")
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
            
            # Select top 25 interactions by magnitude rank
            top_plot = relevant.sort_values('magnitude_rank').head(25)
            top_plot['Comm_Score'] = -np.log10(top_plot['magnitude_rank'] + 1e-10)
            
            pivot_df = top_plot.groupby(['source', 'target'])['Comm_Score'].max().unstack().fillna(0)
            
            plt.figure(figsize=(10, 8))
            sns.heatmap(pivot_df, annot=True, fmt=".2f", cmap="viridis", cbar_kws={'label': 'Signaling Strength (-log10 Rank)'})
            plt.title("Metabolite-Mediated Cell-Cell Communication Heatmap (LIANA+)")
            plt.ylabel("Source Cell Type (Sender)")
            plt.xlabel("Target Cell Type (Receiver)")
            plt.tight_layout()
            plt.show()
            
            # ── Beautiful standard LIANA+ Ligand-Receptor Bubble Plot ──
            print("\\n📊 Generating standard Ligand-Receptor Bubble Plot...")
            plt.figure(figsize=(12, 8))
            plot_data = top_plot.sort_values('Comm_Score', ascending=True)
            
            x_vals = plot_data['source'] + " → " + plot_data['target']
            y_vals = plot_data['ligand_complex'] + " → " + plot_data['receptor_complex']
            sizes = -np.log10(plot_data['specificity_rank'] + 1e-10) * 40 + 20
            
            bubble_sc = plt.scatter(
                x_vals, y_vals,
                s=sizes,
                c=plot_data['Comm_Score'],
                cmap='plasma',
                alpha=0.85,
                edgecolors='white',
                linewidths=0.5
            )
            cbar = plt.colorbar(bubble_sc)
            cbar.set_label('Signaling Strength (-log10 Rank)')
            plt.xticks(rotation=45, ha='right')
            plt.title("Standard Ligand-Receptor Bubble Plot (LIANA+)\\nBubble Size = Specifity Significance; Color = Signaling Strength")
            plt.xlabel("Cell-Type Communication Axis")
            plt.ylabel("Ligand → Receptor Complexes")
            plt.tight_layout()
            plt.show()
            # ── Interactive Plotly LIANA+ Visualizations ──
            try:
                import plotly.express as px
                import plotly.graph_objects as go
                PLOTLY_AVAILABLE = True
            except ImportError:
                PLOTLY_AVAILABLE = False

            if PLOTLY_AVAILABLE:
                print("\\n✨ Generating Interactive LIANA+ Plotly Visualizations...")
                # 1. Interactive Heatmap
                fig_hm = px.imshow(
                    pivot_df,
                    labels=dict(x="Target Cell Type (Receiver)", y="Source Cell Type (Sender)", color="Signaling Strength"),
                    color_continuous_scale="Viridis",
                    title="Interactive Cell-Cell Communication Heatmap (LIANA+)"
                )
                fig_hm.update_layout(width=700, height=600)
                fig_hm.show()
                
                # 2. Interactive Bubble Plot
                plot_data_plotly = top_plot.copy()
                plot_data_plotly["Axis"] = plot_data_plotly["source"] + " → " + plot_data_plotly["target"]
                plot_data_plotly["Interaction"] = plot_data_plotly["ligand_complex"] + " → " + plot_data_plotly["receptor_complex"]
                plot_data_plotly["Specificity_Significance"] = -np.log10(plot_data_plotly["specificity_rank"] + 1e-10)
                fig_bubble = px.scatter(
                    plot_data_plotly,
                    x="Axis",
                    y="Interaction",
                    size="Specificity_Significance",
                    color="Comm_Score",
                    color_continuous_scale="Plasma",
                    title="Interactive Ligand-Receptor Bubble Plot (LIANA+)",
                    labels={"Comm_Score": "Signaling Strength", "Specificity_Significance": "-log10 Specificity Rank"}
                )
                fig_bubble.update_layout(xaxis_tickangle=-45, width=900, height=500)
                fig_bubble.show()
            else:
                print("\\n💡 Tip: Install plotly (pip install plotly) to enable interactive zooming/hovering plots in Jupyter and exported HTML!")
            
            # ── Beautiful Cell-Cell Signaling Network Analysis (using NetworkX) ──
            print("\\n📊 Performing Directed Cell-Cell Signaling Network Analysis...")
            import networkx as nx
            
            G = nx.DiGraph()
            network_data = top_plot.groupby(['source', 'target'])['Comm_Score'].sum().reset_index()
            for _, r in network_data.iterrows():
                G.add_edge(r['source'], r['target'], weight=r['Comm_Score'])
                
            in_deg = G.in_degree(weight='weight')
            out_deg = G.out_degree(weight='weight')
            
            plt.figure(figsize=(12, 10))
            pos = nx.circular_layout(G)
            
            node_sizes = [(in_deg[n] + out_deg[n]) * 100 + 300 for n in G.nodes]
            node_colors = [out_deg[n] - in_deg[n] for n in G.nodes]
            
            nx.draw_networkx_nodes(
                G, pos,
                node_size=node_sizes,
                node_color=node_colors,
                cmap='coolwarm',
                edgecolors='white',
                linewidths=1.5
            )
            
            weights = [G[u][v]['weight'] for u, v in G.edges]
            max_weight = max(weights) if weights else 1
            edge_widths = [w / max_weight * 5 + 1 for w in weights]
            
            nx.draw_networkx_edges(
                G, pos,
                width=edge_widths,
                edge_color=weights,
                edge_cmap=plt.cm.magma,
                arrowsize=20,
                arrowstyle='-|>',
                connectionstyle="arc3,rad=0.15"
            )
            
            labels_pos = {k: [v[0], v[1]+0.08] for k, v in pos.items()}
            nx.draw_networkx_labels(G, labels_pos, font_size=9, font_weight='bold', font_color='#e2e8f0')
            
            plt.title("Directed Cell-Cell Signaling Network Graph (LIANA+)\\nNode Size = Total Communication Activity; Node Color = Sender Dominance (Red) vs Receiver (Blue)\\nEdge Color = Signaling Strength")
            plt.axis('off')
            plt.tight_layout()
            plt.show()
    except Exception as e:
        print(f"⚠️  LIANA+ online resource fetch failed (offline mode): {e}")
        print("   Falling back to direct Metabolite-Mediated Cell-Cell Communication mapping...")
        LIANA_AVAILABLE = False
else:
    print("⚠️  LIANA+ not available. Install with: pip install liana")
    print("   Using direct Metabolite-Mediated Cell-Cell Communication mapping fallback...")

if not LIANA_AVAILABLE:
    # ── High-Fidelity Direct Metabolic Cell-Cell Communication Fallback Model ──
    print("\\n🔄 Running direct Metabolite-Mediated Cell-Cell Communication Inference...")
    
    # Define key cancer metabolites, their biosynthetic/synthesizing enzymes, and target receptors/genes
    metab_ccc_db = [
        {
            "metabolite": "Kynurenine",
            "enzymes": ["IDO1", "IDO2", "TDO2"],
            "receptors": ["AHR"],
            "description": "Tryptophan-Kynurenine-AHR immunosuppressive axis"
        },
        {
            "metabolite": "Adenosine",
            "enzymes": ["NT5E", "ENTPD1"],
            "receptors": ["ADORA2A", "ADORA2B"],
            "description": "ATP-Adenosine immunosuppressive signaling"
        },
        {
            "metabolite": "Prostaglandin E2 (PGE2)",
            "enzymes": ["PTGS2", "PTGES"],
            "receptors": ["PTGER2", "PTGER4"],
            "description": "COX2-PGE2 inflammatory & immunosuppressive signaling"
        },
        {
            "metabolite": "Lactate",
            "enzymes": ["LDHA"],
            "receptors": ["SLC16A1", "HCAR1"],
            "description": "Glycolytic lactate tumor-microenvironment acidification"
        }
    ]
    
    # Calculate mean expression of each enzyme and receptor across cell types
    cell_types = adata_norm.obs['cell_type'].unique()
    ccc_results = []
    import scipy.sparse as sp
    
    for entry in metab_ccc_db:
        metab = entry["metabolite"]
        # Check which enzymes and receptors are in our dataset
        valid_enzymes = [e for e in entry["enzymes"] if e in adata_norm.var_names]
        valid_receptors = [r for r in entry["receptors"] if r in adata_norm.var_names]
        
        if not valid_enzymes or not valid_receptors:
            continue
            
        print(f"Analyzing {metab} signaling axis (Enzymes: {valid_enzymes} → Receptors: {valid_receptors})...")
        
        # Calculate cell-type specific mean expression
        # 1. Enzyme expression (Source capability)
        source_expr = {}
        for ct in cell_types:
            mask = (adata_norm.obs['cell_type'] == ct)
            if mask.sum() == 0:
                continue
            # Get mean expression of valid enzymes
            expr_vals = adata_norm[mask, valid_enzymes].X
            if sp.issparse(expr_vals):
                expr_vals = expr_vals.toarray()
            source_expr[ct] = float(expr_vals.mean())
            
        # 2. Receptor expression (Target capability)
        target_expr = {}
        for ct in cell_types:
            mask = (adata_norm.obs['cell_type'] == ct)
            if mask.sum() == 0:
                continue
            # Get mean expression of valid receptors
            expr_vals = adata_norm[mask, valid_receptors].X
            if sp.issparse(expr_vals):
                expr_vals = expr_vals.toarray()
            target_expr[ct] = float(expr_vals.mean())
            
        # Infer source-target communication scores
        for src in cell_types:
            for tgt in cell_types:
                score = source_expr.get(src, 0) * target_expr.get(tgt, 0)
                if score > 0:
                    ccc_results.append({
                        "Metabolite": metab,
                        "Source_CellType": src,
                        "Target_CellType": tgt,
                        "Source_Enzyme_Expression": source_expr.get(src, 0),
                        "Target_Receptor_Expression": target_expr.get(tgt, 0),
                        "Communication_Score": score,
                        "Description": entry["description"]
                    })
    
    ccc_df = pd.DataFrame(ccc_results)
    if len(ccc_df) > 0:
        print(f"\\n✅ Successfully inferred {len(ccc_df):,} metabolic cell-cell communication links!")
        display(ccc_df.sort_values('Communication_Score', ascending=False).head(15))
        
        # Let's plot and save a beautiful heatmap/clustermap for each metabolite's network!
        import matplotlib.pyplot as plt
        import seaborn as sns
        import numpy as np
        
        # 1. Plot top metabolite signaling network Heatmap
        top_metab = ccc_df.sort_values('Communication_Score', ascending=False)['Metabolite'].iloc[0]
        plot_df = ccc_df[ccc_df['Metabolite'] == top_metab]
        pivot_df = plot_df.pivot(index='Source_CellType', columns='Target_CellType', values='Communication_Score').fillna(0)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(pivot_df, annot=True, fmt=".4f", cmap="magma", cbar_kws={'label': 'Communication Score'})
        plt.title(f"Metabolite-Mediated Cell-Cell Communication Heatmap: {top_metab}\\n({plot_df['Description'].iloc[0]})")
        plt.ylabel("Source Cell Type (Produces Metabolite)")
        plt.xlabel("Target Cell Type (Expresses Receptors)")
        plt.tight_layout()
        plt.show()
        
        # 2. Standard Bubble Plot for local fallback
        print("\\n📊 Generating standard Metabolite-Receptor Bubble Plot...")
        plt.figure(figsize=(12, 8))
        
        x_vals = ccc_df['Source_CellType'] + " → " + ccc_df['Target_CellType']
        y_vals = ccc_df['Metabolite'] + " (" + ccc_df['Description'] + ")"
        sizes = ccc_df['Target_Receptor_Expression'] * 500 + 50
        
        bubble_sc = plt.scatter(
            x_vals, y_vals,
            s=sizes,
            c=ccc_df['Communication_Score'],
            cmap='plasma',
            alpha=0.85,
            edgecolors='white',
            linewidths=0.5
        )
        cbar = plt.colorbar(bubble_sc)
        cbar.set_label('Communication Score')
        plt.xticks(rotation=45, ha='right')
        plt.title("Standard Metabolite-Receptor Bubble Plot (Direct Fallback Mapping)\\nBubble Size = Target Receptor Expression; Color = Communication Score")
        plt.xlabel("Cell-Type Communication Axis")
        plt.ylabel("Metabolite-Receptor Signaling Pathways")
        plt.tight_layout()
        plt.show()
        
        # ── Interactive Plotly Fallback Visualizations ──
        try:
            import plotly.express as px
            import plotly.graph_objects as go
            PLOTLY_AVAILABLE = True
        except ImportError:
            PLOTLY_AVAILABLE = False

        if PLOTLY_AVAILABLE:
            print("\\n✨ Generating Interactive Fallback Plotly Visualizations...")
            # 1. Interactive Heatmap
            fig_hm = px.imshow(
                pivot_df,
                labels=dict(x="Target Cell Type (Expresses Receptors)", y="Source Cell Type (Produces Metabolite)", color="Communication Score"),
                color_continuous_scale="Magma",
                title=f"Interactive Cell-Cell Communication Heatmap: {top_metab}"
            )
            fig_hm.update_layout(width=700, height=600)
            fig_hm.show()
            
            # 2. Interactive Bubble Plot
            plot_data = ccc_df.copy()
            plot_data["Axis"] = plot_data["Source_CellType"] + " → " + plot_data["Target_CellType"]
            fig_bubble = px.scatter(
                plot_data,
                x="Axis",
                y="Metabolite",
                size="Target_Receptor_Expression",
                color="Communication_Score",
                hover_data=["Source_Enzyme_Expression", "Description"],
                color_continuous_scale="Viridis",
                title="Interactive Metabolite-Receptor Bubble Plot",
                labels={"Communication_Score": "Score", "Target_Receptor_Expression": "Receptor Expression"}
            )
            fig_bubble.update_layout(xaxis_tickangle=-45, width=900, height=500)
            fig_bubble.show()
        else:
            print("\\n💡 Tip: Install plotly (pip install plotly) to enable interactive zooming/hovering plots in Jupyter and exported HTML!")
        
        # 3. Directed Network Graph using NetworkX
        print("\\n📊 Performing Directed Cell-Cell Signaling Network Analysis...")
        import networkx as nx
        
        G = nx.DiGraph()
        network_data = ccc_df.groupby(['Source_CellType', 'Target_CellType'])['Communication_Score'].sum().reset_index()
        for _, r in network_data.iterrows():
            G.add_edge(r['Source_CellType'], r['Target_CellType'], weight=r['Communication_Score'])
            
        in_deg = G.in_degree(weight='weight')
        out_deg = G.out_degree(weight='weight')
        
        plt.figure(figsize=(12, 10))
        pos = nx.circular_layout(G)
        
        node_sizes = [(in_deg[n] + out_deg[n]) * 100 + 300 for n in G.nodes]
        node_colors = [out_deg[n] - in_deg[n] for n in G.nodes]
        
        nx.draw_networkx_nodes(
            G, pos,
            node_size=node_sizes,
            node_color=node_colors,
            cmap='coolwarm',
            edgecolors='white',
            linewidths=1.5
        )
        
        weights = [G[u][v]['weight'] for u, v in G.edges]
        max_weight = max(weights) if weights else 1
        edge_widths = [w / max_weight * 5 + 1 for w in weights]
        
        nx.draw_networkx_edges(
            G, pos,
            width=edge_widths,
            edge_color=weights,
            edge_cmap=plt.cm.magma,
            arrowsize=20,
            arrowstyle='-|>',
            connectionstyle="arc3,rad=0.15"
        )
        
        labels_pos = {k: [v[0], v[1]+0.08] for k, v in pos.items()}
        nx.draw_networkx_labels(G, labels_pos, font_size=9, font_weight='bold', font_color='#e2e8f0')
        
        plt.title("Directed Cell-Cell Signaling Network Graph (Direct Fallback Mapping)\\nNode Size = Total Communication Activity; Node Color = Sender Dominance (Red) vs Receiver (Blue)\\nEdge Color = Signaling Strength")
        plt.axis('off')
        plt.tight_layout()
        plt.show()
    else:
        print("⚠️ No valid metabolic communication axes could be computed from available genes.")"""))

# ── Section 9: Cancer pathways ──
cells.append(md("""---
## 9. Cancer Pathway-Level Analysis

### Purpose
Group our target genes into known **cancer metabolic signaling pathways** and examine their cell-type-specific expression. This provides a pathway-level view of metabolic communication in the TME.

### Key Cancer Metabolic Pathways
| Pathway | Key Metabolite | Key Genes | Cancer Relevance |
|---------|---------------|-----------|-----------------|
| IDO1/Kynurenine | Kynurenine, Tryptophan | IDO1, TDO2, AHR, KMO | Immune suppression via T cell exhaustion |
| xCT/Glutamate | Glutamate, Cystine | SLC7A11, SLC3A2 | Ferroptosis resistance, oxidative stress |
| CD73/Adenosine | Adenosine, ATP | NT5E (CD73), ADORA2A, ENTPD1 (CD39) | Purinergic immune checkpoint |
| COX-2/PGE₂ | PGE₂, Arachidonic acid | PTGS2 (COX-2), PTGER2, PTGER4 | Inflammation, immune evasion |
| SPHK1/S1P | S1P | SPHK1, S1PR1-5 | Tumor angiogenesis, lymphocyte trafficking |

### How to Interpret
- **High IDO1 in tumor cells** with **AHR in T cells** → tumor-mediated immune suppression
- **SLC7A11 overexpression** → ferroptosis resistance mechanism
- **NT5E + ADORA2A axis** → purinergic immune checkpoint"""))

cells.append(code("""# Define cancer metabolic pathway gene sets
CANCER_PATHWAYS = {
    'IDO1/Kynurenine': ['IDO1','IDO2','TDO2','AHR','KMO','KYNU','AFMID'],
    'xCT/Glutamate':   ['SLC7A11','SLC3A2','GLS','GLUL','SLC1A5'],
    'CD73/Adenosine':   ['NT5E','ADORA2A','ADORA2B','ENTPD1','ADA'],
    'COX-2/PGE2':       ['PTGS2','PTGS1','PTGER2','PTGER4','TBXAS1'],
    'SPHK1/S1P':        ['SPHK1','SPHK2','S1PR1','S1PR2','S1PR3','S1PR4','S1PR5'],
}

# Check which pathway genes are in our dataset
available_genes = set(adata_norm.var_names)
pathway_results = {}
for pathway, genes in CANCER_PATHWAYS.items():
    found = [g for g in genes if g in available_genes]
    pathway_results[pathway] = found
    status = "✅" if found else "❌"
    print(f"{status} {pathway}: {len(found)}/{len(genes)} genes found → {found}")

# Dotplot of pathway genes
all_pathway_genes = [g for genes in pathway_results.values() for g in genes]
if len(all_pathway_genes) >= 2:
    sc.pl.dotplot(adata_norm, var_names=all_pathway_genes, groupby='cell_type',
                  standard_scale='var', title='Cancer Metabolic Pathway Genes × Cell Types')
else:
    print("\\n⚠️  Too few pathway genes found for dotplot visualization.")"""))

# ── Section 10: Export ──
cells.append(md("""---
## 10. Export & Comprehensive Summary

### Purpose
Export key result tables and generate a text summary of findings for downstream use.

### Outputs
All exports go to `../output/` for consistency with the existing pipeline."""))

cells.append(code("""import os
output_dir = os.path.join('..', 'output')

# Export communication potential table
if len(comm_df) > 0:
    out_comm = os.path.join(output_dir, 'cellxgene_communication_potential.csv')
    comm_df.to_csv(out_comm, index=False)
    print(f"✅ Exported communication potential: {out_comm}")

# Export LIANA results if available
if LIANA_AVAILABLE and 'liana_res' in dir():
    try:
        out_liana = os.path.join(output_dir, 'cellxgene_liana_results.csv')
        liana_res.to_csv(out_liana, index=False)
        print(f"✅ Exported LIANA results: {out_liana}")
    except:
        pass

# Summary
print("\\n" + "="*60)
print("ANALYSIS SUMMARY")
print("="*60)
print(f"Organism:           {ORGANISM}")
print(f"Diseases queried:   {DISEASE_FILTER}")
print(f"Tissues:            {TISSUE_FILTER or 'all'}")
print(f"Total cells:        {adata.n_obs:,}")
print(f"Cell types:         {adata.obs['cell_type'].nunique()}")
print(f"Target genes in DB: {len(target_genes)}")
print(f"Genes found in data: {adata.n_vars}")
if len(comm_df) > 0:
    print(f"Genes with >10% detection: {len(comm_df)}")
print("="*60)
print("\\n✅ Analysis complete. See ../output/ for exported tables.")"""))

cells.append(code("""# ==========================================
# 📄 FULL NOTEBOOK HTML REPORT EXPORT
# ==========================================
# To export this notebook to a stunning HTML report, set SAVE_AS_HTML = True below:
SAVE_AS_HTML = True

if 'SAVE_AS_HTML' in globals() and SAVE_AS_HTML:
    import subprocess
    import os
    
    notebook_filename = 'cancer_cellxgene_integration.ipynb'
    output_html = '../output/cancer_cellxgene_integration_full_report.html'
    
    print(f"Executing full notebook HTML export for '{notebook_filename}'...")
    
    print("Generating gorgeous, styling-preserved HTML report...")
    cmd_html = [
        'jupyter', 'nbconvert', '--to', 'html', 
        notebook_filename, '--output', output_html
    ]
    res_html = subprocess.run(cmd_html, capture_output=True, text=True)
    
    if res_html.returncode == 0:
        print(f"🎉 SUCCESS: Notebook successfully exported as a styled HTML report!")
        print(f"   -> Saved to: '{output_html}'")
        print("\\n💡 Tip: You can open this HTML file in any web browser to view, share, or print to PDF (Cmd+P)!")
    else:
        print("❌ HTML export failed. Error details:")
        print(res_html.stderr)
else:
    print("Full notebook HTML export is currently disabled. Set SAVE_AS_HTML = True to compile the report!")"""))

# ── Build notebook ──
nb = {
    "nbformat": 4,
    "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {"display_name":"Python 3","language":"python","name":"python3"},
        "language_info": {"name":"python","version":"3.12.0"}
    },
    "cells": cells
}

out_path = os.path.join(os.path.dirname(__file__),
                        'cancer_cellxgene_integration.ipynb')
with open(out_path, 'w') as f:
    json.dump(nb, f, indent=1)
print(f"✅ Generated: {out_path}")
