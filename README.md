# metabConnectomeDB

This repository processes and merges multiple metabolite and protein databases (such as HMDB, CellPhoneDBv5, MEBOCOST, MetaLigand, MRCLinkDB, scCellFie, NeuronChat, and Cellinker2) to generate consolidated, unified human and mouse datasets.

## 📂 Directory Structure

```text
metabConnectomeDB/
├── scripts/                            # Pipeline scripts and utilities
│   ├── annotate_with_hmdb.py
│   ├── execute_and_export_notebooks.py # Exports notebooks to styled HTML reports
│   ├── generate_cellxgene_notebook.py  # Dynamically generates the CellxGene notebook
│   ├── generate_final_outputs.py
│   ├── merge_dbs_claude.py
│   ├── merge_simplify_annotate.sh
│   ├── standardize_categories.py       # Standardizes metabolite classifications
│   ├── unique_metab_data_exploration.ipynb        # Unique metabolite EDA
│   ├── metab_targetPair_analysis.ipynb            # Metabolite-Target pair analysis
│   ├── cancer_cellxgene_integration.ipynb         # CellxGene cancer scRNA-seq integration
│   └── primary_vs_metastasis_comparison.ipynb     # Primary vs metastatic DE & signaling comparison
├── input/                              # Main data directory
│   ├── databases/                      # Raw database sources
│   ├── subclass_mapping.csv            # Mapping dictionary for metabolite subclasses
│   └── superclass_mapping.csv          # Mapping dictionary for metabolite superclasses
├── output/                             # Consolidated output CSV files and HTML reports
├── environment.yml                     # Conda environment definition
├── requirements.txt                    # Pip requirements definition
└── README.md                           # This documentation
```

## 🛠 Environment Setup

To ensure reproducibility and avoid issues with dependencies (like `cellxgene-census` which requires Python < 3.13) and cloud syncing, follow these setup guidelines:

> [!IMPORTANT]
> **OneDrive Sync Warning:** Avoid creating your virtual environment (`venv`) inside the project folder if it is synced with OneDrive. Cloud storage providers intercept thousands of small Python library files, causing Jupyter Lab and scripts to hang indefinitely. Always create virtual environments in a local, unsynced folder like `~/venvs/`.

**Using venv (Recommended - Python 3.12):**

Since `cellxgene-census` requires Python 3.10–3.12, create the environment using a Python 3.12 executable (e.g., from Miniconda or Homebrew):

```bash
# Create the environment in a local unsynced folder
python3.12 -m venv ~/venvs/metabConnectomeDB

# Activate it
source ~/venvs/metabConnectomeDB/bin/activate

# Install requirements
pip install --upgrade pip
pip install -r requirements.txt
pip install jupyterlab ipykernel

# Register the kernel for Jupyter
python -m ipykernel install --user --name metabConnectomeDB --display-name "Python 3.12 (metabConnectomeDB)"
```

**Using Conda:**

```bash
conda env create -f environment.yml
conda activate metabconnectome
# Register the kernel for Jupyter
python -m ipykernel install --user --name metabConnectomeDB --display-name "Python 3.12 (metabConnectomeDB)"
```

## ⏳ Pipeline Scripts & Execution Frequency

To help manage your workflow, scripts are divided into **One-Time/Cached** (infrequently run) and **Recurring** (frequently run) categories.

**One-Time / Cached Setup Scripts (Infrequent Execution):**

- `annotate_enzyme_rhea.py`: Handles UniProt API lookups and Rhea SPARQL queries for enzyme product/substrate enrichment. Because it relies heavily on external APIs, it incrementally saves results to local JSON caches in the `input/` directory. Once the caches are built, this script finishes in seconds. You only need to manually delete the caches and re-run this completely if a massive new batch of target genes is introduced or if UniProt/Rhea data drastically updates.

**Recurring Pipeline Scripts (Frequent Execution):**
These scripts should be re-run whenever your raw input datasets change or when you tweak the database merging logic.

- `merge_simplify_annotate.sh`: The master execution wrapper that sequentially runs the data processing steps.
- `merge_dbs_claude.py`: Re-run whenever raw `.csv`/`.txt` data in `input/databases/` is added or modified.
- `generate_final_outputs.py`: Re-run whenever the dataset filtering logic changes.
- `annotate_with_hmdb.py`: Re-run if you update the `HMDB_metabolites` reference dictionary.
- `annotate_with_databases.py`: The final consolidation step. It rapidly triggers the Rhea enrichment (using the caches) and merges `OtherDB`, `Reactions`, `Rhea`, and `Guide to Pharmacology` into the final master CSV outputs. Re-run this whenever any upstream pipeline script changes the dataset.
- `execute_and_export_notebooks.py`: Run at the very end of your workflow to automatically execute the Jupyter notebooks and export the updated `.html` reports to the `output/` folder.

## 📜 Detailed Scripts Overview

The following details the relationship between scripts, their inputs, internal parameters, and expected outputs.

### `merge_simplify_annotate.sh`

- **Role:** The master execution script that runs the entire pipeline in sequence.
- **Input:** None (it executes python scripts).
- **Parameters:** None.
- **Output:** Triggers the pipeline to generate all final merged and annotated CSVs.

### `merge_dbs_claude.py`

- **Role:** The primary data ingestion and standardization script (optimized version). It crawls through all raw database subfolders, normalizes conflicting column names, standardizes data formats, deduplicates entries, and merges data separated by species (human and mouse).
- **Input:** Raw database files (CSVs, TSVs, txts) from the `input/databases/` subfolders.
- **Parameters:**
  - Dynamic `PROJECT_ROOT` resolution.
  - Strict logic mapping for species identification and database-specific normalization rules.
- **Output:**
  - `merged_human_metabolites.csv`
  - `merged_mouse_metabolites.csv`
  - (and corresponding merged protein datasets).

### `generate_final_outputs.py`

- **Role:** Refines and flattens the merged datasets. It filters out invalid records (e.g., `scCellFie_value == 0`), pairs metabolites with their targets (Receptors, Enzymes, Transporters), and calculates summary statistics.
- **Input:** `merged_{species}_metabolites.csv` generated by `merge_dbs.py`.
- **Parameters:** Dynamic `PROJECT_ROOT` resolution.
- **Output:**
  - `merged_{species}_metabolites_filtered_scCellfie_value_0.csv`
  - `{species}_database_merge_unique_metab_target_pairs.csv`
  - `{species}_database_merge_unique_metab.csv`
  - `{species}_database_merge_metab_statistics.txt`

### `annotate_with_hmdb.py`

- **Role:** Enriches the unified metabolite datasets with comprehensive HMDB annotations (such as full HMDB Names). If an HMDB ID is missing, it attempts a fuzzy name lookup.
- **Input:**
  - `{species}_database_merge_unique_metab.csv`
  - `{species}_database_merge_unique_metab_target_pairs.csv`
  - `HMDB_metabolites` (primary reference file).
- **Parameters:** A dictionary mapping lowercase metabolite names to HMDB IDs.
- **Output:**
  - `{species}_database_merge_unique_metab_with_HMDB_Info.csv`
  - `{species}_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`

### `unique_metab_data_exploration.ipynb`

- **Role:** Exploratory data analysis focused on unique metabolites across all integrated databases.
- **Input:** `output/human_database_merge_unique_metab_with_HMDB_Info.csv`
- **Functionality:** Generates summary statistics and plots for database distribution, exclusive metabolite analysis, and HMDB annotation coverage.

### `metab_targetPair_analysis.ipynb`

- **Role:** Analysis of metabolite-target (Receptor/Enzyme/Transporter) pairs and their annotations.
- **Input:** `output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`
- **Functionality:** Investigates the interactions and provides a provenance breakdown of annotations like Disease and Cell Type (primarily from MRCLinkDB).

### `cancer_cellxgene_integration.ipynb`

- **Role:** Bridges MetabConnectomeDB metabolite-target gene pairs with cancer scRNA-seq data from [CellxGene Census](https://cellxgene.cziscience.com/datasets).
- **Input:** `output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`
- **Dependencies:** `cellxgene-census`, `scanpy`, `anndata`, `liana` (install via `pip install -r requirements.txt`)
- **Functionality:**
  - Parameterized filtering by disease, tissue, and organism (`DISEASE_FILTER`, `TISSUE_FILTER`, `ORGANISM`)
  - CellxGene Census metadata exploration & scRNA-seq data retrieval
  - Cell-type-resolved target gene expression analysis
  - Metabolite-target communication potential mapping
  - Intercellular signaling network inference using **LIANA+** (ligand-receptor analysis framework)
  - Cancer pathway-level analysis (IDO1/Kynurenine, xCT/Glutamate, CD73/Adenosine, etc.)
  - All analysis sections include interpretive markdown cells linking findings to cancer biology

### `primary_vs_metastasis_comparison.ipynb`

- **Role:** Directly compares metabolic signaling targets between primary tumors and metastatic sites across different cancer types.
- **Input:**
  - `output/{PRIMARY_PREFIX}.h5ad`
  - `output/{META_PREFIX}.h5ad`
  - `output/{PRIMARY_PREFIX}_cellxgene_liana_results.csv`
  - `output/{META_PREFIX}_cellxgene_liana_results.csv`
  - `output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`
- **Dependencies:** Requires that `cancer_cellxgene_integration.ipynb` is run first for both primary and metastatic sites to generate the `.h5ad` and specific LIANA+ CSV outputs.
- **Functionality:**
  - Merges primary and metastatic scRNA-seq datasets.
  - Performs Differential Expression (DE) analysis specifically on the unified `MetabConnectomeDB` vocabulary.
  - Generates interactive Plotly and static Seaborn Volcano Plots mapping Log2 Fold Change against Adjusted P-Values.
  - Integrates LIANA+ output to tag whether enriched targets are actively mediating cell-cell communication (CCC) specifically in primary, metastatic, or both sites.
  - Dynamically exports a cancer-specific final table (e.g., `primary_vs_metastasis_{cancer_type}_DE_metabolic_targets.csv`) detailing targeting metabolites and HMDB IDs.

### `merge_dbs.py`, `test_merge.py`, `test_script.py`

- **Role:** Legacy scripts or scratchpads used for testing specific functionality.

## 📊 Data Provenance & Metadata

A critical aspect of this unified database is the tracking of metadata provenance.

- **Clinical Metadata:** Columns such as `Disease`, `Cell type`, `Effect`, and `Interaction` originate solely from the **MRCLinkDB** dataset (specifically `Metabolite-cell interaction.txt`).
- **Database Overlap:** The pipeline identifies metabolites shared across multiple databases (e.g., HMDB, Celllinker2, MRCLinkDB) and merges them into a single entry with combined provenance labels. This results in 100% overlap between certain sources like Celllinker2 and MRCLinkDB in unique metabolite mappings.
- **Target Pairing:** Metabolites are paired with targets based on curated mapping files within the `input/databases/` directory.

## 🚀 Version Control (Git)

The repository is configured to exclude large data directories:

- `input/`: Contains raw database files (CSVs tracked explicitly).
- `output/`: Contains generated CSV results (ignored) and rendered HTML reports (tracked).
- `venv/`, `__pycache__/`, `.ipynb_checkpoints/`: Local environment artifacts.

Ensure that the `input/` directory is populated locally before running the pipeline.
