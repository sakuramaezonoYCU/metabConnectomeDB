# metabConnectomeDB

This repository processes and merges multiple metabolite and protein databases (such as HMDB, CellPhoneDBv5, MEBOCOST, MetaLigand, MRCLinkDB, scCellFie, NeuronChat, and Cellinker2) to generate consolidated, unified human and mouse datasets.

## 📂 Directory Structure

```text
metabConnectomeDB/
├── scripts/                            # Pipeline scripts and utilities
│   ├── annotate_with_hmdb.py
│   ├── execute_and_export_notebooks.py # Exports notebooks to styled HTML reports
│   ├── run_cancer_pipeline.py          # Executes full single-cancer pipeline
│   ├── run_all_cancers.py              # Queries CellxGene Census and runs multi-cancer pipeline
│   ├── generate_final_outputs.py
│   ├── merge_simplify_annotate.sh
│   ├── parse_md_tables.py              # Utility to extract Markdown tables
│   ├── standardize_categories.py       # Standardizes metabolite classifications
│   ├── cancer_cellxgene_integration.ipynb         # CellxGene cancer scRNA-seq integration
│   ├── metab_targetPair_analysis.ipynb            # Metabolite-Target pair analysis
│   ├── orphan_metabolic_immune_evasion.ipynb      # Orphan metabolite immune evasion mapping
│   ├── primary_vs_metastasis_comparison.ipynb     # Primary vs metastatic DE & signaling comparison
│   └── unique_metab_data_exploration.ipynb        # Unique metabolite EDA
├── input/                              # Main data directory
│   ├── databases/                      # Raw database sources
│   ├── subclass_mapping.csv            # Mapping dictionary for metabolite subclasses
│   └── superclass_mapping.csv          # Mapping dictionary for metabolite superclasses
├── output/                             # Consolidated output CSV files and HTML reports
│   └── AI_summary_and_insights.md      # Comprehensive pipeline insights and next steps
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
- `run_cancer_pipeline.py`: Runs the complete single-cancer pipeline for specified disease filters, executing notebooks and saving reports.
- `run_all_cancers.py`: Queries CellxGene Census to identify top metastatic sites and runs the entire automated pipeline for all configured cancers.
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
- **Input:** `merged_{species}_metabolites.csv` generated by `merge_dbs_claude.py`.
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

### `orphan_metabolic_immune_evasion.ipynb`

- **Role:** Systematically cross-references literature-sparse ("Tier 2/3") metabolic target pairs against immune cell populations (B-cells, Macrophages, Dendritic cells, T-cells) to identify uncharacterized metabolic ligands mediating tumor immune evasion.
- **Input:**
  - `output/human_database_merge_unique_metab_target_pairs_with_HMDB_Info.csv`
  - Cancer-specific differential expression (DE) metabolic target tables or direct scRNA-seq AnnData object.
- **Functionality:**
  - Filters for "Orphan" (Tier 2/3) metabolite-receptor interactions.
  - Maps target gene enrichment in immune cells.
  - Visualizes top orphan candidates with a high-fidelity continuous color gradient (viridis) indicating log2 fold changes.
  - Generates `immune_evasion_orphan_metabolic_candidates.csv` in the cancer subfolder for therapeutic target discovery.

### `run_cancer_pipeline.py`

- **Role:** Automates the complete single-cancer analysis workflow by executing `cancer_cellxgene_integration.ipynb`, `primary_vs_metastasis_comparison.ipynb`, and `orphan_metabolic_immune_evasion.ipynb` in sequence. It dynamically injects filter and output parameters, routes outputs to a cancer-specific directory under `output/<cancer_type>_results/`, and handles automatic notebook execution and HTML report generation.
- **Input:**
  - Disease filters (e.g., `"breast cancer"`)
  - Target tissue general classes (e.g., `"breast,mammary gland,liver,bone,brain"`)
  - Primary tissue classes (e.g., `"breast,mammary gland"`)
- **Parameters:** Injected dynamically into notebook execution environments at runtime, bypassing manual configuration hurdles.
- **Output:**
  - Cancer-specific `.h5ad` AnnData object.
  - Interactive HTML reports for CellxGene integration, Primary vs. Metastasis comparison, and Orphan Immune Evasion.
  - Combined signaling target tables and candidate CSV files.

### `run_all_cancers.py`

- **Role:** The high-level multi-cancer pipeline orchestrator. It queries the CellxGene Census online database dynamically to identify the top 3 metastatic tissues by cell count for each configured cancer type, then launches parallelized headless runs of `run_cancer_pipeline.py` for each cancer.
- **Input:** Online CellxGene Census database.
- **Parameters:** `CELLXGENE_CAP` (set via environment variable, defaults to `100000` cells to ensure robust statistical power).
- **Output:** Fully automated, end-to-end downstream results and HTML reports for all major cancers (Breast, Colorectal, Melanoma, Lung, Ovarian).

### `pan_cancer_meta_analysis.ipynb` & Generators

- **Role:** The capstone multi-cancer meta-analysis. Computes the mathematical intersection of DE results across all 5 cancers to derive a strictly conserved metastatic metabolic signature. Includes Network visualizations, druggability scoring, and predictive biomarker scoring from primary datasets.
- **Key Scripts (Run in this exact order if you change cell counts or parameters in `pan_cancer_config.py`):**
  1. `run_all_cancers.py`: Re-runs the base pipeline for all 5 cancers (generates individual cancer DE results).
  2. `compute_pan_cancer_meta.py`: Aggregates the 5 single-cancer DE results to generate the strictly conserved gene list, UpSet plot, and network edges.
  3. `generate_predictive_notebook_all5.py`: Computes the "Metastatic Metabolic Score" across primary vs. metastatic single cells for all 5 cancers, outputting CSVs and PNGs.
  4. `generate_combined_pan_cancer_notebook.py`: Auto-generates the final `pan_cancer_meta_analysis.ipynb` notebook logic.
  5. Finally, execute `pan_cancer_meta_analysis.ipynb` top-to-bottom (or use `execute_and_export_notebooks.py`) to render the final HTML report.
- **Output:** Master 5-cancer meta-analysis notebook detailing the conserved genes, network graph, druggability analysis, and predictive scores, alongside all underlying CSV data files in `output/pan_cancer_meta_results/`.

### Specialized Investigation Notebooks

Several targeted Jupyter notebooks dive deep into specific biological questions raised by the pan-cancer analysis:

- **`druggability_axis_analysis.ipynb`**: Investigates the clinical actionability of the highly conserved glutamine-sphingolipid-ketone body axis.
- **`oxygen_tension_analysis.ipynb`**: Correlates the magnitude of metabolic shifts against the physical oxygen tension of varying metastatic niches (e.g., hypoxic pleural effusions vs. oxygenated brain).
- **`nr1d2_master_regulator_analysis.ipynb`**: Explores whether the universally upregulated gene NR1D2 acts as the master transcriptional switch for the pan-cancer metastatic signature using ChEA/ENCODE enrichment.
- **`ovarian_serotonin_immune_evasion.ipynb`**: Specifically investigates the role of up-regulated tumor-derived serotonin receptors in the suppression of local T-cells within the ovarian peritoneal metastatic niche.
- **`deepdive_conserved_metabGeneSig.ipynb`**: Deep dive into the conserved pan-cancer metabolic gene signature, integrating STAT3 targets, directional cell-cell communication, and validating signatures in TCGA survival cohorts using permutation null distributions.
- **`mitf_regulon_expansion.ipynb`**: Investigates the expansion of the MITF regulon and its downstream metabolic targets.
- **`predictive_signature_biomarker.ipynb`**: Explores the pan-cancer predictive capability of the strictly conserved metabolic gene signature.
- **`serotonin_axis_spatial_mapping.ipynb`**: Maps the spatial distribution of the serotonin axis within specific tissue microenvironments.
- **`ml_prognostic_classifier.ipynb`**: Trains Cox Proportional Hazards, Random Forest, and MLP Neural Network classifiers on independent clinical cohorts (e.g., METABRIC breast cancer dataset) to evaluate the prognostic power of our derived metabolic gene signatures. Generates risk stratification models, ROC curves, and Kaplan-Meier plots.

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
