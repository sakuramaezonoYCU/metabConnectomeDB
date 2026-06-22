# metabConnectomeDB

This repository processes and merges multiple metabolite and protein databases (such as HMDB, CellPhoneDBv5, MEBOCOST, MetaLigand, MRCLinkDB, scCellFie, NeuronChat, and Cellinker2) to generate consolidated, unified human and mouse datasets called metabConnectomeDB.

Beyond data consolidation, this repository contains a comprehensive suite of computational pipelines and analytical notebooks for multi-omics cancer research. Key capabilities include:

- **Pan-Cancer Meta-Analysis:** Cross-cancer identification of conserved metabolic gene signatures (e.g., the 21-gene Directed Metastatic Signature and the 12-gene STAT3 Core Axis) using large-scale single-cell RNA-seq integration (CellxGene).
- **Prognostic and Predictive Modeling:** Machine learning classifiers (Random Forest, Cox Proportional Hazards, MLP Neural Networks) trained on clinical cohorts (e.g., METABRIC, TCGA) to evaluate the prognostic power of metabolic signatures.
- **Spatial Transcriptomics & Metabolomics Integration:** Validation of metabolic interactions (like the Serotonin-TAM immune evasion axis) using physical mass-spectrometry metabolomics and high-resolution spatial transcriptomics (Visium).
- **In Silico Microenvironment Analysis:** Computationally simulating intratumoral oxygen gradients, mapping macrophage immunometabolism (CAMP Integration), and modeling directional cell-cell communication networks.
- **Therapeutic Target Discovery:** Cross-referencing identified metabolic targets against pharmacological databases (DGIdb, Open Targets) to profile druggability and guide synthetic lethality strategies.

## 📂 Directory Structure

```text
metabConnectomeDB/
├── scripts/                            # Pipeline scripts and utilities
│   ├── merge_dbs_claude.py             # Primary data ingestion and standardization
│   ├── annotate_with_hmdb.py
│   ├── annotate_with_databases.py      # Final consolidation step mapping external databases
│   ├── execute_and_export_notebooks.py # Exports notebooks to styled HTML reports
│   ├── run_cancer_pipeline.py          # Executes full single-cancer pipeline
│   ├── run_all_cancers.py              # Queries CellxGene Census and runs multi-cancer pipeline
│   ├── run_validation_phase.py         # Automates downstream signature validations
│   ├── massspec_metabolomics_analysis.py # Verifies signatures via mass-spectrometry
│   ├── validate_tcga_signature.py      # TCGA Cox Proportional Hazard regressions
│   ├── verify_spatial.py               # Visium spatial enrichment scoring
│   ├── generate_predictive_notebook.py # Pre-metastatic subclone scoring
│   ├── execute_pancancer_notebooks.py  # Downstream HTML report compiler
│   ├── generate_final_outputs.py
│   ├── merge_simplify_annotate.sh
│   ├── parse_md_tables.py              # Utility to extract Markdown tables
│   ├── standardize_categories.py       # Standardizes metabolite classifications
│   ├── query_advanced_analysis.py      # Queries external APIs for gene interactions/druggability
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
│   ├── AI_summary_and_insights.md      # Comprehensive pipeline insights and next steps
│   ├── *_database_merge_unique_metab*.csv  # Final merged, target-paired, and HMDB-annotated datasets for human and mouse
│   ├── *_database_merge_metab_statistics.txt # Summary metrics and dataset distributions
│   ├── [cancer]_results/               # Cancer-specific DE results, CellxGene integration, and orphan candidates
│   ├── pan_cancer_meta_results/        # Pan-cancer intersection results, signatures, networks, and UpSet plots
│   └── master_regulator_results/       # TF enrichment results and barplots for master regulator analysis
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

## ⏳ Pipeline Execution Architecture & Phases

To guarantee 100% reproducible execution, the pipeline is strictly divided into logical phases.

**Phase 0: Core Metabolite Database Curation (Infrequent Execution)**
This phase curates the generic databases (KEGG, HMDB, CellChat) from `input/databases/`. Re-run only when your raw input datasets change or when you tweak the database merging logic.
- `merge_simplify_annotate.sh`: The master execution wrapper that sequentially runs the data processing steps below.
- `fetch_kegg_pathways.py`: Dynamically pulls pathway gene sets from the KEGG REST API based on `pipeline.config.json` specifications.
- `merge_dbs_claude.py`: Ingests, normalizes, and merges raw `.csv`/`.txt` data.
- `generate_final_outputs.py`: Filters invalid records and pairs metabolites with targets.
- `annotate_with_hmdb.py`: Enriches datasets using the `HMDB_metabolites` reference.
- `annotate_with_databases.py`: Final consolidation step utilizing cached Rhea enrichments.

**Phase 1: Database Exploration and Reporting**
- `execute_and_export_notebooks.py`: Provides a general exploratory report of the fully integrated Phase 0 database.

**Phase 2: Single-Cell Transcriptome Integration (Per-Cancer)**
- `run_all_cancers.py`: Queries CellxGene Census online database dynamically to identify the top metastatic sites and runs parallelized single-cancer pipelines.
- `run_cancer_pipeline.py`: Executes the complete single-cancer pipeline (DE, LIANA+, Immune Evasion) for specified cancers, executing notebooks and saving reports into `output/[cancer]_results/`.

**Phase 3: Meta-Analysis & Cross-Cancer Intersections**
- `extract_dataset_metrics.py`: Extracts cell counts across microenvironments.
- `count_met_cells.py`: Identifies uniquely targeted orphaned metabolites.
- `compute_pan_cancer_meta.py`: Computes the target intersection across multiple cancers to identify highly conserved combinations.
- `generate_ai_summary_tables.py`: Generates the consolidated summary CSV files.

**Phase 4: Pan-Cancer Notebook Compilation**
- `generate_combined_pan_cancer_notebook.py`: Auto-generates the `pan_cancer_meta_analysis.ipynb` master notebook.
- `execute_and_export_notebooks.py` (or manual execution): Runs the notebook to visualize the network graph and output standard plots.

**Phase 5: Dynamic Gene Signature Validation (Frequent Execution)**
This phase is orchestrator-led and validates all dynamically identified (N-1)-cancer combinations generated in Phase 4.
- `run_validation_phase.py`: The master execution script that automatically identifies all output combination signatures and routes them through the following tests:
  - `massspec_metabolomics_analysis.py`: Verifies signature genes using mass-spectrometry clinical cohorts.
  - `validate_tcga_signature.py`: Runs Cox Proportional Hazard regressions on TCGA survival datasets.
  - `verify_spatial.py`: Applies spatial enrichment scoring and calculates Moran's I on high-resolution Visium slides.
  - `generate_predictive_notebook.py`: Scores primary tumor cells directly to identify left/right skewed pre-metastatic subclones.
  - `generate_ml_prognostic_classifier_notebook.py`: Generates the Per-Cancer and Pan-Cancer ML Prognostic Classifier notebooks.
  - `create_camp_notebook.py`: Generates the CAMP Pan-Cancer metabolomics integration notebook.
  - `generate_master_regulator_notebook.py`: Generates the Master Regulator TF analysis notebook.

**Phase 6: Advanced Downstream Validation & HTML Export**
- `execute_pancancer_notebooks.py`: The master driver for the downstream notebook compilation. Executes the following into final HTML reports required by the AI generation phase:
  - Pan-Cancer Meta Analysis
  - Predictive Signature Biomarker
  - Druggability Axis Analysis
  - Visium Spatial Validation
  - Deep-Dive Conserved Metab Gene Sig
  - Serotonin Axis Spatial Mapping
  - Ovarian Serotonin Immune Evasion
  - Oxygen Tension Analysis
  - MITF Regulon Expansion
  - Per-Cancer and Pan-Cancer ML Prognostic Classifiers
  - CAMP Pan-Cancer Integration
  - Master Regulator Analysis
  
**Phase 7: AI Reporting**
- `tmp_build_md.py`: Scrapes the results of Phases 1 through 6 and leverages the Gemini API to produce the final `AI_summary_and_insights.md` pipeline document.

## 📜 Detailed Scripts Overview

The following details the relationship between scripts, their inputs, internal parameters, and expected outputs.

### `merge_simplify_annotate.sh`

- **Role:** The master execution script that runs the entire pipeline in sequence.
- **Input:** None (it executes python scripts).
- **Parameters:** None.
- **Output:** Triggers the pipeline to generate all final merged and annotated CSVs.

### `fetch_kegg_pathways.py`

- **Role:** Dynamically pulls core metabolic pathway gene sets (e.g. Glycolysis, OXPHOS, HIF-1) directly from the KEGG REST API to ensure 100% programmatic provenance and prevent hardcoding.
- **Input:** `input/pipeline.config.json` (specifically the `KEGG_PATHWAYS` block).
- **Output:** `input/kegg_{id}_{name}.json` files containing lists of corresponding HGNC gene symbols.

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

### `extract_dataset_metrics.py`

- **Role:** Extracts exactly how many cells were evaluated for the TME and Malignant compartments across primary and metastatic sites.
- **Input:** Single-Cell `.h5ad` objects.
- **Output:** `output/pan_cancer_meta_results/cell_type_counts_{suffix}.csv`

### `count_met_cells.py`

- **Role:** Identifies the total and unique metabolic target candidates implicated in cell-cell communication (CCC) potential across different tissue microenvironments.
- **Output:** Saves results indicating the spread of immune evasion orphan targets.

### `generate_ai_summary_tables.py`

- **Role:** Aggregates all the extracted outputs into formatted CSVs in `output/ai_summary_tables/` for the final AI insights generation.
- **Input:** Upstream results from differential expression and orphan receptor mapping.
- **Output:** `output/ai_summary_tables/*_summary.csv`

### `tmp_build_md.py`

- **Role:** Scrapes the exact numbers directly from the generated HTML reports to dynamically build the `AI_summary_and_insights.md`.
- **Output:** The final `output/AI_summary_and_insights.md` guaranteeing zero hallucinated metrics.

### `pan_cancer_meta_analysis.ipynb` & Generators

- **Role:** The capstone multi-cancer meta-analysis. Computes the mathematical intersection of DE results across all N configured cancers to derive a strictly conserved metastatic metabolic signature. Includes Network visualizations, druggability scoring, and predictive biomarker scoring from primary datasets.
- **Key Scripts (Run in this exact order if you change cell counts or parameters in `pan_cancer_config.py`):**
  1. `run_all_cancers.py`: Re-runs the base pipeline for all configured cancers (generates individual cancer DE results).
  2. `compute_pan_cancer_meta.py`: Intersects the differentially expressed targets across the cancers, outputting `pan_cancer_signature_XXX.csv` combinations using the dynamic >3 intersection algorithm.
  3. `generate_predictive_notebook_all5.py`: Scores the primary tumor cells of each cancer against their **OWN specific** metastatic signature to identify highly metastatic ("Right-Skewed") subclones present before dissemination.
  4. `generate_combined_pan_cancer_notebook.py`: Auto-generates the final `pan_cancer_meta_analysis.ipynb` notebook logic.
  5. Finally, execute `pan_cancer_meta_analysis.ipynb` top-to-bottom to render the final HTML report.
- **Output:** Master pan-cancer meta-analysis notebook detailing the conserved genes, network graph, druggability analysis, and predictive scores, alongside all underlying CSV data files in `output/pan_cancer_meta_results/`.

### `query_advanced_analysis.py`

- **Role:** A utility module used by downstream notebooks to query external databases for advanced target analysis. It enforces strict data integrity by immediately raising `RuntimeError` exceptions rather than returning empty placeholders if an API call to STRING or Open Targets fails.
- **Input:** Takes lists of target gene symbols (e.g., `["STAT3", "HTR7"]`) as arguments to its internal functions (`query_string_ppi`, `query_tractability`).
- **Output:** Returns Pandas DataFrames containing physical interaction scores (from STRING) and tractability assessments (Small Molecule / Antibody druggability from Open Targets).

### `run_validation_phase.py`

- **Role:** Automates the signature validation pipeline (Phase 5/6) across discovered pan-cancer signatures.
- **Input:** Identifies all `pan_cancer_conserved_genes*.csv` and `pan_cancer_signature_*.csv` files within `output/pan_cancer_meta_results/`.
- **Output:** Sequentially routes the signatures to MassSpec, TCGA, and Spatial validation scripts, and subsequently triggers predictive biomarker notebook generation.

### `massspec_metabolomics_analysis.py` & `massspec_cross_cohort_comparison.py`

- **Role:** Validates identified signature genes using mass-spectrometry clinical cohorts and performs cross-cohort comparisons.
- **Input:** `--signature_csv` for analysis, `--signature-name` for cross-cohort comparison.
- **Output:** Intermediate analysis results in `output/pan_cancer_meta_results/` and generation of the final `massspec_metabolomics_integration_*.ipynb` via `generate_massspec_metabolomics_notebook.py`.

### `validate_tcga_signature.py`

- **Role:** Evaluates signature prognostic power by running Cox Proportional Hazard regressions on TCGA survival datasets.
- **Input:** `--signature_csv` (path to a signature file).
- **Output:** TCGA validation metrics and regression summaries.

### `verify_spatial.py`

- **Role:** Applies spatial enrichment scoring and calculates Moran's I on high-resolution Visium slides to validate signatures.
- **Input:** `--signature_csv` (path to a signature file).
- **Output:** Spatial validation metrics mapping intra-tumor heterogeneity.

### `generate_predictive_notebook.py`

- **Role:** Dynamically builds a notebook to score primary tumor cells against their specific metastatic signature, finding pre-metastatic subclones.
- **Input:** Automatically detects signatures.
- **Output:** Generates `predictive_signature_biomarker.ipynb`.

### `execute_pancancer_notebooks.py`

- **Role:** Downstream master script that builds all final HTML reports from the spatial, pan-cancer, and druggability notebooks.
- **Input:** None (runs based on outputs present in the environment).
- **Output:** Generates the final styled HTML reports for AI ingestion.

### `pan_cancer_config.py`

- **Role:** Centralized configuration script to maintain data provenance and prevent hardcoding constants or API endpoints in analysis scripts.
- **Input:** None.
- **Parameters:** Defines constants such as `OPENTARGETS_API_URL`.
- **Output:** Exposes constants for other modules to import.

### Specialized Investigation Notebooks

Several targeted Jupyter notebooks dive deep into specific biological questions raised by the pan-cancer analysis:

- **`druggability_axis_analysis.ipynb`**: Investigates the clinical actionability of the highly conserved glutamine-sphingolipid-ketone body axis.
- **`master_regulator_analysis.ipynb`**: Explores whether dynamically discovered upregulated genes act as the master transcriptional switch for the pan-cancer metastatic signature using ChEA/ENCODE enrichment.
- **`deepdive_conserved_metabGeneSig.ipynb`**: Deep dive into the conserved pan-cancer metabolic gene signature, integrating STAT3 targets, directional cell-cell communication, and validating signatures in TCGA survival cohorts using permutation null distributions.
- **`mitf_regulon_expansion.ipynb`**: Investigates the expansion of the MITF regulon and its downstream metabolic targets.
- **`predictive_signature_biomarker.ipynb`**: Explores the pan-cancer predictive capability of the strictly conserved metabolic gene signature.
- **`oxygen_tension_analysis.ipynb`**: Computationally simulates and models the impact of intratumoral oxygen gradients (e.g., via HIF-1 pathway profiling) on metabolic signatures.
- **`serotonin_axis_spatial_mapping.ipynb`**: Maps the spatial distribution of the serotonin axis within specific tissue microenvironments.
- **`visium_spatial_validation.ipynb`**: Validates the spatial axis involving HTR7+ tumor-associated macrophages (TAMs) and HR-repair genes using Visium spatial transcriptomics via spatial co-localization analysis.
- **`camp_pancancer_integration_*.ipynb`**: Investigates pan-cancer integration for Directed Metastatic Signatures.
- **`massspec_metabolomics_integration_*.ipynb`**: Integrates mass spectrometry metabolomics data with cross-cohort comparisons to validate metabolic signatures.
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

## 📝 Recent Changelog

**[2026-06-20] Dynamic Pan-Cancer Unhardcoding Patch**
- **Issue:** The AI agents had carelessly hardcoded `CANCERS = ['breast', 'colorectal', 'lung', 'melanoma', 'ovarian']` directly into multiple python scripts, which bypassed the unified `pan_cancer_config.py` and completely ignored the 6th cancer (`kidney`) present in `pipeline.config.json`.
- **Fix:** Stripped all hardcoded arrays. `compute_pan_cancer_meta.py`, `generate_pan_cancer_notebook.py`, and `tmp_build_md.py` now explicitly import `CANCERS_TO_RUN` from `pan_cancer_config.py`.
- **Methodology Preservation:** If the N-cancer strict intersection yields 0 conserved genes, the pipeline now dynamically falls back to the (N-1)-cancer combinations to derive the signature, avoiding data falsification while preserving analytical viability.
- **Documentation:** Updated all `README.md` and `pipeline_execution_checklist.md` references from "5-cancer" and "4-cancer" to the algebraic "N-cancer" and "(N-1)-cancer" to reflect the dynamically scalable architecture.

**[2026-06-22] Dynamic KEGG Pathways & Plot Aesthetics Patch**
- **Issue:** Pathway gene lists for immune-metabolic targets (including TCA/Fumarate) were hardcoded in notebooks, causing missing data mappings. `run_all_cancers.py` threw a `KeyError` due to an `"ovary"` vs `"ovarian"` mismatch in `pipeline.config.json`. Plot axis labels (LIANA+ heatmaps and bubble plots) overlapped severely when processing large numbers of interacting cell types.
- **Fix:** Switched to dynamic programmatic fetching via `fetch_kegg_pathways.py` and `CANCER_PATHWAYS` json config mapping. Changed violin plots to `stacked_violin` and dynamically scaled Plotly/Matplotlib figsize widths against cell-type array lengths. Fixed `OUTPUT_DIR` injection redundancy in `run_cancer_pipeline.py`.
