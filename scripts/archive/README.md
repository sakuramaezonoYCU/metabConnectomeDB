# Archived Scripts

This folder contains ad-hoc, legacy, or standalone utility scripts that are **not** part of the automated `metabConnectomeDB` execution pipeline. They have been archived to keep the main `scripts/` directory clean.

## Categories

### 1. Patching & Hotfix Scripts
These scripts were used to perform one-off string replacements, patch Jupyter notebooks, fix markdown headers, or modify specific code cells without needing to re-run the entire pipeline.
- `add_citations.py`
- `add_hif1_plot.py`
- `add_inputs_to_notebooks.py`
- `fix_all_notebooks.py`
- `fix_columns.py`
- `fix_deepdive_code_cells.py`
- `fix_deepdive_paths.py`
- `fix_header_references.py`
- `fix_header_references_again.py`
- `fix_header_references_rest.py`
- `fix_headers.py`
- `fix_notebook_exports.py`
- `fix_notebook_plots.py`
- `fix_novelty_pmids.py`
- `fix_nr1d2.py`
- `fix_pan_cancer.py`
- `fix_predictive.py`
- `patch_aux_scripts.py`
- `patch_checklist.py`
- `patch_comm_modes.py`
- `patch_deepdive.py`
- `patch_druggability.py`
- `patch_druggability_181.py`
- `patch_druggability_remaining.py`
- `patch_generators.py`
- `patch_hardcoded.py`
- `patch_markdown.py`
- `patch_metab_target_analysis.py`
- `patch_notebook_code.py`
- `patch_notebook_code2.py`
- `patch_notebook_code3.py`
- `patch_nr1d2.py`
- `patch_oxygen_notebook.py`
- `patch_suffix.py`
- `remove_execute.py`
- `rename_outputs.py`
- `scrub_23.py`
- `update_all_tables_in_md.py`
- `update_druggability_notebook.py`
- `update_ipynb.py`
- `update_md.py`

### 2. Standalone Analytical Computations
Scripts that compute specific exploratory metrics or analyses that are either experimental or superseded by the main pipeline functions.
- `compute_metabolic_switching.py`
- `compute_nr1d2_enrichment.py`
- `compute_serotonin_proximity.py`
- `count_met_cells.py`

### 3. External API Queries & Downloads
Scripts previously used to download reference data or query external APIs (like OpenTargets, ChEA, or TCGA). The current pipeline likely incorporates these natively or the data is already fetched.
- `download_chea_data.py`
- `download_tcga_data.py`
- `fetch_opentargets.py`
- `fetch_uniprot_roles.py`
- `query_advanced_analysis.py`
- `query_dbs.py`
- `query_depmap.py`

### 4. Custom Notebook Generators & Exporters
Various scripts that generated specific notebook variants (e.g., specific combinations of cancers or targeted outputs) or helped dump notebook cells to text.
- `dump_nb.py`
- `export_annotated_23_genes.py`
- `export_notebook.py`
- `generate_column_definitions.py`
- `generate_nb3.py`
- `generate_ovarian_serotonin.py`
- `generate_pan_cancer_notebook.py`
- `generate_predictive_notebook_all5.py`
- `generate_version_1_metrics.py`

### 5. Execution Wrappers & Tests
Ad-hoc test scripts for verifying APIs, network state, data frames, or running specific slices of the pipeline independently.
- `inspect_ovarian.py`
- `merge_dbs.py` (Legacy version, replaced by `merge_dbs_claude.py`)
- `parse_md_tables.py`
- `rerun_phase6_only.py`
- `run_just_pair_nb.py`
- `run_just_unique_nb.py`
- `run_just_visium.py`
- `test_apis.py`
- `test_b_cell.py`
- `test_db.py`
- `test_network.py`
- `test_opentargets.py`
- `test_subclone.py`
- `test_v1_metrics.py`
