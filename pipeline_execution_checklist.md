# MetabConnectomeDB Pipeline Execution Checklist

This checklist guarantees a 100% reproducible execution of the pipeline from end-to-end. Every script listed below is fully dynamic, pulling variables directly from `pipeline.config.json` via `pan_cancer_config.py`.

> [!NOTE]
> **Data Fetching Automation:** You do **NOT** need to manually run `fetch_reactome_pathways.py`, `fetch_opentargets.py`, or similar fetchers. They are automatically sourced and executed as `STEP 0` within `run_cancer_pipeline.py` to guarantee data integrity before analysis. Furthermore, `serotonin_config.py` contains an auto-healing fail-safe that automatically fetches any missing pathway JSONs on the fly during downstream generation.

> [!CAUTION]
> Do **NOT** modify hardcoded variables within scripts. If parameters need to change (e.g., cell sub-sampling, thresholds), edit `pipeline.config.json` and re-run the pipeline.

> [!NOTE]
> **Jupyter Notebook Execution Policy:** You do **NOT** need to manually run any generated `.ipynb` notebooks. The Python orchestrator scripts (such as `execute_and_export_notebooks.py`, `run_all_cancers.py`, and `execute_pancancer_notebooks.py`) automatically execute all cells headlessly and export the final outputs to `.html` reports. The notebooks are generated in the `scripts/` directory solely so you can optionally open them in Jupyter Lab for interactive data exploration.

> [!NOTE]
> **Scientific Integrity & Reproducibility Policy:** Uses of computational randomness (`np.random`, `random_state=42`) and analytic fallbacks (e.g., analyzing all primary cells if specific subset filtering yields insufficient data, or utilizing alternative baseline ML models) are strictly deployed for mathematical permutation testing, spatial manifold consistency (UMAP), or methodological robustness. They **DO NOT** mock, fabricate, or falsify scientific data. These fallback conditions and algorithms are explicitly documented via callout boxes dynamically injected into the exported `.html` Jupyter notebook reports.

## Phase 0: Core Metabolite Database Curation

*This phase builds the foundational metabolite-target pair catalog that all single-cell scripts query against. It is run once to curate the huge generic databases (KEGG, HMDB, CellChat).*

> [!NOTE]
> Metabolites are mapped to the official HMDB names while genes are mapped to the official HGNC gene symbols.

- `[ ]` **0a. Curate and Merge Raw Databases**
  - **Command:** `bash scripts/merge_simplify_annotate.sh`
  - **Purpose:** The master execution wrapper that ingests raw databases from `input/databases/`, runs `merge_dbs_claude.py` to standardize and deduplicate entries, refines and filters them with `generate_final_outputs.py`, enriches them via `annotate_with_hmdb.py`, and finally merges the annotations into the master CSV outputs via `annotate_with_databases.py`.
  - **Output Location:** `output/human_database_merge_*.csv`

## Phase 1: Database Exploration and Reporting

*This phase provides a general exploratory report of the fully integrated Phase 0 database.*

- `[ ]` **1. Generate Database Reports**
  - **Command:** `python scripts/execute_and_export_notebooks.py`
  - **Purpose:** Executes the exploratory database notebooks (`unique_metab_data_exploration`, `metab_targetPair_analysis`) and converts them into HTML reports.
  - **Output Location:** `output/` (HTML reports and CSV outputs)
    - `output/unique_metab_data_exploration_full_report.html`
    - `output/metab_targetPair_analysis_full_report.html`

## Phase 2: Single-Cell Transcriptome Integration (Per-Cancer)

*This phase uses the curated Phase 1 catalog to query the raw CellxGene `.h5ad` files. (You do NOT need to re-run this heavy processing if your `[cancer]_results/` folders are already generated).*

- `[ ]` **2. CellxGene Integration**
  - **Command:** `python scripts/run_all_cancers.py`
  - **Purpose:** Downloads/streams single-cell data, filters by specific tissue types and cell counts (e.g., 500k as maximum number of cells per cancer defined in `pipeline.config.json`), and outputs the `[cancer]_results` directories containing the sliced `.h5ad` files. It then dynamically executes three core template scripts for each cancer to perform differential expression: `cancer_cellxgene_integration.ipynb`, `primary_vs_metastasis_comparison.ipynb`, and `orphan_metabolic_immune_evasion.ipynb`.
  - **Output Location:** `output/[cancer]_results/` (Contains `.h5ad` files, DE CSVs, and HTML reports)
    - `cancer_cellxgene_integration_[tissue_slug]_[cap].html`
    - `primary_vs_metastasis_[tissue_slug]_[cap].html`
    - `orphan_immune_[tissue_slug]_[cap].html`

## Phase 3: Metrics Extraction & Microenvironment Quantification

*These scripts extract the relevant metrics and targets from the pre-computed Phase 2 `*_results` folders.*

- `[ ]` **3. Extract Single-Cell Dataset Metrics**
  - **Command:** `python scripts/extract_dataset_metrics.py`
  - **Purpose:** Scans the sliced `h5ad` files to extract exactly how many cells were evaluated for the TME and Malignant compartments across primary and metastatic sites.
  - **Output Location:** `output/pan_cancer_meta_results/cell_type_counts_*.csv`
  
- `[ ]` **4. Quantify Immune Evasion Microenvironments**
  - **Command:** `python scripts/quantify_immune_evasion_ccc.py`
  - **Purpose:** Extracts and aggregates the total unique metabolic CCC targets and immune evasion targets from the Phase 2 CSV reports into a comprehensive summary file.
  - **Output Location:** `output/pan_cancer_meta_results/immune_evasion_ccc_quantification_*.csv`

## Phase 4: Pan-Cancer Meta-Analysis & Signature Derivation

*These scripts generate the unified signatures and single-cell subclone scoring.*

- `[ ]` **5. Compute Pan-Cancer Metastatic Signature & CCC Interactome**
  - **Command:** `python scripts/compute_pan_cancer_meta.py` (or through `generate_combined_pan_cancer_notebook.py`)
  - **Purpose:** Intersects the differentially expressed targets across the cancers to identify the strictly conserved signature. **Methodology Note (MAX CANCER - 1 Rule):** If the strict intersection yields 0 genes, the pipeline automatically falls back to utilizing combinations of N-1 cancers to ensure a robust meta-signature is evaluated downstream. 
  - **Fail-Safe Behavior (LIANA):** During Pan-Cancer Conserved CCC Links intersection, if any cancer is missing its `_cellxgene_liana_results.csv`, the notebook automatically invokes `patch_liana_csvs.py` on the fly to regenerate it and prevent pipeline failure.
  - **Output Location:** `output/pan_cancer_meta_results/` (CSVs, upset plot png, network plot png, CCC bar charts, and `pan_cancer_conserved_ccc_links.csv`)
  
- `[ ]` **6. Pre-Metastatic Subclone Resolution**
  - **Command:** `python scripts/generate_predictive_notebook.py`
  - **Purpose:** Parses subset combination signatures and scores primary tumor cells directly to identify pre-metastatic subclones.
  - **Output Location:** Generates `scripts/predictive_signature_biomarker.ipynb` (Executed later in Phase 6)

- `[ ]` **7. Generate Unified Tabular Summaries**
  - **Command:** `python scripts/generate_ai_summary_tables.py`
  - **Purpose:** Aggregates all the extracted outputs into formatted CSVs in `output/ai_summary_tables/`.
  - **Fail-Safe Behavior:** Automatically invokes `fetch_opentargets.py` and `fetch_uniprot_roles.py` if the respective annotation data is missing or incomplete for any newly derived pan-cancer signature genes. This ensures external knowledge graphs stay synchronized without manual intervention.
  - **Output Location:** `output/ai_summary_tables/`

## Phase 5: Dynamic Gene Signature Validation

*This phase is orchestrator-led and validates all dynamically identified (N-1)-cancer combinations generated in Phase 4.*

- `[ ]` **8. Validate Derived Signatures**
  - **Command:** `python scripts/run_validation_phase.py`
  - **Purpose:** The master execution script that automatically identifies all output combination signatures and routes them through the following tests:
    - `massspec_metabolomics_analysis.py`: Verifies signature genes using mass-spectrometry clinical cohorts.
    - `validate_tcga_signature.py`: Runs Cox Proportional Hazard regressions on TCGA survival datasets.
    - `verify_spatial.py`: Applies spatial enrichment scoring and calculates Moran's I on high-resolution Visium slides.
    - `generate_predictive_notebook.py`: Scores primary tumor cells directly to identify left/right skewed pre-metastatic subclones.
    - `generate_ml_prognostic_classifier_notebook.py`: Generates the ML Prognostic Classifier notebook.
    - `create_camp_notebook.py`: Generates the CAMP Pan-Cancer metabolomics integration notebook.
    - `generate_master_regulator_notebook.py`: Generates the Master Regulator TF analysis notebook.
  - **Output Location:** `output/pan_cancer_meta_results/`, `output/tcga_validation/`, `output/spatial_verification/`, and `scripts/` (Validation CSVs, plots, and generated notebooks)

## Phase 6: Advanced Downstream Validation & Export

*This phase pulls advanced API data (Druggability) and exports finalized notebooks using the derived Pan-Cancer signatures.*

- `[ ]` **9. Downstream Analysis & HTML Report Export**
  - **Command:** `python scripts/execute_pancancer_notebooks.py`
  - **Purpose:** This script acts as the master driver for the entire downstream notebook compilation. It executes and exports the spatial, pan-cancer meta, and druggability notebooks into the final HTML reports required by the AI generation phase.
  - **Output Location:** `output/pan_cancer_meta_results/`, `output/druggability/`, and `output/` (HTML reports)
    - `output/pan_cancer_meta_results/pan_cancer_meta_analysis_report.html`
    - `output/pan_cancer_meta_results/predictive_signature_biomarker_[SUFFIX].html`
    - `output/druggability/druggability_axis_[SUFFIX].html`
    - `output/visium_spatial_validation_report.html`
    - `output/deepdive_conserved_metabGeneSig/deepdive_conserved_metabGeneSig_report.html`
    - `output/serotonin_axis_spatial_mapping_report.html`
    - `output/ovarian_serotonin_immune_evasion_report.html`
    - `output/oxygen_tension_analysis_report.html` (Note: Requires `scripts/compute_metabolic_switching.py` to be present; do not archive this dependency)
    - `output/mitf_regulon_expansion_report.html`
    - `output/[cancer]_ml_prognostic_classifier_report.html`
    - `output/pancancer_ml_prognostic_classifier_report.html`
    - `output/camp_pancancer_integration_report.html`
    - `output/master_regulator_analysis_report.html`

## Phase 7: Dynamic AI Insights Report Generation

*This final phase reads the pre-computed CSVs and validated outputs from all prior phases to build a summarized Markdown document with AI-generated interpretations per section.*

- `[ ]` **10. Build Dynamic AI Summary Document**
  - **Command:** `python scripts/tmp_build_md.py --phase all --reset`
  - **Purpose:** Aggregates and **summarizes** (not raw-dumps) the pipeline outputs into `AI_summary_and_insights.md`. For each phase section, it:
    - Reads the actual CSV/output files directly (no HTML scraping for validated data)
    - Computes summary metrics (e.g., Top 5 signatures by significance, median C-Index, detection rates)
    - Includes methodology notes explaining how each metric was derived
    - Uploads reference PDFs from `papers/` directory and phase-specific Jupyter HTML reports to Gemini via the File API for deep biological cross-referencing
    - Sends the condensed summary + uploaded files to Gemini API for AI interpretation
    - Post-processes every AI response through a **PMID Verification Pipeline** that: (1) extracts all cited PMIDs, (2) verifies them via NCBI E-Fetch API, (3) runs a secondary semantic verification LLM call to confirm the paper's abstract supports the claim, and (4) formats verified PMIDs as clickable `[PMID:X - Title](url)` links while flagging hallucinated citations
    - Includes 503/UNAVAILABLE retry logic with exponential backoff (up to 5 retries)
  - **Flags:**
    - `--phase all` — Builds all phases (1–6) sequentially
    - `--phase 6` — Builds only Phase 6 (useful for testing)
    - `--reset` — Deletes the existing `AI_summary_and_insights.md` before rebuilding
  - **Output Location:** `output/AI_summary_and_insights.md`

> [!NOTE]
> **Phase 7 Data Flow:** `tmp_build_md.py` reads directly from:
> - `papers/*.pdf` (Reference literature uploaded via Gemini File API for cross-referencing)
> - `output/**/*.html` (Phase-specific Jupyter notebook HTML reports uploaded via Gemini File API)
> - `output/ai_summary_tables/*.csv` (Phases 1–5 summaries from `generate_ai_summary_tables.py`)
> - `output/tcga_validation/*/true_signature_metrics.csv` (Phase 6.1 TCGA)
> - `output/pan_cancer_meta_results/pre_metastatic_subclone_summary*.csv` (Phase 6.2 Subclones)
> - `output/spatial_verification/*/visium_spatial_clustering_summary.csv` (Phase 6.3 Spatial)
> - `output/massspec_metabolomics/*/per_gene_metabolite_profile.csv` (Phase 6.4 MassSpec)
> - `output/ml_prognostic_results/tcga/*/ml_metrics.csv` (Phase 6.6 ML Classifiers)

---

## Recent Script Modifications (Changelog)

| Date | Script | Change |
| --- | --- | --- |
| 2026-06-20 | `scripts/tmp_build_md.py` | **Phase 6 Summarization Refactor:** Replaced raw table dumping with aggregated Top-5/Top-3 summary tables for TCGA (6.1), Subclones (6.2), Spatial (6.3), MassSpec (6.4), and ML (6.6) sections |
| 2026-06-20 | `scripts/tmp_build_md.py` | **Methodology Notes:** Added `**Methodology:**` context blocks to all Phase 6 subsections explaining metric derivation for AI summarizer guidance |
| 2026-06-20 | `scripts/tmp_build_md.py` | **API Resilience:** Added `503 UNAVAILABLE` to retry loop with exponential backoff (up to 5 retries) |
| 2026-06-20 | `scripts/tmp_build_md.py` | **Phase 3 Disease Counts:** Filtered out zero-count rows, stripped string artifacts (`...`), formatted as clean table |
| 2026-06-20 | `scripts/generate_ai_summary_tables.py` | **pO2 Fix:** Added missing `normalize_cancer_name` import; fixed Tumor/Normal pO2 lookup logic via `CANCER_PO2_CSV_MAPPING` |
| 2026-06-20 | `scripts/generate_ml_prognostic_classifier_notebook.py` | **CSV Export:** Notebook already saves `ml_metrics.csv` with schema: Signature, Train_Size, Test_Size, Optimal_Penalizer, CV_C_Index, Test_C_Index |
| 2026-06-22 | `pipeline.config.json`, `cancer_cellxgene_integration.ipynb`, `run_cancer_pipeline.py` | **Dynamic KEGG & Plot Aesthetics Fix:** Replaced hardcoded immune-metabolic gene targets with dynamic config JSON mapping. Fixed `'ovary'` key mismatch. Implemented `sc.pl.stacked_violin` and mathematically scaled bubble plots/heatmaps to completely eliminate axis label overlap. Fixed redundant `OUTPUT_DIR` script injections. Added TCA Cycle (Fumarate) pathways. |
| 2026-06-23 | `cancer_cellxgene_integration.ipynb`, `primary_vs_metastasis_comparison.ipynb`, `orphan_metabolic_immune_evasion.ipynb` | **Zero-Tolerance Provenance & Visual Exports Patch:** Removed hardcoded KEGG fallback in integration notebook to enforce hard crash. Split Section 9 squished dotplots into dynamically-sized, pathway-specific plots. Added explicit standalone PDF and CSV exports for previously inline-only visuals (Static Volcano, Niche UpSet, Orphan Multi-Panel Bar). |
| 2026-06-23 | `generate_combined_pan_cancer_notebook.py` | **Phase 4 LIANA Integration:** Hooked up `immune_evasion_ccc_quantification` and `_cellxgene_liana_results.csv`. Implemented automatic fail-safe to execute `patch_liana_csvs.py` if LIANA results are missing during Pan-Cancer intersection. |
| 2026-06-23 | `generate_final_outputs.py` | **HGNC Target Canonicalization Patch:** Added canonical resolution mapping via `input/hgnc_approved_genes.json` to process `Target_original` into official `Target` names. This completely eliminates obsolete Uniprot TrEMBL fragmentation downstream. |
| 2026-06-24 | `run_cancer_pipeline.py`, `serotonin_config.py` | **Auto-Fetching & Self-Healing:** Configured the pipeline orchestrator to automatically discover and source all `fetch_*.py` scripts as STEP 0 with a hard-fail condition. Implemented an auto-healing subroutine in `serotonin_config.py` to trigger Reactome/API fetches dynamically if config JSONs are missing, completely removing the need for manual API script executions. |
| 2026-06-25 | `tmp_build_md.py`, `compute_metabolic_switching.py`, `master_regulator_analysis.ipynb` | **Downstream Extraction Integrity Patch:** Rewrote the `tmp_build_md.py` scraper to explicitly extract statistical metrics (P-values, HRs, target counts) from actual CSV outputs rather than falsely reporting simple file existence. Fixed the computation script to properly save the oxygen tension results, and reverted the Master Regulator notebook back to its dynamically un-hardcoded state to ensure complete data provenance. |
| 2026-06-25 | `tmp_build_md.py` | **Deep Research & PMID Verification Patch:** Integrated Gemini File API to upload `papers/*.pdf` reference literature and phase-specific Jupyter HTML reports for deep biological cross-referencing. Implemented a 3-stage PMID Verification Pipeline: (1) NCBI E-Fetch existence check, (2) semantic verification via secondary LLM call comparing the claim against the paper's abstract, (3) formatting verified PMIDs as clickable markdown links with titles. Hallucinated or unsupported citations are automatically flagged and removed. |
